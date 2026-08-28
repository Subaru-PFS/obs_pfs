"""Tests for measuring a ramp's illumination window from its flux.

The interesting cases are the ones that bit in practice: a read only partly
illuminated must place the edge inside that read rather than at a boundary,
and a ramp whose lamp is still on at the end must report that rather than
inventing a trailing dark read.
"""
import sys
import unittest

import numpy as np

import lsst.utils.tests
from lsst.daf.base import PropertyList

from lsst.obs.pfs.h4utils.rampTiming import measureIllumination, lampsOn


class FakeCube:
    """Minimal stand-in for an ImageCube: cumulative planes plus metadata."""

    def __init__(self, rates, frameTime=10.0, exptime=30.0, read0=0,
                 shape=(4, 4)):
        self.metadata = PropertyList()
        self.metadata["W_H4FRMT"] = frameTime
        self.metadata["H4READ0"] = read0
        self.metadata["EXPTIME"] = exptime
        self.metadata["DARKTIME"] = (len(rates) + 1)*frameTime
        self.metadata["MJD-STR"] = 60000.0
        self._planes = np.cumsum(np.asarray(rates, float)*frameTime)
        self._shape = shape

    def getNumReads(self):
        return len(self._planes)

    def getReadArray(self, index):
        return np.full(self._shape, self._planes[index], dtype=np.float32)


class RampTimingTestCase(lsst.utils.tests.TestCase):
    def testFullyLitReadsGiveBoundaryEdges(self):
        """A lamp switching exactly at read boundaries lands on them."""
        cube = FakeCube([0.0, 100.0, 100.0, 0.0])
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertAlmostEqual(window.on, 20.0, places=6)
        self.assertAlmostEqual(window.offFlux, 40.0, places=6)
        self.assertFalse(window.lampStillOn)

    def testPartialReadPlacesEdgeInsideIt(self):
        """A 20%-lit read puts lamp-on 80% of the way through that read.

        Treating such a read as dark -- which a 50%-of-plateau threshold does
        -- would truncate the answer to the read boundary and misplace the
        edge by most of a read.
        """
        cube = FakeCube([0.0, 20.0, 100.0, 100.0, 0.0])
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertAlmostEqual(window.on, 30.0 - 0.2*10.0, places=6)
        self.assertAlmostEqual(window.litFraction(1), 0.2, places=6)

    def testTrailingEdgeIsBracketed(self):
        """Prompt persistence in the straddling read makes off a lower bound."""
        cube = FakeCube([0.0, 100.0, 100.0, 16.0])
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertAlmostEqual(window.offFlux, 40.0, places=6)
        self.assertGreater(window.offUpper, window.offFlux)
        self.assertAlmostEqual(window.offUpper, 40.0 + 0.16*10.0, places=6)

    def testLampStillOnAtRampEnd(self):
        """No trailing dark read is reported as such, not invented."""
        cube = FakeCube([0.0, 100.0, 100.0, 100.0])
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertTrue(window.lampStillOn)
        self.assertEqual(window.trailingDark, 0)

    def testLeadingAndTrailingDarkCounts(self):
        cube = FakeCube([0.0, 0.0, 100.0, 100.0, 0.0, 0.0])
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertGreaterEqual(window.leadingDark, 1)
        self.assertGreaterEqual(window.trailingDark, 1)

    def testAnomalousFirstReadIsNotAPlateau(self):
        """A single hot read must not be mistaken for the illuminated level.

        Read 0 carries reset and settling behaviour and on a faint flat can
        exceed every illuminated read -- as on 146285/n2, where read 0 ran at
        7.0 e-/s against ~3.1 for the lit reads. Taking the peak rate as the
        plateau then reports the lamp coming on in read 0, which it did not.
        """
        cube = FakeCube([7.0, 3.0, 3.1, 2.3, 0.4])
        # The single hot read cannot hold a plateau across two reads, so no
        # window is reported rather than a bogus one. Read 0 is no longer
        # skipped -- the lamp genuinely does arrive in the first interval on
        # real ramps -- so this guard is what protects against it.
        self.assertIsNone(measureIllumination(cube, minPlateau=1.0))

    def testPlateauHeldAcrossTwoReadsIsAccepted(self):
        """The guard must not reject a genuine short plateau."""
        cube = FakeCube([0.0, 100.0, 100.0, 0.0])
        self.assertIsNotNone(measureIllumination(cube, minPlateau=1.0))

    def testTooFaintReturnsNone(self):
        """Arcs are often too faint to locate an edge; say so rather than guess."""
        cube = FakeCube([0.0, 1.0, 1.0, 0.0])
        self.assertIsNone(measureIllumination(cube, minPlateau=5.0))

    def testDurationAgainstCommanded(self):
        cube = FakeCube([0.0, 100.0, 100.0, 100.0, 0.0], exptime=30.0)
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertAlmostEqual(window.duration, 30.0, places=6)
        self.assertAlmostEqual(window.off, window.on + 30.0, places=6)

    def testScanEdgeRecoversMidReadTransition(self):
        """A lamp switching mid-read leaves a ramp across the scan direction.

        Built as a synthetic frame whose illuminated fraction rises linearly
        with image column, which is what the row pointer produces: columns
        scanned before the switch see nothing, those after see progressively
        more.
        """
        from lsst.obs.pfs.h4utils.rampTiming import measureEdgeByScan

        ncol, nrow, frame = 4096, 64, 10.0
        onColumn = 1600.0
        ramp = np.clip((np.arange(ncol) - onColumn)/ncol, 0.0, 1.0)
        full = np.tile(np.full(ncol, 100.0), (nrow, 1))*frame
        partial = np.tile(ramp*100.0, (nrow, 1))*frame

        class ScanCube:
            def __init__(self):
                self.metadata = PropertyList()
                self.metadata["W_H4FRMT"] = frame
                self.metadata["H4READ0"] = 0
                self.metadata["EXPTIME"] = 30.0
                self.metadata["DARKTIME"] = 4*frame
                self.metadata["MJD-STR"] = 60000.0
                self._p = [partial, partial + full, partial + 2*full]

            def getNumReads(self):
                return len(self._p)

            def getReadArray(self, i):
                return self._p[i]

        result = measureEdgeByScan(ScanCube(), litPercentile=50.0)
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result["onColumn"], onColumn, delta=60)
        self.assertAlmostEqual(result["slope"]/result["expectedSlope"], 1.0,
                               delta=0.1)

    def testTrussLampsAreReported(self):
        """getLamps knows nothing of the truss lamps, which are voltages."""
        md = {"DATA-TYP": "test", "W_TFF2VV": 40.0}
        self.assertEqual(lampsOn(md), {"truss2"})

    def testQuartzStillReported(self):
        md = {"DATA-TYP": "test", "W_AITQTH": True}
        self.assertIn("Quartz", lampsOn(md))


class ReadSelectionTestCase(lsst.utils.tests.TestCase):
    """Read selection strings, which are the part users type by hand."""

    def testSelections(self):
        from lsst.obs.pfs.h4utils.displayReads import parseReads
        self.assertEqual(parseReads("all", 5), [0, 1, 2, 3, 4])
        self.assertEqual(parseReads(None, 5), [0, 1, 2, 3, 4])
        self.assertEqual(parseReads("0-2", 5), [0, 1, 2])
        self.assertEqual(parseReads("1,3", 5), [1, 3])
        self.assertEqual(parseReads("2:", 5), [2, 3, 4])
        self.assertEqual(parseReads("-1", 5), [4])
        self.assertEqual(parseReads("-2:", 5), [3, 4])

    def testOutOfRangeIsClipped(self):
        from lsst.obs.pfs.h4utils.displayReads import parseReads
        self.assertEqual(parseReads("3-99", 5), [3, 4])
        self.assertEqual(parseReads("99", 5), [])

    def testDuplicatesCollapse(self):
        from lsst.obs.pfs.h4utils.displayReads import parseReads
        self.assertEqual(parseReads("1,1,0-1", 5), [1, 0])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
