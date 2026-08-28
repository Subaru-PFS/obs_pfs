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
        self.assertAlmostEqual(window.off, 40.0, places=6)
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
        self.assertAlmostEqual(window.off, 40.0, places=6)
        self.assertGreater(window.offUpper, window.off)
        self.assertAlmostEqual(window.offUpper, 40.0 + 0.16*10.0, places=6)

    def testLampStillOnAtRampEnd(self):
        """No trailing dark read is reported as such, not invented."""
        cube = FakeCube([0.0, 100.0, 100.0, 100.0])
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertTrue(window.lampStillOn)
        self.assertEqual(window.trailingDark, 0)

    def testLeadingAndTrailingDarkCounts(self):
        cube = FakeCube([0.0, 0.0, 100.0, 0.0, 0.0])
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertGreaterEqual(window.leadingDark, 1)
        self.assertGreaterEqual(window.trailingDark, 1)

    def testTooFaintReturnsNone(self):
        """Arcs are often too faint to locate an edge; say so rather than guess."""
        cube = FakeCube([0.0, 1.0, 1.0, 0.0])
        self.assertIsNone(measureIllumination(cube, minPlateau=5.0))

    def testDurationAgainstCommanded(self):
        cube = FakeCube([0.0, 100.0, 100.0, 100.0, 0.0], exptime=30.0)
        window = measureIllumination(cube, minPlateau=1.0)
        self.assertAlmostEqual(window.duration, 30.0, places=6)

    def testTrussLampsAreReported(self):
        """getLamps knows nothing of the truss lamps, which are voltages."""
        md = {"DATA-TYP": "test", "W_TFF2VV": 40.0}
        self.assertEqual(lampsOn(md), {"truss2"})

    def testQuartzStillReported(self):
        md = {"DATA-TYP": "test", "W_AITQTH": True}
        self.assertIn("Quartz", lampsOn(md))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
