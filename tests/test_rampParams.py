"""Tests for ``PfsIsrTask.rampParams`` — the CDS-vs-UTR and read-range dispatch.

Two behaviours are pinned down here:

1. The read range follows the *resolved* reduction mode, not the observation
   type alone. A science ramp reduced UTR drops the shutter-closed reads
   (``1:-4``); the same ramp reduced CDS wants those endpoint reads and so uses
   the whole ramp (``0:-1``).
2. The resolved range is checked against the ramp's length. A range that does
   not fit — a 2-read ramp asked for ``lastRead=-4`` — is a fatal error with an
   explanatory message, not a bare ``IndexError`` from deep inside the read
   accessors.
"""
import unittest

import lsst.utils.tests

from testUtils import HAS_DRP_STELLA, requireDrpStella

if HAS_DRP_STELLA:
    # isrTask imports pfs.drp.stella.crosstalk.
    from lsst.obs.pfs import isrTask as pfsIsrTask


class _FakeObsInfo:
    def __init__(self, observationType):
        self.observation_type = observationType


class _FakeRaw:
    """Minimal ``PfsRaw`` surface ``rampParams`` needs: the observation type and
    the ramp length (via ``getNumReads`` and the ``positiveIndex`` built on it).
    """

    def __init__(self, observationType, numReads):
        self.obsInfo = _FakeObsInfo(observationType)
        self._numReads = numReads

    def getNumReads(self):
        return self._numReads

    def positiveIndex(self, n):
        # Same semantics as PfsRaw.positiveIndex: out-of-range raises IndexError.
        return range(0, self._numReads)[n]


def _makeIsrTask(**h4):
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doFlat = False
    config.doDark = True
    config.doDefect = True
    config.doSaturationInterpolation = False
    config.h4.doIPC = False
    config.h4.doLinearize = False
    config.h4.doCR = False
    for name, value in h4.items():
        setattr(config.h4, name, value)
    config.validate()
    return pfsIsrTask.PfsIsrTask(config=config)


@requireDrpStella
class RampParamsDispatchTestCase(lsst.utils.tests.TestCase):
    """Per-obstype defaults, and the config overrides layered on top."""

    def testScienceDefaultsToTrimmedUTR(self):
        task = _makeIsrTask()
        self.assertEqual(task.rampParams(_FakeRaw("science", 90)), (False, 1, -4))

    def testScienceCDSUsesWholeRamp(self):
        # CDS differences the endpoint frames, so the shutter-closed reads the
        # UTR trim drops are exactly the ones it needs: no trim.
        task = _makeIsrTask(quickCDS=True)
        self.assertEqual(task.rampParams(_FakeRaw("science", 90)), (True, 0, -1))

    def testShortScienceRampReducesUnderCDS(self):
        # The ticket's motivating case: a 2-read simulation ramp.
        task = _makeIsrTask(quickCDS=True)
        self.assertEqual(task.rampParams(_FakeRaw("science", 2)), (True, 0, -1))

    def testArcsAndFlatsUseCDSOverWholeRamp(self):
        task = _makeIsrTask()
        for obsType in ("comparison", "flat"):
            self.assertEqual(task.rampParams(_FakeRaw(obsType, 90)), (True, 0, -1))

    def testDarkAndUnknownUseFullRampUTR(self):
        task = _makeIsrTask()
        for obsType in ("dark", "bias", "", None):
            self.assertEqual(task.rampParams(_FakeRaw(obsType, 90)), (False, 0, -1))

    def testObsTypeIsCaseInsensitive(self):
        task = _makeIsrTask()
        self.assertEqual(task.rampParams(_FakeRaw("SCIENCE", 90)), (False, 1, -4))
        self.assertEqual(task.rampParams(_FakeRaw("Flat", 90)), (True, 0, -1))

    def testUTRForcedOnArc(self):
        # quickCDS=False on an arc: UTR, and the arc's whole-ramp range stands.
        task = _makeIsrTask(quickCDS=False)
        self.assertEqual(task.rampParams(_FakeRaw("comparison", 90)), (False, 0, -1))

    def testExplicitReadsOverrideScienceUTRDefault(self):
        task = _makeIsrTask(firstRead=2, lastRead=-2)
        self.assertEqual(task.rampParams(_FakeRaw("science", 90)), (False, 2, -2))

    def testExplicitReadsOverrideCDSDefault(self):
        task = _makeIsrTask(quickCDS=True, firstRead=3, lastRead=-5)
        self.assertEqual(task.rampParams(_FakeRaw("science", 90)), (True, 3, -5))

    def testSingleExplicitReadOverridesOneEnd(self):
        # Only lastRead set: firstRead keeps its dispatched science-UTR value.
        task = _makeIsrTask(lastRead=-1)
        self.assertEqual(task.rampParams(_FakeRaw("science", 90)), (False, 1, -1))


@requireDrpStella
class RampParamsValidationTestCase(lsst.utils.tests.TestCase):
    """The resolved range must fit the ramp, with at least one delta in it."""

    def testShortScienceRampUTRRaises(self):
        task = _makeIsrTask()
        with self.assertRaises(RuntimeError):
            task.rampParams(_FakeRaw("science", 2))

    def testScienceUTRBoundary(self):
        # 1:-4 needs r1 = n - 4 > r0 = 1, i.e. n >= 6.
        task = _makeIsrTask()
        self.assertEqual(task.rampParams(_FakeRaw("science", 6)), (False, 1, -4))
        with self.assertRaises(RuntimeError):
            task.rampParams(_FakeRaw("science", 5))

    def testSingleReadRampRaises(self):
        # Even the whole-ramp default has no delta in a 1-read ramp.
        task = _makeIsrTask(quickCDS=True)
        with self.assertRaises(RuntimeError):
            task.rampParams(_FakeRaw("science", 1))

    def testExplicitOutOfRangeReadsRaise(self):
        for kwargs in ({"lastRead": -20}, {"firstRead": 40}, {"firstRead": -30}):
            task = _makeIsrTask(quickCDS=True, **kwargs)
            with self.assertRaises(RuntimeError):
                task.rampParams(_FakeRaw("science", 10))

    def testInvertedExplicitRangeRaises(self):
        task = _makeIsrTask(quickCDS=True, firstRead=8, lastRead=3)
        with self.assertRaises(RuntimeError):
            task.rampParams(_FakeRaw("science", 10))

    def testErrorMessageIsInformative(self):
        task = _makeIsrTask()
        with self.assertRaises(RuntimeError) as cm:
            task.rampParams(_FakeRaw("science", 2))
        message = str(cm.exception)
        for expected in ("science", "2", "-4", "quickCDS"):
            self.assertIn(expected, message)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
