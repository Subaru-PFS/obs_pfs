"""Tests for ``PfsIsrTask.checkNirDark`` — the guard that rejects a NIR dark
that does not match the ramp about to be dark-subtracted.

Two independent hard requirements:

1. the dark's IRP ratio (``W_H4IRPN``) equals the ramp's, and
2. the dark has enough reads to cover the ramp's read range.

Either mismatch must raise, rather than surfacing as a cryptic
``KeyError: 'Extension IMAGE_<n> not found'`` deep in ``subtractDarkCube``.
"""
import unittest

import lsst.utils.tests
from lsst.obs.pfs import isrTask as pfsIsrTask


class _FakeRaw:
    """Minimal ``PfsRaw`` surface ``checkNirDark`` needs: the ramp IRP ratio."""

    def __init__(self, irpN):
        self.irpN = irpN


class _FakeDark:
    """Minimal ``ImageCube`` surface: ``metadata`` (with ``W_H4IRPN``) and
    ``getNumReads()``."""

    def __init__(self, irpN, nreads):
        self.metadata = {"W_H4IRPN": irpN}
        self._nreads = nreads

    def getNumReads(self):
        return self._nreads


def _makeIsrTask():
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doFlat = False
    config.doDark = True
    config.doDefect = True
    config.doSaturationInterpolation = False
    config.h4.quickCDS = False
    config.h4.doIPC = False
    config.h4.doWriteRawCube = True
    config.h4.doLinearize = False
    config.h4.doCR = False
    config.validate()
    return pfsIsrTask.PfsIsrTask(config=config)


class CheckNirDarkTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.task = _makeIsrTask()

    def testRaisesOnIrpRatioMismatch(self):
        # irp4 ramp, irp1 dark: different per-read cadence -> must reject.
        raw = _FakeRaw(irpN=4)
        dark = _FakeDark(irpN=1, nreads=200)
        with self.assertRaises(RuntimeError):
            self.task.checkNirDark(raw, dark, nReadsNeeded=90)

    def testRaisesWhenDarkTooShort(self):
        # Ratios match, but the dark has fewer reads than the ramp needs.
        raw = _FakeRaw(irpN=1)
        dark = _FakeDark(irpN=1, nreads=50)
        with self.assertRaises(RuntimeError):
            self.task.checkNirDark(raw, dark, nReadsNeeded=90)

    def testPassesWhenRatioMatchesAndEnoughReads(self):
        # Matching ratio and enough reads: no error.
        raw = _FakeRaw(irpN=1)
        dark = _FakeDark(irpN=1, nreads=100)
        self.task.checkNirDark(raw, dark, nReadsNeeded=90)

    def testExactReadCountIsEnough(self):
        # nReadsNeeded reads cover indices 0..nReadsNeeded-1 exactly.
        raw = _FakeRaw(irpN=1)
        dark = _FakeDark(irpN=1, nreads=90)
        self.task.checkNirDark(raw, dark, nReadsNeeded=90)

    def testOldDarkWithoutRatioAcceptedForIrp1(self):
        # An old dark predating the W_H4IRPN card (ratio absent) is irp1;
        # an irp1 ramp must accept it (backward compatibility).
        raw = _FakeRaw(irpN=1)
        dark = _FakeDark(irpN=None, nreads=88)
        self.task.checkNirDark(raw, dark, nReadsNeeded=86)

    def testOldDarkRejectedForIrp4Ramp(self):
        # The same old (implicitly irp1) dark must NOT serve an irp4 ramp.
        raw = _FakeRaw(irpN=4)
        dark = _FakeDark(irpN=None, nreads=200)
        with self.assertRaises(RuntimeError):
            self.task.checkNirDark(raw, dark, nReadsNeeded=90)


class _GainDark:
    """Minimal ``ImageCube`` surface for the GAIN guard: just ``metadata``."""

    def __init__(self, gain=None):
        self.metadata = {} if gain is None else {"GAIN": gain}


class DarkGainGuardTestCase(lsst.utils.tests.TestCase):
    """``_darkGain`` rejects a non-physical nirDark GAIN (e.g. the 9999 raw
    placeholder) rather than silently mis-scaling the dark subtraction.
    """

    def setUp(self):
        self.task = _makeIsrTask()

    def testPhysicalGainReturned(self):
        self.assertAlmostEqual(self.task._darkGain(_GainDark(3.097)), 3.097)

    def testMissingGainIsUnity(self):
        # No GAIN header -> 1.0 (dark already in ADU); allowed, no division.
        self.assertEqual(self.task._darkGain(_GainDark()), 1.0)

    def testPlaceholder9999Raises(self):
        with self.assertRaises(RuntimeError):
            self.task._darkGain(_GainDark(9999.0))

    def testNonPhysicalGainRaises(self):
        for bad in (0.0, -3.0, 1.0e6):
            with self.assertRaises(RuntimeError):
                self.task._darkGain(_GainDark(bad))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
