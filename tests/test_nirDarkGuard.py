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


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
