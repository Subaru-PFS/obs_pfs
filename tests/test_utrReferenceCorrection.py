"""``PfsIsrTask.makeUTRcumulative`` must reference-correct each read the same
way ``makeCDS`` does: subtract the IRP reference when ``h4.useIRP`` is set and
interleaved reference pixels are present, but otherwise fall back to the
Teledyne refPixel4 **border** reference (``borderCorrect``).

The bug: the UTR loop subtracted the IRP reference unconditionally, ignoring
``useIRP``. So a ``useIRP=False`` ramp (e.g. a no-IRP dark build) silently got
IRP-corrected anyway -- identical to ``IRPfilter``'s default -- instead of the
border correction it asked for.
"""
import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs import isrTask as pfsIsrTask


class _FakeRaw:
    """Minimal ``PfsRaw`` surface ``makeUTRcumulative`` needs."""

    def __init__(self, irpN, nreads):
        self.irpN = irpN
        self._nreads = nreads

    def positiveIndex(self, i):
        return i if i >= 0 else self._nreads + i


class _RecordingTask(pfsIsrTask.PfsIsrTask):
    """Drive the real ``makeUTRcumulative`` loop but replace the leaf I/O and
    the two reference primitives with deterministic, call-recording stand-ins.

    - ``makeRawDataArray``  -> a 4x4 plane whose value is ``10 * readNum``
    - ``makeRawIrpArray``   -> a 4x4 plane whose value is ``readNum``
    - ``getFinalDiffIrp``   -> constant 3.0 (records that IRP was used)
    - ``borderCorrect``     -> ``image * 2`` (records that border was used)
    """

    SHAPE = (4, 4)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.irpCalls = 0
        self.borderCalls = 0

    def makeRawDataArray(self, pfsRaw, readNum, fromArray=None):
        return np.full(self.SHAPE, 10.0 * readNum, dtype="f4")

    def makeRawIrpArray(self, pfsRaw, readNum, forceIrp1=True, fromArray=None):
        return np.full(self.SHAPE, float(readNum), dtype="f4")

    def applyIRPcrosstalk(self, pfsRaw, irp, data):
        pass

    def getFinalDiffIrp(self, pfsRaw, rawDiffIrp, useFft=True):
        self.irpCalls += 1
        return np.full_like(rawDiffIrp, 3.0)

    def borderCorrect(self, pfsRaw, image, colWindow=4, doRows=True, doCols=True):
        self.borderCalls += 1
        return image * 2.0


def _makeTask(useIRP):
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doFlat = False
    config.doDark = False
    config.doDefect = False
    config.doSaturationInterpolation = False
    config.h4.quickCDS = False
    config.h4.doIPC = False
    config.h4.doWriteRawCube = False
    config.h4.doLinearize = False
    config.h4.doCR = False
    config.h4.useIRP = useIRP
    config.validate()
    return _RecordingTask(config=config)


class MakeUTRcumulativeReferenceTestCase(lsst.utils.tests.TestCase):
    def testUsesIrpWhenUseIrpTrue(self):
        task = _makeTask(useIRP=True)
        raw = _FakeRaw(irpN=1, nreads=3)
        stack = task.makeUTRcumulative(raw, r0=0, r1=2)
        # IRP path: ddata = (data1 - data0) - getFinalDiffIrp = (10r - 0) - 3
        self.assertGreater(task.irpCalls, 0)
        self.assertEqual(task.borderCalls, 0)
        np.testing.assert_allclose(stack[0], 10.0 - 3.0)  # read 1
        np.testing.assert_allclose(stack[1], 20.0 - 3.0)  # read 2

    def testUsesBorderWhenUseIrpFalse(self):
        task = _makeTask(useIRP=False)
        raw = _FakeRaw(irpN=1, nreads=3)
        stack = task.makeUTRcumulative(raw, r0=0, r1=2)
        # Border path: ddata = borderCorrect(data1) - borderCorrect(data0)
        #            = 2*(10r) - 2*0 = 20r ; IRP must NOT be used.
        self.assertEqual(task.irpCalls, 0)
        self.assertGreater(task.borderCalls, 0)
        np.testing.assert_allclose(stack[0], 20.0)  # read 1: 2*10
        np.testing.assert_allclose(stack[1], 40.0)  # read 2: 2*20

    def testUsesBorderWhenNoInterleavedReferencePixels(self):
        # useIRP requested but the ramp has no IRP planes (irpN == 0):
        # must also fall back to border, like makeCDS.
        task = _makeTask(useIRP=True)
        raw = _FakeRaw(irpN=0, nreads=3)
        stack = task.makeUTRcumulative(raw, r0=0, r1=2)
        self.assertEqual(task.irpCalls, 0)
        self.assertGreater(task.borderCalls, 0)
        np.testing.assert_allclose(stack[0], 20.0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
