"""``PfsIsrTask.getSimpleDiffIrp`` must replace each ASIC channel with its
per-column median (broadcast down the channel's rows).

The original implementation tiled the per-column vector by the full image
height instead of the channel height, so the assignment raised a broadcast
``ValueError`` for any real (nchan > 1) image — i.e. the ``IRPfilter=-1`` path
never worked.
"""
import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs import isrTask as pfsIsrTask


class _FakeRaw:
    """Minimal surface ``getSimpleDiffIrp`` needs: channel count and IRP ratio."""

    def __init__(self, nchan, irpN=1):
        self.nchan = nchan
        self.irpN = irpN


class _FakeBadRefPixels:
    """Minimal ``badRefPixels`` calib surface: a list of bad reference rows."""

    def __init__(self, pixels):
        self.pixels = np.asarray(pixels)


def _makeIsrTask():
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doDark = False
    config.doFlat = False
    config.doDefect = False
    config.doSaturationInterpolation = False
    config.h4.quickCDS = False
    config.h4.doIPC = False
    config.h4.doWriteRawCube = False
    config.h4.doLinearize = False
    config.h4.doCR = False
    config.validate()
    return pfsIsrTask.PfsIsrTask(config=config)


class GetSimpleDiffIrpTestCase(lsst.utils.tests.TestCase):
    def testPerChannelPerColumnMedian(self):
        task = _makeIsrTask()
        raw = _FakeRaw(nchan=2)
        # 4 rows = two 2-row channels; 3 columns.
        img = np.array([[1., 2., 3.],
                        [3., 4., 5.],     # channel 0: rows 0-1
                        [10., 20., 30.],
                        [40., 60., 80.]])  # channel 1: rows 2-3
        out = task.getSimpleDiffIrp(raw, img)
        # per-channel, per-column median, broadcast down the channel's rows
        expected = np.array([[2., 3., 4.],
                             [2., 3., 4.],
                             [25., 40., 55.],
                             [25., 40., 55.]])
        self.assertEqual(out.shape, img.shape)
        np.testing.assert_array_equal(out, expected)

    def testExcludesBadRefRows(self):
        # Bad reference rows must be left out of the per-column median so they
        # can't bias it. Channel 0 rows 2,3 are outliers and flagged bad; the
        # median must come from the good rows (0,1) only.
        task = _makeIsrTask()
        task._badRefPixels = _FakeBadRefPixels([2, 3])
        raw = _FakeRaw(nchan=2)  # 8 rows -> two 4-row channels
        img = np.array([[10.], [10.], [999.], [999.],   # channel 0: rows 2,3 bad
                        [20.], [20.], [20.], [20.]])     # channel 1: clean
        out = task.getSimpleDiffIrp(raw, img)
        # Good-rows median is 10, not the all-rows median 504.5.
        np.testing.assert_allclose(out[0:4, 0], 10.0)
        np.testing.assert_allclose(out[4:8, 0], 20.0)

    def testDoesNotMaskForIrp4(self):
        # badRefPixels is an IRP1-frame calib, so masking is disabled for
        # irpN > 1: the median runs over all rows (bad rows included).
        task = _makeIsrTask()
        task._badRefPixels = _FakeBadRefPixels([2, 3])
        raw = _FakeRaw(nchan=2, irpN=4)
        img = np.array([[10.], [10.], [999.], [999.],
                        [20.], [20.], [20.], [20.]])
        out = task.getSimpleDiffIrp(raw, img)
        # No masking -> all-rows median of channel 0 = 504.5.
        np.testing.assert_allclose(out[0:4, 0], 504.5)
        np.testing.assert_allclose(out[4:8, 0], 20.0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
