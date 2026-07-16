"""Per-pixel gate for ASIC-glitch correction: a pixel is *correctable*
(kept rather than masked) only when it is in a correction-scoped channel (all
channels except the known-bad ones, e.g. n3 ch24), not already flagged bad, and
every glitch is a clean two-delta event -- an above-threshold outlier with a
roughly-equal, opposite return. Lone persistent spikes (CR-like) and entangled
multi-outlier pixels stay masked; small sub-threshold strays are ignored. Kept
glitches are flagged in the published (non-BAD) ASIC_GLITCH plane.
"""
from __future__ import annotations

import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs import isrTask as pfsIsrTask


def _makeTask(irpFilter=None):
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doFlat = False
    config.doDark = False
    config.doDefect = False
    config.doSaturationInterpolation = False
    config.h4.quickCDS = False
    config.h4.doIPC = False
    config.h4.doWriteRawCube = False
    config.h4.doLinearize = False
    if irpFilter is not None:
        config.h4.IRPfilter = irpFilter
    config.validate()
    return pfsIsrTask.PfsIsrTask(config=config)


class GlitchCorrectionTestCase(lsst.utils.tests.TestCase):
    def test_loadCorrectGlitchChannels(self):
        task = _makeTask()
        allChannels = tuple(range(32))
        # every channel is in correction scope ...
        self.assertEqual(task.loadCorrectGlitchChannels("n4"), allChannels)
        self.assertEqual(task.loadCorrectGlitchChannels("n1"), allChannels)
        # ... except the known-bad ASIC channels under Hann (default IRPfilter):
        # n3 ch24's IRP is bypassed, injecting digital glitches -> excluded.
        self.assertEqual(task.loadCorrectGlitchChannels("n3"), tuple(c for c in range(32) if c != 24))

    def test_loadCorrectGlitchChannelsMedianIncludesBadChannel(self):
        # Under IRPfilter=-1 the per-column median cleans ch24's IRP, so its
        # glitches are correctable and the channel is no longer excluded.
        task = _makeTask(irpFilter=-1)
        self.assertEqual(task.loadCorrectGlitchChannels("n3"), tuple(range(32)))
        self.assertEqual(task.loadCorrectGlitchChannels("n4"), tuple(range(32)))

    def test_correctableGlitchMask(self):
        task = _makeTask()
        H, W, nDeltas, nChan = 64, 10, 24, 32
        ch = H // nChan            # 2 rows per channel
        r18, r16 = 18 * ch, 16 * ch
        glitch = np.zeros((H, W, nDeltas), dtype=bool)
        cr = np.zeros((H, W, nDeltas), dtype=bool)
        deltas = np.zeros((H, W, nDeltas), dtype=np.float32)
        bad = np.zeros((H, W), dtype=bool)
        # deltas are otherwise flat -> IQR sigma 0 -> threshold = nSigma*floor = 5*8 = 40

        def pair(y, x, k, a, ret):
            glitch[y, x, k] = True
            glitch[y, x, k + 1] = True
            deltas[y, x, k] = a
            deltas[y, x, k + 1] = ret

        # A: a clean outlier+return pair -> correctable
        pair(r18, 0, 5, 60.0, -55.0)
        # B: several clean pairs (any number) -> correctable
        pair(r18, 1, 3, 60.0, -55.0)
        pair(r18, 1, 12, -58.0, 52.0)
        # T: outlier + return + a small spurious third flag (the triple-run case) -> correctable
        glitch[r18, 2, 5:8] = True
        deltas[r18, 2, 5:8] = (60.0, -52.0, 10.0)
        # R: outlier + a sub-threshold (below 40) but ~equal return -> correctable
        pair(r18, 3, 5, 60.0, -36.0)
        # L: a lone outlier with no return (CR-like) -> not correctable
        glitch[r18, 4, 5] = True
        deltas[r18, 4, 5] = 60.0
        # S: an outlier whose only glitch neighbour is same-sign (no valid return) -> not correctable
        glitch[r18, 5, 5:7] = True
        deltas[r18, 5, 5:7] = (60.0, 55.0)
        # F: a clean pair with an adjacent CR (ambiguous) -> not correctable
        pair(r18, 6, 5, 60.0, -55.0)
        cr[r18, 6, 7] = True
        deltas[r18, 6, 7] = 45.0
        # E: a clean pair but the pixel is already bad -> not correctable
        pair(r18, 8, 5, 60.0, -55.0)
        bad[r18, 8] = True
        # D: a clean pair on another channel (ch16) -- all channels are in scope -> correctable
        pair(r16, 0, 5, 60.0, -55.0)

        mask = task.correctableGlitchMask("n4", glitch, cr, deltas, bad, nChannels=nChan)
        self.assertEqual(mask.shape, (H, W))
        self.assertTrue(mask[r18, 0], "clean outlier+return pair should be correctable")
        self.assertTrue(mask[r18, 1], "several clean pairs should be correctable")
        self.assertTrue(mask[r18, 2], "outlier + return + a small spurious flag should be correctable")
        self.assertTrue(mask[r18, 3], "outlier + a sub-threshold ~equal return should be correctable")
        self.assertFalse(mask[r18, 4], "a lone outlier with no return should not be correctable")
        self.assertFalse(mask[r18, 5], "an outlier with no opposite-sign return should not be correctable")
        self.assertFalse(mask[r18, 6], "a CR adjacent to a glitch should not be correctable")
        self.assertFalse(mask[r18, 8], "already-bad pixel should not be correctable")
        self.assertTrue(mask[r16, 0], "a clean pair on any in-scope channel should be correctable")

    def test_correctGlitchScopeExcludesBadChannel(self):
        task = _makeTask()
        H, W, nDeltas, nChan = 64, 4, 24, 32
        ch = H // nChan
        glitch = np.zeros((H, W, nDeltas), dtype=bool)
        cr = np.zeros((H, W, nDeltas), dtype=bool)
        deltas = np.zeros((H, W, nDeltas), dtype=np.float32)
        bad = np.zeros((H, W), dtype=bool)

        def pair(y, x, k, a, ret):
            glitch[y, x, k] = True
            glitch[y, x, k + 1] = True
            deltas[y, x, k] = a
            deltas[y, x, k + 1] = ret

        pair(24 * ch, 0, 5, 60.0, -55.0)   # clean pair on n3 ch24 (digital bad-list channel)
        pair(0, 0, 5, 60.0, -55.0)         # clean pair on n3 ch0
        mask = task.correctableGlitchMask("n3", glitch, cr, deltas, bad, nChannels=nChan)
        self.assertFalse(mask[24 * ch, 0], "n3/ch24 is excluded from the correction scope")
        self.assertTrue(mask[0, 0], "other n3 channels are in the correction scope")


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(__import__("sys").modules["__main__"])
    unittest.main(failfast=True)
