"""Tests for ``PfsIsrTask.calcUTRrates`` on an asymmetric ramp with
position-dependent per-pixel slopes.

The UTR weights are the standard linear-LS slope estimator on evenly
spaced reads: with ``cube[i, y, x] = (i + 1) * slope[y, x]`` (no
offset) and exact arithmetic, ``calcUTRrates(cube)[y, x] == slope[y, x]``.
We use a position-encoded slope ``slope[y, x] = y * 10 + x`` and
asymmetric ``(N=15, H=4, W=7)`` so that:

- An ``H <-> W`` swap inside ``calcUTRrates`` would either raise (shape
  mismatch between cube[0] and the rate buffer) or produce a transposed
  rate plane whose values are observably wrong at every asymmetric
  ``(y, x)``.
- An axis-0 vs axis-1 swap (treating reads as height) would produce a
  shape mismatch immediately.
- A read-direction swap (cube[N-1-i] iterated instead of cube[i]) would
  invert the sign of the recovered slope and is also caught.
"""
import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs import isrTask as pfsIsrTask
from lsst.obs.pfs.h4Linearity import cr


N_READS = 15
H = 4
W = 7


def _slopePerPixel(H, W, dtype=np.float32):
    y = np.arange(H, dtype=dtype)[:, None]
    x = np.arange(W, dtype=dtype)[None, :]
    return y * 10.0 + x


def _linearRamp(nReads, H, W, dtype=np.float32):
    """``cube[i, y, x] = (i + 1) * slope[y, x]``."""
    slope = _slopePerPixel(H, W, dtype=dtype)
    i = np.arange(1, nReads + 1, dtype=dtype)[:, None, None]
    return (i * slope[None, :, :]).astype(dtype)


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


class CalcUTRRatesTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.task = _makeIsrTask()

    def testRecoversPerPixelSlope(self):
        cube = _linearRamp(N_READS, H, W)
        rates = self.task.calcUTRrates(cube)
        expected = _slopePerPixel(H, W)
        self.assertEqual(rates.shape, (H, W))
        np.testing.assert_allclose(
            rates, expected, atol=1e-3,
            err_msg="calcUTRrates must recover the per-pixel linear slope",
        )

    def testValuesAtAsymmetricPixels(self):
        # Spot-check a few asymmetric (y, x) so a silent H<->W swap shows
        # up as a value mismatch even if the rate plane shape "looks
        # right" after the swap.
        cube = _linearRamp(N_READS, H, W)
        rates = self.task.calcUTRrates(cube)
        for y, x in [(0, 0), (3, 0), (0, 6), (2, 5), (3, 6)]:
            self.assertAlmostEqual(
                float(rates[y, x]),
                y * 10.0 + x,
                places=3,
                msg=f"rate at (y={y}, x={x}) must equal y*10 + x",
            )

    def testReadOrderMatters(self):
        # The UTR weights are anti-symmetric about the ramp midpoint:
        # reversing the read order along axis 0 negates the recovered
        # slope. Confirms axis 0 is actually being treated as the read
        # index (not the spatial axis).
        cube = _linearRamp(N_READS, H, W)
        cubeRev = cube[::-1].copy()
        rates = self.task.calcUTRrates(cube)
        ratesRev = self.task.calcUTRrates(cubeRev)
        np.testing.assert_allclose(
            ratesRev, -rates, atol=1e-3,
            err_msg="reversing reads must negate the recovered slope",
        )

    def testNreadsTruncatesAxisZero(self):
        # ``nreads=k`` should only use the first k reads (axis 0). With
        # cube[i, y, x] = (i+1) * slope, the first k reads still satisfy
        # the linear model exactly, so the recovered slope is unchanged.
        cube = _linearRamp(N_READS, H, W)
        rates_full = self.task.calcUTRrates(cube)
        rates_truncated = self.task.calcUTRrates(cube, nreads=N_READS - 5)
        np.testing.assert_allclose(
            rates_truncated, rates_full, atol=1e-3,
            err_msg="truncating reads (still a linear ramp) must give "
                    "the same slope",
        )

    def testShapeMatchesSpatialPlane(self):
        # The output shape must equal cube[0].shape (i.e. (H, W)) — not
        # some transposition. If axes 1 and 2 were swapped inside the
        # function the shape would be (W, H) here.
        cube = _linearRamp(N_READS, H, W)
        rates = self.task.calcUTRrates(cube)
        self.assertEqual(rates.shape, cube[0].shape)
        self.assertEqual(rates.shape, (H, W))

    def testThreeRateFunctionsAgreeOnCleanRamp(self):
        """Cross-check: ``calcUTRrates`` (read-space), ``calcUTRrateFrom-
        Deltas`` (delta-space closed form), and the CR detector's
        ``IterativeRepairResult.rate`` are mathematically tied —
        ``u[j] = Σ_{i>j} w[i] = 6(j+1)(N-1-j) / (N(N-1)(N+1))``. On a
        clean ramp with no flags they must produce identical per-pixel
        rates, within float32 round-off. Catches silent drift if anyone
        re-derives one of the formulas.
        """
        cube = _linearRamp(N_READS, H, W)
        cubeHWN = np.ascontiguousarray(cube.transpose(1, 2, 0))
        deltas = np.diff(cubeHWN, axis=-1)

        rRead = self.task.calcUTRrates(cube)
        rDelta = self.task.calcUTRrateFromDeltas(deltas)
        result = cr.iterativeUtrDetectAndRepair(
            deltas.copy(),
            goodPixelMask=np.ones(cubeHWN.shape[:-1], dtype=bool),
            glitchPixelMask=None,
        )
        rCR = result.rate

        np.testing.assert_allclose(
            rDelta, rRead, rtol=1e-5,
            err_msg="calcUTRrateFromDeltas must match calcUTRrates on a "
                    "clean linear ramp",
        )
        np.testing.assert_allclose(
            rCR, rRead, rtol=1e-5,
            err_msg="iterResult.rate must match calcUTRrates on a clean "
                    "linear ramp (no flags fired)",
        )


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
