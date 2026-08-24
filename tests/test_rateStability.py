import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs.h4Linearity import chunking, rateStability


def _deltas(rate1, rate2, H=4, W=4, nDeltas=10):
    """Construct a noiseless delta cube with two-half rates rate1 and rate2.

    Returns ``(H, W, nDeltas)`` -- time axis last to match the H4 ISR
    cube convention.
    """
    nHalf = nDeltas // 2
    d = np.empty((H, W, nDeltas), dtype=np.float32)
    d[..., :nHalf] = rate1
    d[..., nHalf:] = rate2
    return d


class DetectRateInstabilityTestCase(lsst.utils.tests.TestCase):

    def testFlatRampNotRejected(self):
        deltas = np.full((4, 4, 10), 100.0, dtype=np.float32)
        flags = np.zeros_like(deltas, dtype=bool)
        result = rateStability.detectRateInstability(deltas, flags)
        self.assertEqual(result.nRejected, 0)
        finite = np.isfinite(result.fraction)
        self.assertTrue(np.all(result.fraction[finite] < 1e-6))

    def testSmallRateChangeNotRejected(self):
        # rate1=100, rate2=85 -> fraction = 0.15 < 0.20 default threshold.
        deltas = _deltas(100.0, 85.0)
        flags = np.zeros_like(deltas, dtype=bool)
        result = rateStability.detectRateInstability(deltas, flags)
        self.assertTrue(np.allclose(result.fraction, 0.15, atol=1e-5))
        self.assertEqual(result.nRejected, 0)

    def testLargeRateChangeRejected(self):
        # rate1=100, rate2=50 -> fraction = 0.50 > 0.20.
        deltas = _deltas(100.0, 50.0)
        flags = np.zeros_like(deltas, dtype=bool)
        result = rateStability.detectRateInstability(deltas, flags)
        self.assertTrue(np.allclose(result.fraction, 0.50, atol=1e-5))
        self.assertEqual(result.nRejected, 4 * 4)
        self.assertTrue(np.all(result.rejectMask))

    def testSignCrossingRejected(self):
        # rate1=20, rate2=-10 -> |diff|=30, max(20,10,5)=20 -> fraction=1.5.
        deltas = _deltas(20.0, -10.0)
        flags = np.zeros_like(deltas, dtype=bool)
        result = rateStability.detectRateInstability(deltas, flags)
        self.assertTrue(np.allclose(result.fraction, 1.5, atol=1e-5))
        self.assertTrue(np.all(result.rejectMask))

    def testDarkPixelFloorProtection(self):
        # rate1=0.1, rate2=-0.5 -> |diff|=0.6, denom clamped to rateFloorADU=5.
        # fraction = 0.12 < 0.20 -> not rejected.
        deltas = _deltas(0.1, -0.5)
        flags = np.zeros_like(deltas, dtype=bool)
        result = rateStability.detectRateInstability(deltas, flags)
        self.assertTrue(np.allclose(result.fraction, 0.12, atol=1e-5))
        self.assertEqual(result.nRejected, 0)

    def testUntestableWhenSegmentTooShort(self):
        # 10 deltas, default minDeltasPerSegment=3. Flag deltas 0..3 of pixel
        # (1,1) so segment 1 has 1 un-flagged delta -> untestable.
        deltas = _deltas(100.0, 100.0)
        flags = np.zeros_like(deltas, dtype=bool)
        flags[1, 1, 0:4] = True
        result = rateStability.detectRateInstability(deltas, flags)
        self.assertLess(int(result.nTestable[1, 1]), 2)
        self.assertFalse(result.rejectMask[1, 1])
        self.assertTrue(np.isnan(result.fraction[1, 1]))
        self.assertGreaterEqual(result.nUntestable, 1)

    def testFlaggedDeltaExcludedFromMean(self):
        # Two large spikes in segment 1 at pixel (2,2): unflagged -> drags
        # mean high enough that adjusted fraction > threshold -> rejected.
        # Flagged -> excluded -> rates stay equal -> not rejected.
        # Two spikes are needed because a single spike in a uniform background
        # yields sem == diff exactly, collapsing the noise-adjusted numerator
        # to zero.
        deltas = _deltas(100.0, 100.0, nDeltas=20)
        deltas[2, 2, 2] = 1000.0
        deltas[2, 2, 3] = 1000.0
        flags = np.zeros_like(deltas, dtype=bool)
        rejected = rateStability.detectRateInstability(deltas, flags)
        self.assertTrue(rejected.rejectMask[2, 2])
        flags[2, 2, 2] = True
        flags[2, 2, 3] = True
        repaired = rateStability.detectRateInstability(deltas, flags)
        self.assertFalse(repaired.rejectMask[2, 2])

    def testGoodPixelMaskExcludes(self):
        deltas = _deltas(100.0, 50.0)
        flags = np.zeros_like(deltas, dtype=bool)
        good = np.ones((4, 4), dtype=bool)
        good[2, 2] = False
        result = rateStability.detectRateInstability(
            deltas, flags, goodPixelMask=good)
        self.assertFalse(result.rejectMask[2, 2])
        self.assertTrue(np.isnan(result.fraction[2, 2]))

    def testCustomThreshold(self):
        # rate1=100, rate2=85 -> fraction=0.15. threshold=0.10 -> reject;
        # threshold=0.20 (default) -> not rejected.
        deltas = _deltas(100.0, 85.0)
        flags = np.zeros_like(deltas, dtype=bool)
        low = rateStability.detectRateInstability(
            deltas, flags, threshold=0.10)
        self.assertTrue(np.all(low.rejectMask))
        high = rateStability.detectRateInstability(deltas, flags)
        self.assertEqual(high.nRejected, 0)

    def testNoisyStableBrightPixelNotRejected(self):
        # A bright (~100 ADU/read) pixel whose half-rate diff is purely
        # noise must not be rejected: quadrature subtraction should
        # collapse the numerator to ~0.
        rng = np.random.default_rng(seed=20260523)
        nDeltas = 40
        H, W = 4, 4
        deltas = rng.normal(loc=100.0, scale=10.0,
                            size=(H, W, nDeltas)).astype(np.float32)
        flags = np.zeros_like(deltas, dtype=bool)
        result = rateStability.detectRateInstability(deltas, flags)
        # No pixel should be rejected; all fractions are well below 0.20.
        self.assertEqual(result.nRejected, 0)
        finite = np.isfinite(result.fraction)
        self.assertTrue(np.all(result.fraction[finite] < 0.20))

    def testNoisyBrokenPixelRejected(self):
        # A pixel whose first-half mean is 100 and second-half mean is
        # 50 (with small noise) must be rejected: signal dominates the
        # noise subtraction.
        rng = np.random.default_rng(seed=20260523)
        nDeltas = 40
        H, W = 4, 4
        deltas = np.empty((H, W, nDeltas), dtype=np.float32)
        nHalf = nDeltas // 2
        deltas[..., :nHalf] = rng.normal(
            loc=100.0, scale=5.0, size=(H, W, nHalf))
        deltas[..., nHalf:] = rng.normal(
            loc=50.0, scale=5.0, size=(H, W, nDeltas - nHalf))
        flags = np.zeros_like(deltas, dtype=bool)
        result = rateStability.detectRateInstability(deltas, flags)
        # All pixels are clearly broken; fractions are well above 0.20.
        self.assertEqual(result.nRejected, H * W)
        self.assertTrue(np.all(result.fraction > 0.20))

    def testSegmentRatesArePlainMeanOfUnflaggedDeltas(self):
        # Plain-mean estimator, production-consistent. With a flag, the mean
        # is recomputed over the un-flagged deltas only.
        deltas = _deltas(100.0, 60.0)
        flags = np.zeros_like(deltas, dtype=bool)
        flags[0, 0, 0] = True  # flag one delta in the first half at pixel (0,0)
        result = rateStability.detectRateInstability(deltas, flags)
        # Pixel (0,0) first half [index 0]: 4 unflagged deltas of 100 -> mean 100.
        self.assertAlmostEqual(float(result.segmentRates[0, 0, 0]), 100.0,
                               places=3)
        # Pixel (1,1) (no flags) first half [index 0]: 5 deltas of 100 -> mean 100.
        self.assertAlmostEqual(float(result.segmentRates[0, 1, 1]), 100.0,
                               places=3)
        # Second half [index 1]: 5 deltas of 60 -> mean 60.
        self.assertAlmostEqual(float(result.segmentRates[1, 0, 0]), 60.0,
                               places=3)

    def testGoodPixelMaskEchoedBack(self):
        # The result echoes back the input goodPixelMask so a downstream
        # diagnostic can rerun the metric at different knobs against the
        # exact pixel set the pipeline used.
        deltas = _deltas(100.0, 50.0)
        flags = np.zeros_like(deltas, dtype=bool)
        good = np.ones((4, 4), dtype=bool)
        good[2, 2] = False
        result = rateStability.detectRateInstability(
            deltas, flags, goodPixelMask=good)
        self.assertEqual(result.goodPixelMask.shape, (4, 4))
        self.assertTrue(np.array_equal(result.goodPixelMask, good))

    def testRejectsNon3DInput(self):
        flat = np.zeros((10, 10), dtype=np.float32)
        with self.assertRaises(ValueError):
            rateStability.detectRateInstability(flat, flat.astype(bool))

    def testRejectsFewerDeltasThanRequired(self):
        # nDeltas=4 < 2*minDeltasPerSegment=6.
        deltas = np.zeros((4, 4, 4), dtype=np.float32)
        with self.assertRaises(ValueError):
            rateStability.detectRateInstability(
                deltas, np.zeros_like(deltas, dtype=bool))

    def testRejectsFlagMaskShapeMismatch(self):
        deltas = np.zeros((4, 4, 10), dtype=np.float32)
        flags = np.zeros((4, 5, 10), dtype=bool)
        with self.assertRaises(ValueError):
            rateStability.detectRateInstability(deltas, flags)

    def testRejectsGoodPixelMaskShapeMismatch(self):
        deltas = np.zeros((4, 4, 10), dtype=np.float32)
        flags = np.zeros_like(deltas, dtype=bool)
        good = np.ones((4, 5), dtype=bool)
        with self.assertRaises(ValueError):
            rateStability.detectRateInstability(
                deltas, flags, goodPixelMask=good)


class RowChunkingTestCase(lsst.utils.tests.TestCase):
    """Banding the segment accumulation must not change a single bit.

    The sums reduce along the delta axis per pixel, so how many image rows are
    handled together is invisible to the result. Banding is what keeps the
    float64 cast and its square -- two copies of half the delta cube, ~18.5 GB
    at production size -- bounded by the band.
    """

    def _inputs(self, H=29, W=7, nDeltas=17, seed=1885):
        rng = np.random.RandomState(seed)
        deltas = rng.normal(120.0, 15.0, size=(H, W, nDeltas)).astype(np.float32)
        # A rate step in the second half of some pixels, plus flagged deltas.
        deltas[3, 2, nDeltas // 2:] += 400.0
        deltas[20, 5, nDeltas // 2:] -= 250.0
        flagMask = np.zeros((H, W, nDeltas), dtype=bool)
        flagMask[7, 1, 2] = True
        flagMask[7, 1, 3] = True
        flagMask[25, 6, nDeltas - 2] = True
        goodPixelMask = np.ones((H, W), dtype=bool)
        goodPixelMask[11, 4] = False
        return deltas, flagMask, goodPixelMask

    def _run(self, chunkBytes):
        deltas, flagMask, goodPixelMask = self._inputs()
        original = chunking.CHUNK_BYTES
        chunking.CHUNK_BYTES = chunkBytes
        try:
            return rateStability.detectRateInstability(
                deltas, flagMask, goodPixelMask=goodPixelMask)
        finally:
            chunking.CHUNK_BYTES = original

    def testEveryChunkSizeAgreesExactly(self):
        reference = self._run(1 << 30)
        # 29 rows: none of these band sizes divides it evenly.
        for chunkBytes in (1 << 30, 1 << 12, 1 << 9, 1):
            with self.subTest(chunkBytes=chunkBytes):
                got = self._run(chunkBytes)
                np.testing.assert_array_equal(got.rejectMask,
                                              reference.rejectMask)
                np.testing.assert_array_equal(got.nTestable,
                                              reference.nTestable)
                np.testing.assert_array_equal(got.segmentRates,
                                              reference.segmentRates)
                # fraction carries NaN at untested pixels; compare bitwise.
                self.assertTrue(
                    np.array_equal(got.fraction, reference.fraction,
                                   equal_nan=True))
                self.assertEqual(got.nRejected, reference.nRejected)
                self.assertEqual(got.nUntestable, reference.nUntestable)

    def testFixtureActuallyRejectsSomething(self):
        # Otherwise the comparison above could pass on an all-empty result.
        result = self._run(1 << 30)
        self.assertGreater(result.nRejected, 0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
