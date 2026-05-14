import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs.h4Linearity import cr


def _flatRamp(nReads=16, H=8, W=8, rate=10.0, dtype=np.float32):
    """Build a noiseless cumulative ramp: cumulative[k, y, x] = k * rate."""
    return (
        np.arange(nReads, dtype=np.float64)[:, None, None]
        * rate
        * np.ones((1, H, W), dtype=np.float64)
    ).astype(dtype)


def _runCR(flux, **kwargs):
    """Test helper: run the deltas-based detector on a cumulative cube.

    ``cr.iterativeUtrDetectAndRepair`` now takes a delta cube — the
    production caller in ``isrTask.makeNirExposure`` does the
    ``np.diff`` once and the reconstructing ``np.cumsum`` once (only
    when needed). Tests still construct cubes (it's the natural way to
    express "inject a CR at read k") and check the post-repair cube,
    so we wrap the diff/cumsum here to keep the existing assertions.

    Modifies ``flux`` in place when ``repair=True`` (default).
    """
    repair = kwargs.get("repair", True)
    deltas = np.diff(flux, axis=0)
    read0 = flux[0:1].copy() if repair else None
    result = cr.iterativeUtrDetectAndRepair(deltas, **kwargs)
    if repair:
        flux[0:1] = read0
        np.cumsum(deltas, axis=0, out=flux[1:])
        flux[1:] += read0
    return result


def _injectCR(flux, y, x, k, amount):
    """Add a CR contribution at delta index k: cumulative reads >= k+1 jump by amount."""
    flux[k + 1:, y, x] += amount


def _injectGlitchPair(flux, y, x, k, amount):
    """Add a single-read digital glitch: read k is offset by `amount`.

    In delta space this creates an up/down pair at delta indices k-1 and k:
    delta[k-1] = rate + amount, delta[k] = rate - amount.
    """
    flux[k, y, x] += amount


class IterativeUtrDetectAndRepairTestCase(lsst.utils.tests.TestCase):

    def testNoDefectsReturnsClean(self):
        flux = _flatRamp(nReads=20, rate=10.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        glitchMask = np.zeros_like(good)
        before = flux.copy()

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
        )

        self.assertEqual(result.nCRs, 0)
        self.assertEqual(result.nGlitchPairs, 0)
        self.assertFalse(result.crFlagMask.any())
        self.assertFalse(result.glitchFlagMask.any())
        np.testing.assert_array_equal(flux, before)

    def testSingleCRDetectedAndRepaired(self):
        flux = _flatRamp(nReads=20, rate=10.0)
        _injectCR(flux, y=3, x=4, k=7, amount=300.0)
        good = np.ones(flux.shape[1:], dtype=bool)

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=None,
        )

        self.assertEqual(result.nCRs, 1)
        self.assertEqual(result.nGlitchPairs, 0)
        self.assertTrue(result.crFlagMask[7, 3, 4])
        # Repair should restore the ramp to the clean linear line.
        np.testing.assert_allclose(flux, _flatRamp(nReads=20, rate=10.0), atol=1.0)

    def testGlitchOnActiveChannelClassifiedAsGlitch(self):
        flux = _flatRamp(nReads=20, rate=10.0)
        # ASIC glitch on (y=2, x=2): read 6 is offset by +400. Delta 5
        # jumps up, delta 6 jumps down by the same amount.
        _injectGlitchPair(flux, y=2, x=2, k=6, amount=400.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        glitchMask = np.zeros_like(good)
        glitchMask[2, :] = True   # row 2 is glitch-active

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
        )

        self.assertEqual(result.nCRs, 0)
        self.assertEqual(result.nGlitchPairs, 1)
        self.assertTrue(result.glitchFlagMask[5, 2, 2])
        self.assertTrue(result.glitchFlagMask[6, 2, 2])
        # Repair should remove the digital offset.
        np.testing.assert_allclose(flux, _flatRamp(nReads=20, rate=10.0), atol=1.0)

    def testGlitchOnInactiveChannelClassifiedAsCR(self):
        """Same digital pair as above but on a row NOT in glitchPixelMask.

        With no glitch classification active, both deltas pass the single-read
        check; only the positive-residual one is flagged as a CR.
        """
        flux = _flatRamp(nReads=20, rate=10.0)
        _injectGlitchPair(flux, y=2, x=2, k=6, amount=400.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        glitchMask = np.zeros_like(good)   # row 2 is NOT glitch-active

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
        )

        self.assertEqual(result.nGlitchPairs, 0)
        self.assertGreaterEqual(result.nCRs, 1)
        # The positive-residual delta gets the CR flag.
        self.assertTrue(result.crFlagMask[5, 2, 2])

    def testGlitchMaskNoneDisablesGlitchDetection(self):
        flux = _flatRamp(nReads=20, rate=10.0)
        _injectGlitchPair(flux, y=2, x=2, k=6, amount=400.0)
        good = np.ones(flux.shape[1:], dtype=bool)

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=None,
        )

        self.assertEqual(result.nGlitchPairs, 0)

    def testRepairFalseLeavesCubeUnchanged(self):
        flux = _flatRamp(nReads=20, rate=10.0)
        _injectCR(flux, y=3, x=4, k=7, amount=300.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        before = flux.copy()

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=None, repair=False,
        )

        # Still flagged ...
        self.assertEqual(result.nCRs, 1)
        # ... but cube untouched.
        np.testing.assert_array_equal(flux, before)

    def testIterationsConverge(self):
        """nByIteration should monotonically diminish and stop before maxIterations."""
        flux = _flatRamp(nReads=20, rate=10.0)
        _injectCR(flux, y=3, x=4, k=7, amount=300.0)
        _injectGlitchPair(flux, y=5, x=6, k=10, amount=400.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        glitchMask = np.zeros_like(good)
        glitchMask[5, :] = True

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
            maxIterations=10,
        )

        self.assertGreaterEqual(result.nIterations, 1)
        self.assertLessEqual(result.nIterations, 10)
        self.assertEqual(len(result.nByIteration), result.nIterations)
        # All defects accounted for in iter 1 on a clean synthetic ramp.
        nc1, ng1 = result.nByIteration[0]
        self.assertGreaterEqual(nc1, 1)
        self.assertGreaterEqual(ng1, 1)
        # Subsequent iterations add nothing.
        for nc, ng in result.nByIteration[1:]:
            self.assertEqual(nc, 0)
            self.assertEqual(ng, 0)

    def testAsymmetricGlitchPairPaired(self):
        """Glitch where the recovery delta is just under threshold still pairs.

        The matched-pair criterion only requires ONE of the adjacent
        deltas to be independently flagged, plus opposite-sign residuals
        that cancel within threshold.

        Threshold here = nSigma * sigmaFloor = 5 * 8 = 40 ADU. Construct
        a pair with delta[5] residual = +50 (flagged) and delta[6]
        residual = -35 (NOT flagged), opposite signs, sum = +15 < 40
        (cancels). This must be classified as a glitch, not a CR.
        """
        flux = _flatRamp(nReads=20, rate=10.0)
        # Read 6 offset by +50 → delta[5] = rate+50, delta[6] = rate-50.
        # Bring reads 7+ back by 15 → delta[6] = rate-35.
        flux[6, 4, 4] += 50.0
        flux[7:, 4, 4] += 15.0
        good = np.ones(flux.shape[1:], dtype=bool)
        glitchMask = np.ones_like(good)

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
        )

        self.assertEqual(result.nCRs, 0)
        self.assertEqual(result.nGlitchPairs, 1)
        self.assertTrue(result.glitchFlagMask[5, 4, 4])
        self.assertTrue(result.glitchFlagMask[6, 4, 4])

    def testShortRampGlitchDoesNotDriveRateRunaway(self):
        """A single big glitch on a short (22-read) ramp must NOT bias
        the rate estimator into the runaway regime that false-flags
        every other delta.

        With LSQ slope, a single −5000 ADU outlier near the center of a
        22-read ramp drags the slope by ~−10 ADU/read, which is enough
        to push typical (rate≈0) deltas to near the 5σ threshold,
        triggering iter-2 false flags and a divergent rate. Median is
        robust to a few outliers; rate stays near 0.
        """
        flux = _flatRamp(nReads=22, rate=0.0)
        flux[11:, 4, 4] -= 5000.0  # big downward glitch at read 11
        flux[12:, 4, 4] += 5000.0  # recovery at read 12
        good = np.ones(flux.shape[1:], dtype=bool)
        glitchMask = np.ones_like(good)

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
        )

        # Final rate must be near the true 0, not the runaway-divergent
        # value the LSQ slope would produce.
        self.assertLess(abs(float(result.rate[4, 4])), 5.0)
        # Exactly one glitch pair detected; no false CRs.
        self.assertEqual(result.nGlitchPairs, 1)
        self.assertEqual(result.nCRs, 0)

    def testGlitchRepairPreservesPostGlitchCumulative(self):
        """Repair of a glitch pair must leave cube[k+2] (and beyond) unchanged.

        A digital ASIC glitch corrupts only ONE cumulative value (cube[k+1]);
        the recovery delta has restored the true cumulative at read k+2.
        Replacing both deltas with `rate` (naive repair) erases real signal
        in the second delta and shifts everything that follows. The correct
        repair fixes only cube[k+1].
        """
        flux = _flatRamp(nReads=20, rate=5.0)
        # Inject some real signal evolution after the glitch — a real flux
        # event at reads 10..14 — that the repair must preserve.
        flux[10:15, 4, 4] += np.arange(5, dtype=np.float32) * 3.0
        flux[15:, 4, 4] += 12.0
        # Inject a downward glitch at read 6 (large offset −500):
        flux[6, 4, 4] -= 500.0
        expectedAfterGlitch = flux[7:, 4, 4].copy()
        # Repaired cube[6] should be cube[5] + rate (smooth interpolation).
        expectedRead6 = flux[5, 4, 4] + 5.0

        good = np.ones(flux.shape[1:], dtype=bool)
        glitchMask = np.ones_like(good)

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
        )

        self.assertEqual(result.nGlitchPairs, 1)
        # The single corrupted read is restored to the smooth ramp.
        self.assertAlmostEqual(float(flux[6, 4, 4]), float(expectedRead6), delta=1.0)
        # Reads after the glitch are byte-preserved (the recovery delta
        # already restored the true cumulative).
        np.testing.assert_allclose(
            flux[7:, 4, 4], expectedAfterGlitch, atol=0.01,
            err_msg="Reads after a glitch pair must match the original cube.",
        )

    def testGlitchAmplitudeMinADUFloor(self):
        """A small glitch pair stops being classified as glitch when the
        amplitude floor exceeds its residual size.

        With default (floor=0), a small but still-above-threshold pair
        IS classified as glitch. Raising the floor above its amplitude
        suppresses the glitch classification (the positive-residual
        delta falls back to CR).
        """
        good = np.ones((4, 4), dtype=bool)
        glitchMask = np.ones_like(good)

        # Small symmetric glitch: read 6 offset by +60 → delta[5] resid +50,
        # delta[6] resid −50. Threshold (5*sigmaFloor=40) is comfortably
        # exceeded, sum cancels to 0.
        flux = _flatRamp(nReads=20, H=4, W=4, rate=10.0)
        flux[6, 1, 1] += 60.0
        # Floor=0: classified as glitch.
        result = _runCR(
            flux.copy(), goodPixelMask=good, glitchPixelMask=glitchMask,
        )
        self.assertEqual(result.nGlitchPairs, 1)
        self.assertEqual(result.nCRs, 0)

        # Floor=80 ADU (above the 50-ADU residual): the pair is no longer
        # classified, so the positive delta falls back to CR.
        result = _runCR(
            flux.copy(), goodPixelMask=good, glitchPixelMask=glitchMask,
            glitchAmplitudeMinADU=80.0,
        )
        self.assertEqual(result.nGlitchPairs, 0)
        self.assertGreaterEqual(result.nCRs, 1)

    def testBoundaryGlitchClassifiedAsGlitch(self):
        """A flagged delta at index 0 or N-2 (no neighbor to pair with)
        should be marked as ASIC_GLITCH, not CR or unclassified.

        Digital glitches that hit the first or last read look like a
        single huge delta at the ramp boundary and have no pair partner.
        The boundary heuristic catches them anyway.
        """
        good = np.ones((4, 4), dtype=bool)
        glitchMask = np.ones_like(good)

        # Glitch at FIRST read: corrupt cube[0] by +200. delta[0] becomes
        # rate - 200 → huge negative residual.
        flux = _flatRamp(nReads=20, H=4, W=4, rate=5.0)
        flux[0, 1, 1] += 200.0
        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
        )
        self.assertTrue(result.glitchFlagMask[0, 1, 1],
                        'first-read glitch should hit the ASIC_GLITCH plane')
        self.assertFalse(result.crFlagMask[0, 1, 1],
                         'first-read glitch must NOT also be CR-flagged')

        # Glitch at LAST read: corrupt cube[N-1] by -200. delta[N-2]
        # becomes rate - 200 → huge negative residual at the very last
        # delta index.
        flux = _flatRamp(nReads=20, H=4, W=4, rate=5.0)
        flux[-1, 2, 2] -= 200.0
        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=glitchMask,
        )
        self.assertTrue(result.glitchFlagMask[-1, 2, 2],
                        'last-read glitch should hit the ASIC_GLITCH plane')
        self.assertFalse(result.crFlagMask[-1, 2, 2],
                         'last-read glitch must NOT also be CR-flagged')

    def testGoodPixelMaskExcludesPixels(self):
        flux = _flatRamp(nReads=20, rate=10.0)
        _injectCR(flux, y=3, x=4, k=7, amount=300.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        good[3, 4] = False   # exclude the CR pixel

        result = _runCR(
            flux, goodPixelMask=good, glitchPixelMask=None,
        )

        self.assertEqual(result.nCRs, 0)


class RampQRTestCase(lsst.utils.tests.TestCase):
    """Sanity tests for ``cr._rampQR``.

    Equivalence to ``np.percentile(deltas, [25, 50, 75], axis=0)`` with
    method='linear' is the contract; if it ever drifts, the CR detector's
    threshold changes silently.
    """

    def _checkMatchesNumpy(self, deltas):
        p25, median, p75 = cr._rampQR(deltas)
        ref = np.percentile(deltas, [25, 50, 75], axis=0)
        np.testing.assert_allclose(p25, ref[0], rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(median, ref[1].astype(np.float32),
                                   rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(p75, ref[2], rtol=1e-6, atol=1e-6)
        self.assertEqual(median.dtype, np.float32)

    def testOddNMatchesNumpy(self):
        # 87 deltas = 88-read dark, the production case (odd N -> exact median).
        rng = np.random.default_rng(0)
        deltas = rng.standard_normal((87, 4, 5)).astype(np.float32)
        self._checkMatchesNumpy(deltas)

    def testEvenNMatchesNumpy(self):
        rng = np.random.default_rng(1)
        deltas = rng.standard_normal((88, 4, 5)).astype(np.float32)
        self._checkMatchesNumpy(deltas)

    def testSmallNMatchesNumpy(self):
        # Edge: minimum N for which the CR loop is even allowed (3 reads
        # -> 2 deltas) and a couple of small N where the rank set
        # collapses.
        rng = np.random.default_rng(2)
        for N in (2, 3, 4, 5):
            deltas = rng.standard_normal((N, 3, 3)).astype(np.float32)
            with self.subTest(N=N):
                self._checkMatchesNumpy(deltas)

    def testSubsetShapeMatchesNumpy(self):
        # The iter 2+ subset path calls _rampQR with shape (N-1, 1, M).
        rng = np.random.default_rng(3)
        deltas = rng.standard_normal((87, 1, 17)).astype(np.float32)
        self._checkMatchesNumpy(deltas)

    def testConstantArrayIsExact(self):
        # All ranks equal -> p25 == median == p75 == the constant.
        deltas = np.full((87, 3, 3), 7.5, dtype=np.float32)
        p25, median, p75 = cr._rampQR(deltas)
        np.testing.assert_array_equal(p25, 7.5)
        np.testing.assert_array_equal(median, 7.5)
        np.testing.assert_array_equal(p75, 7.5)

    def testDoesNotMutateInput(self):
        # np.partition is not in place when called as a function; verify
        # _rampQR preserves the caller's deltas array (the residual pass
        # reads it after _rampQR returns).
        rng = np.random.default_rng(4)
        deltas = rng.standard_normal((87, 4, 5)).astype(np.float32)
        snapshot = deltas.copy()
        cr._rampQR(deltas)
        np.testing.assert_array_equal(deltas, snapshot)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
