import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs.h4Linearity import cr, sim


class RampGenerationTestCase(lsst.utils.tests.TestCase):

    def testCleanRampIsLinearOnAverage(self):
        params = sim.RampParams(nReads=30, H=8, W=8, rate=100.0,
                                readNoise=1.0, poisson=False)
        cube = sim.makeRamp(params, rng=0)
        self.assertEqual(cube.shape, (30, 8, 8))
        self.assertEqual(cube.dtype, np.float32)
        # Per-pixel slope ~ rate (low-noise case).
        for y in range(8):
            for x in range(8):
                slope = float(np.polyfit(np.arange(30), cube[:, y, x], 1)[0])
                self.assertAlmostEqual(slope, 100.0, delta=0.3)

    def testNoisePropagates(self):
        params = sim.RampParams(nReads=50, H=32, W=32, rate=0.0,
                                readNoise=5.0, poisson=False)
        cube = sim.makeRamp(params, rng=42)
        # Empirical per-read std should match readNoise; cumulative-space
        # std of read 0 (vs the cube-mean prediction zero) is just readNoise.
        std_read0 = float(cube[0].std())
        self.assertAlmostEqual(std_read0, 5.0, delta=0.5)

    def testCRPersists(self):
        params = sim.RampParams(nReads=20, H=4, W=4, rate=10.0,
                                readNoise=0.0, poisson=False)
        crs = [sim.CR(y=2, x=3, read=10, amount=500.0)]
        cube = sim.makeRamp(params, crs=crs, rng=0)
        # New cumulative convention: cube[k] is the cumulative AFTER (k+1)
        # reads, i.e. (k+1)*rate. CR injected at read 10 adds to read 10's
        # per-read flux, so cube[10..] are bumped up by +500.
        for k in range(10):
            self.assertAlmostEqual(float(cube[k, 2, 3]),
                                   (k + 1) * 10.0, delta=0.5)
        for k in range(10, 20):
            self.assertAlmostEqual(float(cube[k, 2, 3]),
                                   (k + 1) * 10.0 + 500.0, delta=0.5)

    def testAsicGlitchIsSingleRead(self):
        params = sim.RampParams(nReads=20, H=4, W=4, rate=10.0,
                                readNoise=0.0, poisson=False)
        glitches = [sim.AsicGlitch(y=1, x=1, read=8, amount=2048.0)]
        cube = sim.makeRamp(params, glitches=glitches, rng=0)
        # cube[k] = (k+1)*rate; the glitch offsets ONE read (k=8).
        self.assertAlmostEqual(float(cube[7, 1, 1]), 8 * 10.0, delta=0.5)
        self.assertAlmostEqual(float(cube[8, 1, 1]),
                               9 * 10.0 + 2048.0, delta=0.5)
        self.assertAlmostEqual(float(cube[9, 1, 1]), 10 * 10.0, delta=0.5)
        # In delta space, this shows as +2048 then −2048.
        deltas = np.diff(cube[:, 1, 1])
        # delta[7] = read[8] - read[7] ≈ rate + glitch
        self.assertAlmostEqual(float(deltas[7]), 10.0 + 2048.0, delta=0.5)
        # delta[8] = read[9] - read[8] ≈ rate − glitch
        self.assertAlmostEqual(float(deltas[8]), 10.0 - 2048.0, delta=0.5)

    def testDigitalGlitchAmountsArePowersOf2(self):
        amounts = sim.digitalGlitchAmounts(
            n=200, bits=(10, 11, 12), signed=True, rng=42,
        )
        self.assertEqual(amounts.shape, (200,))
        for a in amounts:
            self.assertIn(abs(float(a)), {1024.0, 2048.0, 4096.0})
            self.assertIn(np.sign(float(a)), {-1.0, 1.0})

    def testRandomGenerators(self):
        crs = sim.randomCRs(n=20, H=8, W=8, nReads=30, rng=0)
        self.assertEqual(len(crs), 20)
        for c in crs:
            self.assertTrue(0 <= c.y < 8)
            self.assertTrue(0 <= c.x < 8)
            self.assertTrue(1 <= c.read < 29)
            self.assertTrue(50.0 <= c.amount <= 5000.0)

        glitches = sim.randomAsicGlitches(n=10, H=8, W=8, nReads=30, rng=0)
        self.assertEqual(len(glitches), 10)
        for g in glitches:
            self.assertTrue(0 <= g.y < 8)
            self.assertTrue(0 <= g.x < 8)
            self.assertTrue(1 <= g.read < 28)
            self.assertIn(abs(g.amount), {1024.0, 2048.0, 4096.0, 8192.0})


class NonlinearityTestCase(lsst.utils.tests.TestCase):

    def testLinearCaseIsIdentity(self):
        nl = sim.Nonlinearity(alpha=0.0)
        t = np.linspace(0, 50000, 1000, dtype=np.float32)
        np.testing.assert_array_equal(nl.forward(t), t)
        np.testing.assert_array_equal(nl.inverse(t), t)

    def testCompressedAtFullWell(self):
        nl = sim.Nonlinearity(alpha=0.05, qMax=60000.0)
        # forward at t=qMax is qMax * (1 − alpha) = 0.95 * qMax = 57000.
        m = nl.forward(np.array([60000.0], dtype=np.float32))
        self.assertAlmostEqual(float(m[0]), 57000.0, delta=1.0)

    def testRoundTripIsExact(self):
        nl = sim.Nonlinearity(alpha=0.05, qMax=60000.0)
        t = np.linspace(0.0, 50000.0, 5000, dtype=np.float32)
        recovered = nl.inverse(nl.forward(t))
        # Round-trip in float32 should be tight; allow a few ADU slop.
        np.testing.assert_allclose(recovered, t, atol=2.0)

    def testRoundTripWith3DCube(self):
        nl = sim.Nonlinearity(alpha=0.05, qMax=60000.0)
        cube = np.linspace(0.0, 50000.0, 100, dtype=np.float32).reshape(10, 5, 2)
        recovered = nl.inverse(nl.forward(cube))
        np.testing.assert_allclose(recovered, cube, atol=2.0)


class RawAndLinearRampsTestCase(lsst.utils.tests.TestCase):

    def testRawIsCompressedRelativeToLinearWhenAlphaPositive(self):
        params = sim.RampParams(nReads=40, H=4, W=4, rate=1500.0,
                                readNoise=0.0, poisson=False)
        nl = sim.Nonlinearity(alpha=0.05, qMax=60000.0)
        raw, lin = sim.makeRawAndLinearRamps(params, nonlinearity=nl,
                                              rng=0)
        # With cum[k] = (k+1)*rate, last read is 40*1500 = 60000 (full well).
        # raw at full well ≈ 60000 * (1 − 0.05) = 57000.
        self.assertAlmostEqual(float(lin[-1].mean()), 60000.0, delta=200.0)
        self.assertAlmostEqual(float(raw[-1].mean()), 57000.0, delta=300.0)
        # Linearized rate (slope) should recover the true rate.
        rate_est = sim.utrRate(lin)
        self.assertAlmostEqual(float(rate_est.mean()), params.rate, delta=10.0)

    def testCRInRawAndLinearMatchTruth(self):
        params = sim.RampParams(nReads=30, H=4, W=4, rate=100.0,
                                readNoise=0.0, poisson=False)
        nl = sim.Nonlinearity(alpha=0.05, qMax=60000.0)
        crs = [sim.CR(y=2, x=3, read=15, amount=2000.0)]
        raw, lin = sim.makeRawAndLinearRamps(params, nonlinearity=nl,
                                              crs=crs, rng=0)
        # CR injected in true-linear space, so the linearized cube should
        # show a clean +2000 step starting at read 15.
        delta_lin = float(lin[15, 2, 3] - lin[14, 2, 3])
        self.assertAlmostEqual(delta_lin, params.rate + 2000.0, delta=10.0)

    def testAsicGlitchInRawIsCleanSingleRead(self):
        params = sim.RampParams(nReads=20, H=4, W=4, rate=100.0,
                                readNoise=0.0, poisson=False)
        nl = sim.Nonlinearity(alpha=0.05, qMax=60000.0)
        g = sim.AsicGlitch(y=1, x=1, read=10, amount=2048.0)
        raw = sim.makeRawRamp(params, nonlinearity=nl, asicGlitches=[g], rng=0)
        # In raw: only read 10 is offset; reads 9 and 11 are at the
        # nominal compressed ramp. With cum[k] = (k+1)*rate convention.
        nominal_9 = float(nl.forward(np.array([10 * 100.0]))[0])
        nominal_10 = float(nl.forward(np.array([11 * 100.0]))[0])
        nominal_11 = float(nl.forward(np.array([12 * 100.0]))[0])
        self.assertAlmostEqual(float(raw[9, 1, 1]), nominal_9, delta=1.0)
        self.assertAlmostEqual(float(raw[10, 1, 1]),
                               nominal_10 + 2048.0, delta=1.0)
        self.assertAlmostEqual(float(raw[11, 1, 1]), nominal_11, delta=1.0)


class UtrRateClosureTestCase(lsst.utils.tests.TestCase):

    def testRecoversTrueRateOnCleanRamp(self):
        """utrRate(clean linear ramp) should recover params.rate within noise."""
        params = sim.RampParams(nReads=45, H=16, W=16, rate=87.0,
                                readNoise=5.0, poisson=False)
        cube = sim.makeRamp(params, rng=0)
        rate = sim.utrRate(cube)
        self.assertEqual(rate.shape, (16, 16))
        # Expected precision: LSQ slope std = readNoise * sqrt(12 / (N*(N-1)*(N+1)))
        # ≈ 5 * sqrt(12 / (45 * 44 * 46)) ≈ 0.12 ADU/read. So 3-sigma ~ 0.4.
        self.assertAlmostEqual(float(rate.mean()), 87.0, delta=0.5)
        self.assertLess(float(rate.std()), 0.5)

    def testRecoversRateWithPoissonNoise(self):
        params = sim.RampParams(nReads=45, H=16, W=16, rate=100.0,
                                readNoise=5.0, poisson=True)
        cube = sim.makeRamp(params, rng=0)
        rate = sim.utrRate(cube)
        # With Poisson, the variance grows along the ramp; the unweighted
        # LSQ slope is slightly less optimal but still close.
        self.assertAlmostEqual(float(rate.mean()), 100.0, delta=1.0)

    def testReadMask1DExcludesReads(self):
        """Excluding some reads should still give the right slope."""
        params = sim.RampParams(nReads=30, H=4, W=4, rate=50.0,
                                readNoise=0.0, poisson=False)
        cube = sim.makeRamp(params, rng=0)
        # Inject a CR at read 10 that pollutes reads 10..29.
        cube[10:, :, :] += 5000.0
        # All-reads fit will be wrong; masking out reads 10..29 should recover
        # the original slope.
        rate_all = sim.utrRate(cube)
        mask = np.ones(30, dtype=bool)
        mask[10:] = False
        rate_masked = sim.utrRate(cube, readMask=mask)
        self.assertNotAlmostEqual(float(rate_all.mean()), 50.0, delta=10.0)
        self.assertAlmostEqual(float(rate_masked.mean()), 50.0, delta=0.5)

    def testReadMask3DPerPixel(self):
        """Per-pixel read masks should work."""
        params = sim.RampParams(nReads=30, H=4, W=4, rate=50.0,
                                readNoise=0.0, poisson=False)
        cube = sim.makeRamp(params, rng=0)
        # CR at one pixel only.
        cube[15:, 2, 3] += 3000.0
        # Per-pixel mask: exclude reads 15.. for that one pixel; keep all for others.
        mask3D = np.ones((30, 4, 4), dtype=bool)
        mask3D[15:, 2, 3] = False
        rate = sim.utrRate(cube, readMask=mask3D)
        self.assertAlmostEqual(float(rate[2, 3]), 50.0, delta=0.5)
        self.assertAlmostEqual(float(rate[0, 0]), 50.0, delta=0.5)

    def testRoundTripRateThroughNonlinearityAndInverse(self):
        """sim(rate) -> rawRamp -> nl.inverse -> utrRate should recover rate."""
        params = sim.RampParams(nReads=45, H=8, W=8, rate=200.0,
                                readNoise=5.0, poisson=True)
        nl = sim.Nonlinearity(alpha=0.05, qMax=60000.0)
        _, linearized = sim.makeRawAndLinearRamps(params, nl, rng=0)
        rate = sim.utrRate(linearized)
        # Closure: recovered rate should equal the true rate within noise.
        self.assertAlmostEqual(float(rate.mean()), 200.0, delta=2.0)


class IterativeUtrDetectAndRepairTestCase(lsst.utils.tests.TestCase):
    """Tests for cr.iterativeUtrDetectAndRepair against simulated truth."""

    def _makeCube(self, params, crs=None, glitches=None, rng=0):
        # Note: the iterative detector works on the linearized cube. We
        # inject CRs in true space and glitches in measured space (the
        # standard sim chain) and feed the linearized cube into detect.
        nl = sim.Nonlinearity(alpha=0.0)  # tests focus on detection, not lin model
        raw = sim.makeRawRamp(params, nonlinearity=nl, crs=crs,
                              asicGlitches=glitches, rng=rng)
        return raw  # alpha=0 → raw == linearized

    def testBrightCRs(self):
        params = sim.RampParams(nReads=30, H=16, W=16, rate=50.0,
                                readNoise=5.0, poisson=True)
        truth = [sim.CR(y=3, x=4, read=10, amount=500.0),
                 sim.CR(y=8, x=8, read=20, amount=1500.0)]
        cube = self._makeCube(params, crs=truth, rng=0)
        good = np.ones((16, 16), dtype=bool)
        result = cr.iterativeUtrDetectAndRepair(
            cube, goodPixelMask=good, nSigma=5.0,
        )
        self.assertEqual(result.nGlitchPairs, 0)
        self.assertGreaterEqual(result.nCRs, len(truth))
        # Each truth CR should have a flag at the corresponding (delta_index, y, x).
        for c in truth:
            self.assertTrue(result.crFlagMask[c.read - 1, c.y, c.x],
                            f"missed CR at {(c.read, c.y, c.x)}")

    def testHighAmplitudeGlitches(self):
        params = sim.RampParams(nReads=30, H=16, W=16, rate=50.0,
                                readNoise=5.0, poisson=True)
        truth = [sim.AsicGlitch(y=2, x=2, read=10, amount=+2048.0),
                 sim.AsicGlitch(y=5, x=5, read=15, amount=-1024.0)]
        cube = self._makeCube(params, glitches=truth, rng=0)
        good = np.ones((16, 16), dtype=bool)
        glitchMask = np.ones((16, 16), dtype=bool)  # opt in across the cube
        result = cr.iterativeUtrDetectAndRepair(
            cube, goodPixelMask=good, glitchPixelMask=glitchMask,
            nSigma=5.0,
        )
        # Glitch at read m → pair at delta indices (m-1, m). Both should flag.
        for g in truth:
            self.assertTrue(result.glitchFlagMask[g.read - 1, g.y, g.x],
                            f"missed glitch at delta {g.read - 1}, ({g.y}, {g.x})")
            self.assertTrue(result.glitchFlagMask[g.read, g.y, g.x],
                            f"missed glitch at delta {g.read}, ({g.y}, {g.x})")
        # No CRs in this scenario.
        self.assertEqual(result.nCRs, 0)

    def testGlitchDetectionOffByDefault(self):
        """With glitchPixelMask=None (default), glitches are not detected."""
        params = sim.RampParams(nReads=30, H=8, W=8, rate=50.0,
                                readNoise=5.0, poisson=True)
        truth = [sim.AsicGlitch(y=3, x=3, read=10, amount=+2048.0)]
        cube = self._makeCube(params, glitches=truth, rng=0)
        good = np.ones((8, 8), dtype=bool)
        result = cr.iterativeUtrDetectAndRepair(
            cube, goodPixelMask=good, nSigma=5.0,
        )
        self.assertEqual(result.nGlitchPairs, 0,
                         "no glitches should be detected without glitchPixelMask")
        # The +half of the glitch (positive residual) gets classified as a CR
        # because pair detection is disabled. The negative half is flagged
        # but neither CR nor glitch.
        self.assertGreaterEqual(result.nCRs, 1,
                                "positive glitch half should be misclassified as CR")

    def testGlitchDetectionRestrictedToChannel(self):
        """Pair detection runs only where glitchPixelMask is True."""
        # Two identical glitches: one in the "good" channel (x=5), one outside.
        params = sim.RampParams(nReads=30, H=4, W=12, rate=50.0,
                                readNoise=5.0, poisson=True)
        truth = [
            sim.AsicGlitch(y=2, x=5, read=12, amount=+2048.0),   # inside channel
            sim.AsicGlitch(y=2, x=10, read=15, amount=+2048.0),  # outside
        ]
        cube = self._makeCube(params, glitches=truth, rng=0)
        good = np.ones((4, 12), dtype=bool)
        # Mask covers x in [4, 8) — channel for glitch detection.
        glitchMask = np.zeros((4, 12), dtype=bool)
        glitchMask[:, 4:8] = True
        result = cr.iterativeUtrDetectAndRepair(
            cube, goodPixelMask=good, glitchPixelMask=glitchMask,
            nSigma=5.0,
        )
        # In-channel glitch: detected as a pair.
        self.assertTrue(result.glitchFlagMask[11, 2, 5],
                        "in-channel glitch first half not flagged")
        self.assertTrue(result.glitchFlagMask[12, 2, 5],
                        "in-channel glitch second half not flagged")
        # Out-of-channel glitch: NOT detected as a pair.
        self.assertFalse(result.glitchFlagMask[14, 2, 10],
                         "out-of-channel glitch wrongly tagged as pair")
        self.assertFalse(result.glitchFlagMask[15, 2, 10],
                         "out-of-channel glitch wrongly tagged as pair")

    def testLowAmplitudeGlitchesProbeNoiseFloor(self):
        """N=3,4,5 glitches (8/16/32 ADU) at readNoise=5 → threshold ~35 ADU.

        Expectation: amplitudes well above the threshold are caught;
        amplitudes at or below it are not — that's the design limit.
        """
        params = sim.RampParams(nReads=30, H=24, W=8, rate=50.0,
                                readNoise=5.0, poisson=False)
        glitches = []
        # Three glitches per amplitude, at distinct pixels. Skip read 0 and
        # the last read (atEnd-style; glitches must be interior).
        amps_by_bit = {3: 8.0, 4: 16.0, 5: 32.0, 6: 64.0, 8: 256.0, 10: 1024.0}
        y = 0
        for bit, amp in amps_by_bit.items():
            for x in range(3):
                glitches.append(sim.AsicGlitch(y=y, x=x, read=10 + 2 * x,
                                                 amount=+amp))
            y += 1

        cube = self._makeCube(params, glitches=glitches, rng=42)
        good = np.ones((24, 8), dtype=bool)
        glitchMask = np.ones((24, 8), dtype=bool)
        # With poisson=False the only noise source is read noise; per-delta
        # std ≈ sqrt(2)·readNoise = sqrt(2)·5 ≈ 7 ADU. IQR sigma ≈ 7, so
        # threshold = 5·7 = 35 ADU. Bits with amp > 35 detect, amp ≤ 35 don't.
        result = cr.iterativeUtrDetectAndRepair(
            cube, goodPixelMask=good, glitchPixelMask=glitchMask,
            sigmaFloorADU=0.0,  # disable floor for this experiment
            nSigma=5.0,
        )

        detected_by_bit = {}
        for bit, amp in amps_by_bit.items():
            y_for_bit = list(amps_by_bit.keys()).index(bit)
            # Each amp has 3 glitches at row y_for_bit
            hits = 0
            for x in range(3):
                read = 10 + 2 * x
                if result.glitchFlagMask[read - 1, y_for_bit, x]:
                    hits += 1
            detected_by_bit[bit] = hits

        # Above the noise wall: at least one of three glitches detected
        # (allow some randomness from per-pixel rate jitter).
        self.assertGreaterEqual(detected_by_bit[10], 3,
                                f"bit-10 (1024 ADU) hits: {detected_by_bit[10]}")
        self.assertGreaterEqual(detected_by_bit[8], 3,
                                f"bit-8 (256 ADU) hits: {detected_by_bit[8]}")
        self.assertGreaterEqual(detected_by_bit[6], 2,
                                f"bit-6 (64 ADU) hits: {detected_by_bit[6]}")
        # Below the noise wall (32, 16, 8 ADU at threshold ~50): hits should
        # be small (mostly missed). Don't require zero — chance pickups OK.
        self.assertLessEqual(detected_by_bit[5], 2,
                             f"bit-5 (32 ADU) over-detected: {detected_by_bit[5]}")
        self.assertLessEqual(detected_by_bit[4], 1,
                             f"bit-4 (16 ADU) over-detected: {detected_by_bit[4]}")
        self.assertLessEqual(detected_by_bit[3], 1,
                             f"bit-3 (8 ADU) over-detected: {detected_by_bit[3]}")

    def testRepairRecoversCleanRamp(self):
        """After repair, utrRate of the cube should match the true rate."""
        params = sim.RampParams(nReads=40, H=8, W=8, rate=80.0,
                                readNoise=4.0, poisson=False)
        crs = [sim.CR(y=2, x=3, read=15, amount=1500.0)]
        glitches = [sim.AsicGlitch(y=5, x=5, read=20, amount=+2048.0)]
        cube = self._makeCube(params, crs=crs, glitches=glitches, rng=0)
        good = np.ones((8, 8), dtype=bool)
        glitchMask = np.ones((8, 8), dtype=bool)
        result = cr.iterativeUtrDetectAndRepair(
            cube, goodPixelMask=good, glitchPixelMask=glitchMask,
            nSigma=5.0,
        )
        self.assertEqual(result.nCRs, 1)
        self.assertEqual(result.nGlitchPairs, 1)
        # Post-repair rate should match params.rate.
        rate_after = sim.utrRate(cube)
        self.assertAlmostEqual(float(rate_after.mean()), params.rate, delta=0.5)
        # Specifically the CR pixel and the glitch pixel should also match.
        self.assertAlmostEqual(float(rate_after[2, 3]), params.rate, delta=1.0)
        self.assertAlmostEqual(float(rate_after[5, 5]), params.rate, delta=1.0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
