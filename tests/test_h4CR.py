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


def _injectCR(flux, y, x, k, amount):
    """Add a CR contribution at delta index k: cumulative reads >= k+1 jump by amount."""
    flux[k + 1:, y, x] += amount


class CRDetectAndRepairTestCase(lsst.utils.tests.TestCase):

    def testNoCRReturnsClean(self):
        flux = _flatRamp()
        good = np.ones(flux.shape[1:], dtype=bool)
        before = flux.copy()

        result = cr.detectAndRepair(flux, goodPixelMask=good)

        self.assertEqual(result.nFlagged, 0)
        self.assertFalse(result.flagMask.any())
        np.testing.assert_array_equal(flux, before)

    def testSingleCRDetectedAndRepaired(self):
        flux = _flatRamp(nReads=16, rate=10.0)
        _injectCR(flux, y=3, x=4, k=7, amount=200.0)
        good = np.ones(flux.shape[1:], dtype=bool)

        result = cr.detectAndRepair(flux, goodPixelMask=good)

        self.assertEqual(result.nFlagged, 1)
        self.assertTrue(result.flagMask[3, 4])
        self.assertEqual(int(result.crReadIdx[3, 4]), 7)
        self.assertAlmostEqual(float(result.excess[3, 4]), 200.0, delta=1.0)
        # Repair should restore the ramp to (approximately) the clean line.
        np.testing.assert_allclose(flux, _flatRamp(nReads=16, rate=10.0), atol=1.0)

    def testFaintPixelCRStillDetected(self):
        """A small CR on a near-zero-rate pixel must still trip the threshold."""
        flux = _flatRamp(nReads=20, rate=1.0)
        _injectCR(flux, y=2, x=2, k=10, amount=40.0)  # well above the 15-ADU excess floor
        good = np.ones(flux.shape[1:], dtype=bool)

        result = cr.detectAndRepair(flux, goodPixelMask=good)

        self.assertEqual(result.nFlagged, 1)
        self.assertTrue(result.flagMask[2, 2])

    def testBelowExcessFloorNotFlagged(self):
        """Spikes smaller than the absolute excess floor stay unflagged."""
        flux = _flatRamp(nReads=20, rate=1.0)
        _injectCR(flux, y=2, x=2, k=10, amount=10.0)  # below default excessFloorADU=15
        good = np.ones(flux.shape[1:], dtype=bool)

        result = cr.detectAndRepair(flux, goodPixelMask=good)

        self.assertEqual(result.nFlagged, 0)

    def testMaskedPixelSkipped(self):
        """A pixel marked not-good is never flagged, even with a clear CR."""
        flux = _flatRamp(nReads=16, rate=10.0)
        _injectCR(flux, y=3, x=4, k=7, amount=200.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        good[3, 4] = False
        before = flux.copy()

        result = cr.detectAndRepair(flux, goodPixelMask=good)

        self.assertEqual(result.nFlagged, 0)
        # And the flux at that pixel is left untouched.
        np.testing.assert_array_equal(flux[:, 3, 4], before[:, 3, 4])

    def testShortRampSkipped(self):
        """Ramps shorter than minReads are no-ops (no detection, no repair)."""
        flux = _flatRamp(nReads=4, rate=10.0)
        _injectCR(flux, y=3, x=4, k=2, amount=200.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        before = flux.copy()

        result = cr.detectAndRepair(flux, goodPixelMask=good, minReads=8)

        self.assertEqual(result.nFlagged, 0)
        np.testing.assert_array_equal(flux, before)

    def testMultipleCRsIndependent(self):
        flux = _flatRamp(nReads=20, rate=10.0)
        crLocs = [(0, 0, 5, 100.0), (3, 4, 10, 250.0), (7, 7, 15, 80.0)]
        for y, x, k, amount in crLocs:
            _injectCR(flux, y, x, k, amount)
        good = np.ones(flux.shape[1:], dtype=bool)

        result = cr.detectAndRepair(flux, goodPixelMask=good)

        self.assertEqual(result.nFlagged, len(crLocs))
        for y, x, k, _ in crLocs:
            self.assertTrue(result.flagMask[y, x], f"missing flag at ({y},{x})")
            self.assertEqual(int(result.crReadIdx[y, x]), k)
        np.testing.assert_allclose(flux, _flatRamp(nReads=20, rate=10.0), atol=1.0)

    def testRepairOnlyAffectsReadsAtOrAfterCR(self):
        """Repair must subtract only from cumulative reads >= k_cr+1."""
        flux = _flatRamp(nReads=12, rate=5.0)
        _injectCR(flux, y=1, x=1, k=4, amount=300.0)
        before = flux.copy()
        good = np.ones(flux.shape[1:], dtype=bool)

        result = cr.detectAndRepair(flux, goodPixelMask=good)

        self.assertEqual(result.nFlagged, 1)
        # Reads 0..k_cr are untouched; reads k_cr+1..N-1 are reduced by ~excess.
        np.testing.assert_array_equal(flux[: 4 + 1, 1, 1], before[: 4 + 1, 1, 1])
        diffs = before[5:, 1, 1] - flux[5:, 1, 1]
        np.testing.assert_allclose(diffs, np.full(diffs.shape, 300.0), atol=1.0)

    def testInPlaceModification(self):
        """The flux argument is the array that gets mutated; no copy is returned."""
        flux = _flatRamp(nReads=16, rate=10.0)
        _injectCR(flux, y=2, x=3, k=6, amount=150.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        flux_id = id(flux)

        cr.detectAndRepair(flux, goodPixelMask=good)

        # Same object, modified in place.
        self.assertEqual(id(flux), flux_id)
        np.testing.assert_allclose(flux, _flatRamp(nReads=16, rate=10.0), atol=1.0)

    def testGoodPixelMaskShapeChecked(self):
        flux = _flatRamp()
        with self.assertRaises(ValueError):
            cr.detectAndRepair(flux, goodPixelMask=np.ones((4, 4), dtype=bool))

    def testFluxDimsChecked(self):
        with self.assertRaises(ValueError):
            cr.detectAndRepair(
                np.zeros((4, 4), dtype=np.float32),
                goodPixelMask=np.ones((4, 4), dtype=bool),
            )

    def testWithoutDiagnosticsReturnsNoneFields(self):
        flux = _flatRamp(nReads=16, rate=10.0)
        _injectCR(flux, y=3, x=4, k=7, amount=200.0)
        good = np.ones(flux.shape[1:], dtype=bool)

        result = cr.detectAndRepair(flux, goodPixelMask=good)

        self.assertIsNone(result.candidateMask)
        self.assertIsNone(result.medianDelta)
        self.assertIsNone(result.sigma)

    def testReturnDiagnosticsPopulatesCandidateOnlyArrays(self):
        flux = _flatRamp(nReads=16, rate=10.0)
        _injectCR(flux, y=3, x=4, k=7, amount=200.0)
        good = np.ones(flux.shape[1:], dtype=bool)

        result = cr.detectAndRepair(flux, goodPixelMask=good, returnDiagnostics=True)

        self.assertIsNotNone(result.candidateMask)
        self.assertIsNotNone(result.medianDelta)
        self.assertIsNotNone(result.sigma)
        # The CR pixel passed pre-screen and was flagged.
        self.assertTrue(result.candidateMask[3, 4])
        self.assertTrue(result.flagMask[3, 4])
        # Non-CR pixels (max delta = rate = 10 < excessFloorADU=15) are not candidates.
        self.assertFalse(result.candidateMask[0, 0])
        # Per-candidate fields populated where candidate is True.
        self.assertAlmostEqual(float(result.medianDelta[3, 4]), 10.0, delta=1.0)
        # And zero elsewhere.
        self.assertEqual(float(result.medianDelta[0, 0]), 0.0)
        self.assertEqual(float(result.sigma[0, 0]), 0.0)

    def testPreScreenSkipsBelowFloorPixels(self):
        """Flat ramp with rate above excessFloorADU should produce candidates but no flags;
        flat ramp at rate=1 should produce no candidates at all."""
        # Rate above floor: every pixel is a candidate, but no CR was injected.
        flux = _flatRamp(nReads=16, rate=20.0)
        good = np.ones(flux.shape[1:], dtype=bool)
        result = cr.detectAndRepair(flux, goodPixelMask=good, returnDiagnostics=True)
        self.assertEqual(result.nFlagged, 0)
        # Some pixels passed the pre-screen (max_delta == 20 > 15).
        self.assertGreater(int(result.candidateMask.sum()), 0)

        # Rate below floor: no pixel passes pre-screen.
        flux = _flatRamp(nReads=16, rate=5.0)
        result = cr.detectAndRepair(flux, goodPixelMask=good, returnDiagnostics=True)
        self.assertEqual(int(result.candidateMask.sum()), 0)
        self.assertEqual(result.nFlagged, 0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
