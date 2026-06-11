"""Integration tests for ``validate.processRamp`` on an asymmetric ramp.

Uses the ``reuseExposure`` path (pre-built exposure + cube, butler
not touched) so we can drive the CR-detection portion of
``processRamp`` with deterministic synthetic inputs and verify:

- The cube returned has shape ``(H, W, N)`` (axis order preserved
  through the in-helper diff -> CR detector -> cumsum reconstruction).
- The exposure mask is stamped at ``(y, x)`` for CR / ASIC_GLITCH;
  axes 0 and 1 of the cube map onto ``mask.array[y, x]``, not
  ``mask.array[x, y]``.
- ``crResult.crFlagMask`` and ``glitchFlagMask`` have shape
  ``(H, W, N-1)`` with hits at the expected delta indices.
- The pre-seeded ``intermediates`` opt-in captures ``'crCorrected'``
  only when requested.

Asymmetric shape (N=20, H=5, W=7) is used throughout so any axis
permutation surfaces as a shape error rather than a silent transposed
result.
"""
import unittest

import numpy as np

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.utils.tests

from lsst.obs.pfs.h4Linearity import validate as pfsValidate


N_READS = 20
H = 5
W = 7
RATE = 10.0  # cumulative[y, x, k] = (k+1) * RATE


def _flatCube(nReads=N_READS, h=H, w=W, rate=RATE, dtype=np.float32):
    """``cumulative[y, x, k] = (k+1) * rate`` over a flat (H, W) plane,
    returned in the project's ``(H, W, N)`` time-axis-last layout."""
    k = np.arange(1, nReads + 1, dtype=np.float64) * rate
    return np.broadcast_to(k, (h, w, nReads)).astype(dtype, copy=True)


def _makeExposure(h=H, w=W):
    """Build an ExposureF with mask + image planes at (H, W)."""
    # Extent2I is (width, height); mask.array.shape becomes (H, W).
    return afwImage.ExposureF(geom.Extent2I(w, h))


class ProcessRampAxisOrderTestCase(lsst.utils.tests.TestCase):
    """Verify processRamp preserves the (H, W, N) convention and stamps
    masks at the right ``(y, x)``."""

    def testCubeShapePreservedNoDefects(self):
        # No injected defects: the cube comes back the same shape.
        cube = _flatCube()
        exp = _makeExposure()
        _, cubeOut, crResult = pfsValidate.processRamp(
            butler=None, dataId={}, cam="n3",
            exposure=exp, cube=cube.copy(),
            doLinearize=False, doCR=True, repairCR=True,
        )
        self.assertEqual(cubeOut.shape, (H, W, N_READS))
        self.assertEqual(crResult.crFlagMask.shape, (H, W, N_READS - 1))
        self.assertEqual(crResult.glitchFlagMask.shape, (H, W, N_READS - 1))
        self.assertEqual(int(crResult.nCRs), 0)
        self.assertEqual(int(crResult.nGlitchPairs), 0)

    def testCRPixelStampedAtCorrectYX(self):
        # Inject a CR at (y=3, x=4, read 10): cube[3, 4, 10:] += 200.
        # The single delta delta[3, 4, 9] = RATE + 200, far above
        # threshold; result.crFlagMask[3, 4, 9] should fire and the
        # exposure mask should get the CR bit at (3, 4).
        cube = _flatCube()
        cube[3, 4, 10:] += 200.0
        exp = _makeExposure()
        _, cubeOut, crResult = pfsValidate.processRamp(
            butler=None, dataId={}, cam="n3",
            exposure=exp, cube=cube.copy(),
            doLinearize=False, doCR=True, repairCR=True,
        )
        self.assertTrue(
            bool(crResult.crFlagMask[3, 4, 9]),
            "CR flag must land at delta index k-1 of the injected read",
        )
        # The exposure mask carries the CR bit at the spatial (y, x) =
        # (3, 4) — not at (4, 3) or any transposed location.
        crBit = exp.mask.getPlaneBitMask("CR")
        self.assertTrue(
            bool(exp.mask.array[3, 4] & crBit),
            "exposure.mask.array[y, x] must have CR bit set at (y=3, x=4)",
        )
        # And NO CR bit at the symmetric (4, 3) position — fixture
        # shape is (5, 7) so (4, 3) is valid; if axes were swapped the
        # mask would be at (4, 3) instead.
        self.assertFalse(
            bool(exp.mask.array[4, 3] & crBit),
            "no CR should be flagged at (y=4, x=3); axis (y, x) "
            "would be swapped if this fires",
        )

    def testRepairRestoresFluxAtRepairedPixel(self):
        # With repairCR=True, the cube returned has the CR contribution
        # subtracted; ``cube[3, 4, 10]`` should be ~11*RATE again.
        cube = _flatCube()
        cubeOrig = cube.copy()
        cube[3, 4, 10:] += 200.0
        exp = _makeExposure()
        _, cubeOut, _ = pfsValidate.processRamp(
            butler=None, dataId={}, cam="n3",
            exposure=exp, cube=cube,
            doLinearize=False, doCR=True, repairCR=True,
        )
        # The repaired ramp at (3, 4) matches the original flat ramp
        # to within ~0.5 ADU (the detector's iterative residual).
        np.testing.assert_allclose(
            cubeOut[3, 4, :], cubeOrig[3, 4, :], atol=2.0,
            err_msg="repaired ramp at (3, 4) must align with the flat baseline",
        )

    def testNoRepairLeavesCubeUntouched(self):
        # With repairCR=False, the cube isn't modified — the spike at
        # cube[3, 4, 10:] survives; only the mask gets stamped.
        cube = _flatCube()
        cube[3, 4, 10:] += 200.0
        cubeIn = cube.copy()
        exp = _makeExposure()
        _, cubeOut, crResult = pfsValidate.processRamp(
            butler=None, dataId={}, cam="n3",
            exposure=exp, cube=cube,
            doLinearize=False, doCR=True, repairCR=False,
        )
        np.testing.assert_array_equal(cubeOut, cubeIn)
        # But the flag still fires.
        self.assertTrue(bool(crResult.crFlagMask[3, 4, 9]))


class ProcessRampIntermediatesTestCase(lsst.utils.tests.TestCase):
    """The pre-seeded ``intermediates`` opt-in API: only requested keys
    are populated; unrelated keys remain absent."""

    def testOptInCrCorrectedOnly(self):
        cube = _flatCube()
        cube[3, 4, 10:] += 200.0
        exp = _makeExposure()
        intermediates = {"crCorrected": None}
        _, cubeOut, _ = pfsValidate.processRamp(
            butler=None, dataId={}, cam="n3",
            exposure=exp, cube=cube.copy(),
            doLinearize=False, doCR=True, repairCR=True,
            intermediates=intermediates,
        )
        self.assertIn("crCorrected", intermediates)
        self.assertIsNotNone(intermediates["crCorrected"])
        self.assertEqual(intermediates["crCorrected"].shape, (H, W, N_READS))
        # Sanity: the captured crCorrected matches cubeOut elementwise
        # (both are the post-repair cumulative cube).
        np.testing.assert_allclose(intermediates["crCorrected"], cubeOut)
        # Other keys did not appear.
        self.assertNotIn("raw", intermediates)
        self.assertNotIn("darkSubbed", intermediates)
        self.assertNotIn("linearized", intermediates)

    def testNoneIntermediatesCaptureNothing(self):
        cube = _flatCube()
        exp = _makeExposure()
        _, _, _ = pfsValidate.processRamp(
            butler=None, dataId={}, cam="n3",
            exposure=exp, cube=cube.copy(),
            doLinearize=False, doCR=True, repairCR=True,
            intermediates=None,
        )
        # Nothing to assert beyond "doesn't crash"; this exercises
        # the ``intermediates is None`` short-circuit in processRamp.


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
