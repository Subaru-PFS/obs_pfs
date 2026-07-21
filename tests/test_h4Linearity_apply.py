"""End-to-end integration tests for ``h4Linearity.apply`` with an
asymmetric, position-encoded ramp.

The polynomial models themselves are exercised in ``test_polynomial.py``;
this file's job is to confirm the cube-level orchestration in ``apply``
preserves axis order and pixel identity. Each per-pixel Chebyshev
coefficient depends on ``(y, x)`` so an axis swap would yield a
visibly wrong output value at any pixel, not just at the corners.

Conventions verified here:
- Input ``Ramp.reads`` is ``(H, W, N)`` with the time axis last;
  output ``cumulativeLinear`` is the same shape.
- Output ``badPixelMask`` is ``(H, W)`` with axes matching the spatial
  plane of ``reads``.
- ``BELOW_VALID_RANGE`` / ``ABOVE_VALID_RANGE`` bits fire on the
  per-pixel ``fitMin`` / ``fitMax`` interval, not a global one.
- Bad pixels (``correction.badPixelMask != 0`` or ``ramp.validMask``)
  pass their input value through unchanged.
"""
import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs.h4Linearity.apply import apply
from lsst.obs.pfs.h4Linearity.models.polynomial import PolynomialModel
from lsst.obs.pfs.h4Linearity.types import (
    ABOVE_VALID_RANGE,
    BELOW_VALID_RANGE,
    Diagnostics,
    LinearityCorrection,
    MASKED_BY_INPUT,
    Ramp,
)


N_READS = 11
H = 5
W = 7
FIT_MIN_VALUE = 0.0
FIT_MAX_VALUE = 100.0


# Per-pixel linear-in-Chebyshev-domain slope. With c1[y, x] = y * 10 + x
# and c0 = 0, the linearized t = c1[y, x] * cheb_x; mismatched axes
# produce a different value at every pixel.
def _expectedC1(H, W):
    y = np.arange(H, dtype=np.float32)[:, None]
    x = np.arange(W, dtype=np.float32)[None, :]
    return y * 10.0 + x


def _makeCorrection(H, W, fitMin=FIT_MIN_VALUE, fitMax=FIT_MAX_VALUE,
                    badPixelMask=None):
    """Build a degree-1 polynomial LinearityCorrection.

    coefs[0] = 0 (no offset), coefs[1] = y * 10 + x (per-pixel slope).
    fitMin/fitMax are uniform unless overridden.
    """
    coefs = np.zeros((2, H, W), dtype=np.float32)
    coefs[1] = _expectedC1(H, W)
    if badPixelMask is None:
        badPixelMask = np.zeros((H, W), dtype=np.uint8)
    diag = Diagnostics(
        residualRms=np.zeros((H, W), dtype=np.float32),
        maxAbsResidual=np.zeros((H, W), dtype=np.float32),
        nPointsUsed=np.full((H, W), N_READS, dtype=np.int32),
        monotonic=np.ones((H, W), dtype=bool),
        conditionNumber=np.ones((H, W), dtype=np.float32),
        summary={},
    )
    return LinearityCorrection(
        model=PolynomialModel(order=1),
        coefficients=coefs,
        fitMin=np.full((H, W), fitMin, dtype=np.float32),
        fitMax=np.full((H, W), fitMax, dtype=np.float32),
        badPixelMask=badPixelMask,
        diagnostics=diag,
    )


RAMP_RATE = 5.0
# cumulative[k] = (k+1)*RATE; chosen so even N=11 stays inside the
# default [FIT_MIN_VALUE, FIT_MAX_VALUE] window.


def _makeRamp(N, H, W, *, fillValue=None):
    """Build a ``(H, W, N)`` cumulative ramp.

    ``fillValue=None`` (default): cumulative[y, x, k] = (k+1) * RAMP_RATE
    (uniform across pixels, varies along reads). Easy to predict in the
    Chebyshev-rescaled domain.
    """
    if fillValue is None:
        k = np.arange(1, N + 1, dtype=np.float32) * RAMP_RATE
        return np.broadcast_to(k, (H, W, N)).astype(np.float32, copy=True)
    return np.full((H, W, N), float(fillValue), dtype=np.float32)


def _expectedLinearized(reads, c1, fitMin, fitMax):
    """Reference: t = c1 * cheb_x with cheb_x = 2*(m-fitMin)/(fitMax-fitMin) - 1.

    ``reads`` is ``(H, W, N)``; ``c1`` is ``(H, W)``; ``fitMin`` / ``fitMax``
    are scalars in this test. Broadcasts ``c1[..., None]`` along the
    time axis.
    """
    chebX = 2.0 * (reads - fitMin) / (fitMax - fitMin) - 1.0
    return (c1[..., None] * chebX).astype(np.float32)


class ApplyAxisOrderTestCase(lsst.utils.tests.TestCase):
    """The cube comes back in the same axis order it went in, with each
    pixel linearized through *its own* coefficients."""

    def testOutputShapeMatchesInput(self):
        correction = _makeCorrection(H, W)
        ramp = Ramp(reads=_makeRamp(N_READS, H, W))
        result = apply(correction, ramp)
        self.assertEqual(result.cumulativeLinear.shape, (H, W, N_READS))
        self.assertEqual(result.badPixelMask.shape, (H, W))

    def testValuesMatchPerPixelCoefficients(self):
        correction = _makeCorrection(H, W)
        readsOrig = _makeRamp(N_READS, H, W)
        ramp = Ramp(reads=readsOrig.copy())  # apply mutates ramp.reads
        result = apply(correction, ramp)
        expected = _expectedLinearized(
            readsOrig, _expectedC1(H, W), FIT_MIN_VALUE, FIT_MAX_VALUE
        )
        np.testing.assert_allclose(
            result.cumulativeLinear, expected, atol=1e-4,
            err_msg="apply() must linearize each pixel through its own coefficients",
        )

    def testPositionAsymmetricSpotChecks(self):
        # Pin the result at three specific asymmetric (y, x, k) so a
        # transposition of (y, x) — silent under any H==W fixture — is
        # caught by direct value comparison.
        correction = _makeCorrection(H, W)
        ramp = Ramp(reads=_makeRamp(N_READS, H, W))
        result = apply(correction, ramp)
        # cheb_x at read k for m=(k+1)*RAMP_RATE: 2*m/100 - 1 = (k+1)*0.1 - 1.
        # t at (y, x) = c1[y, x] * cheb_x = (y*10 + x) * ((k+1)*0.1 - 1).
        for k, y, x in [(0, 2, 3), (5, 0, 6), (10, 4, 1)]:
            chebX = (k + 1) * 0.1 - 1.0
            expected = (y * 10.0 + x) * chebX
            self.assertAlmostEqual(
                float(result.cumulativeLinear[y, x, k]),
                float(expected),
                places=4,
                msg=f"apply mismatch at (y={y}, x={x}, k={k})",
            )


class ApplyBadPixelTestCase(lsst.utils.tests.TestCase):
    """Bad pixels (fit-time mask, caller validMask, both) pass input
    through unchanged and the merged mask flags them in the output."""

    def testFitTimeBadPixelsPassThrough(self):
        badMask = np.zeros((H, W), dtype=np.uint8)
        badMask[2, 3] = 0x04  # arbitrary fit-time bad bit
        correction = _makeCorrection(H, W, badPixelMask=badMask)
        readsOrig = _makeRamp(N_READS, H, W)
        ramp = Ramp(reads=readsOrig.copy())
        result = apply(correction, ramp)
        # Bad pixel passed through unchanged on every read.
        np.testing.assert_array_equal(
            result.cumulativeLinear[2, 3, :],
            readsOrig[2, 3, :],
            err_msg="fit-time bad pixel must pass through unchanged",
        )
        # The fit-time bit survives into the output mask.
        self.assertEqual(int(result.badPixelMask[2, 3]) & 0x04, 0x04)
        # And other pixels are NOT bad.
        self.assertEqual(int(result.badPixelMask[0, 0]), 0)
        self.assertEqual(int(result.badPixelMask[4, 6]), 0)

    def testValidMaskMergedAsMaskedByInput(self):
        correction = _makeCorrection(H, W)
        validMask = np.zeros((H, W), dtype=np.uint8)
        validMask[1, 5] = 1  # caller-supplied defect
        readsOrig = _makeRamp(N_READS, H, W)
        ramp = Ramp(reads=readsOrig.copy(), validMask=validMask)
        result = apply(correction, ramp)
        np.testing.assert_array_equal(
            result.cumulativeLinear[1, 5, :], readsOrig[1, 5, :],
            err_msg="caller-flagged defect must pass through unchanged",
        )
        self.assertEqual(
            int(result.badPixelMask[1, 5]) & MASKED_BY_INPUT,
            MASKED_BY_INPUT,
        )


class ApplyRangeFlagsTestCase(lsst.utils.tests.TestCase):
    """Per-pixel fitMin/fitMax range checks stamp the right bits."""

    def testAboveRangeFlagged(self):
        # Push pixel (3, 4) over fitMax via a final read above 100.
        correction = _makeCorrection(H, W)
        reads = _makeRamp(N_READS, H, W)
        reads[3, 4, -1] = FIT_MAX_VALUE + 1.0
        result = apply(correction, Ramp(reads=reads.copy()))
        self.assertEqual(
            int(result.badPixelMask[3, 4]) & ABOVE_VALID_RANGE,
            ABOVE_VALID_RANGE,
        )
        # Neighbouring pixels stay un-flagged.
        self.assertEqual(int(result.badPixelMask[3, 3]) & ABOVE_VALID_RANGE, 0)
        self.assertEqual(int(result.badPixelMask[4, 4]) & ABOVE_VALID_RANGE, 0)

    def testBelowRangeFlagged(self):
        correction = _makeCorrection(H, W)
        reads = _makeRamp(N_READS, H, W)
        reads[1, 2, 0] = FIT_MIN_VALUE - 1.0
        result = apply(correction, Ramp(reads=reads.copy()))
        self.assertEqual(
            int(result.badPixelMask[1, 2]) & BELOW_VALID_RANGE,
            BELOW_VALID_RANGE,
        )

    def testPerPixelFitRangeNotGlobal(self):
        # Pixel (2, 5) gets a tighter range; reads that are within the
        # global range but outside this pixel's range should still flag.
        correction = _makeCorrection(H, W)
        # Mutate via construction since LinearityCorrection is frozen.
        fitMax = correction.fitMax.copy()
        fitMax[2, 5] = 50.0
        correction = LinearityCorrection(
            model=correction.model,
            coefficients=correction.coefficients,
            fitMin=correction.fitMin,
            fitMax=fitMax,
            badPixelMask=correction.badPixelMask,
            diagnostics=correction.diagnostics,
        )
        # Default ramp has reads up to N_READS * RAMP_RATE = 55, comfortably
        # above the per-pixel fitMax of 50 at (2, 5).
        reads = _makeRamp(N_READS, H, W)
        result = apply(correction, Ramp(reads=reads.copy()))
        self.assertEqual(
            int(result.badPixelMask[2, 5]) & ABOVE_VALID_RANGE,
            ABOVE_VALID_RANGE,
        )
        # A neighbour with the (default) fitMax=100 stays clean.
        self.assertEqual(
            int(result.badPixelMask[2, 4]) & ABOVE_VALID_RANGE, 0
        )


class ApplyInputValidationTestCase(lsst.utils.tests.TestCase):
    """Shape / dimensionality errors surface as ValueError."""

    def testRejectsNon3D(self):
        correction = _makeCorrection(H, W)
        with self.assertRaises(ValueError):
            apply(correction, Ramp(reads=np.zeros((H, W), dtype=np.float32)))

    def testRejectsZeroReads(self):
        correction = _makeCorrection(H, W)
        with self.assertRaises(ValueError):
            apply(correction, Ramp(reads=np.zeros((H, W, 0), dtype=np.float32)))

    def testRejectsHWMismatch(self):
        correction = _makeCorrection(H, W)
        # A (H+1, W) ramp spatial plane must be rejected, not silently broadcast.
        with self.assertRaises(ValueError):
            apply(
                correction,
                Ramp(reads=np.zeros((H + 1, W, N_READS), dtype=np.float32)),
            )


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
