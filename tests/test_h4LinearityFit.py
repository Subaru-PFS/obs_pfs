"""Tests for the top-level fit() function."""

from __future__ import annotations

import numpy as np

from lsst.obs.pfs.h4Linearity.fit import fit
from lsst.obs.pfs.h4Linearity.models import PolynomialModel
from lsst.obs.pfs.h4Linearity.types import MASKED_BY_INPUT, Ramp


def test_fitSingleRampRecoversTarget(smallSyntheticRamp):
    ramp, truth = smallSyntheticRamp
    # Disable HIGH_FIT_RESIDUAL flagging: the synthetic noise floor between
    # otherwise-identical pixels can trip the 5×median rule on a tiny fixture.
    correction = fit([ramp], badLinearityMedianMultiplier=None)
    assert correction.coefficients.shape == (5, 4, 5)
    assert (correction.badPixelMask == 0).all()
    # Evaluate at the fit points: map m → x, then evaluate.
    m = ramp.reads.astype(np.float32)
    denom = correction.fitMax - correction.fitMin
    denom = np.where(denom > 0, denom, 1.0)
    x = 2.0 * (m - correction.fitMin[None]) / denom[None] - 1.0
    # Upstream evaluate expects spatial axes leading: (H, W, ...) — bridge until
    # the PIPE2D-1843 transpose is finished in evaluateMonomial too.
    tPred = np.moveaxis(
        correction.model.evaluate(correction.coefficients, np.moveaxis(x, 0, -1)),
        -1, 0,
    )
    fitRate = float(np.median(ramp.reads[2] - ramp.reads[1]))
    N = ramp.reads.shape[0]
    expected = fitRate * np.arange(N, dtype=np.float32)
    expectedBroad = np.broadcast_to(expected[:, None, None], tPred.shape)
    np.testing.assert_allclose(tPred, expectedBroad, rtol=1e-3, atol=1.0)
    # Summary should carry percentiles.
    assert "residualRmsP50" in correction.diagnostics.summary
    assert "residualRmsP95" in correction.diagnostics.summary
    assert "residualRmsP99" in correction.diagnostics.summary


def test_fitTilingIsDeterministic(smallSyntheticRamp):
    """Fitting with different block sizes must yield identical coefficients
    (within float32 precision)."""
    ramp, _ = smallSyntheticRamp
    c1 = fit([ramp], blockSize=(4, 5), badLinearityMedianMultiplier=None).coefficients
    c2 = fit([ramp], blockSize=(2, 3), badLinearityMedianMultiplier=None).coefficients
    np.testing.assert_allclose(c1, c2, rtol=1e-5, atol=1e-5)


def test_fitPropagatesInputMask(tinyLinearRamp):
    ramp, _ = tinyLinearRamp
    mask = np.zeros(ramp.reads.shape[1:], dtype=np.uint8)
    mask[0, 0] = 1  # Mark pixel (0, 0) as invalid
    maskedRamp = Ramp(reads=ramp.reads, validMask=mask)
    correction = fit(
        [maskedRamp], model=PolynomialModel(order=1), badLinearityMedianMultiplier=None
    )
    assert correction.badPixelMask[0, 0] & MASKED_BY_INPUT
    assert correction.badPixelMask[0, 1] == 0


def test_fitMultipleRampsConcatenates():
    """Two ramps of different lengths combine per-pixel."""
    H, W = 3, 4
    # Pixel-linear: t = m for every pixel.
    # Ramp 1: 9 reads (1 implicit zero + 8), rate 100.
    # Ramp 2: 13 reads (1 implicit zero + 12), rate 200.
    rate1 = 100.0
    rate2 = 200.0
    N1, N2 = 9, 13
    reads1 = np.full((1, H, W), rate1, dtype=np.float32) * np.arange(N1, dtype=np.float32)[:, None, None]
    reads2 = np.full((1, H, W), rate2, dtype=np.float32) * np.arange(N2, dtype=np.float32)[:, None, None]
    correction = fit(
        [Ramp(reads=reads1), Ramp(reads=reads2)],
        model=PolynomialModel(order=2),
        badLinearityMedianMultiplier=None,
    )
    # Verify via evaluation: for ramp1, evaluate at its m values → should match targets.
    m1 = reads1
    denom = correction.fitMax - correction.fitMin
    denom = np.where(denom > 0, denom, 1.0)
    x1 = 2.0 * (m1 - correction.fitMin[None]) / denom[None] - 1.0
    tPred = np.moveaxis(
        correction.model.evaluate(correction.coefficients, np.moveaxis(x1, 0, -1)),
        -1, 0,
    )
    # Target for ramp1: rate1 * n (with target[0] = 0).
    expected1 = rate1 * np.arange(N1, dtype=np.float32)
    expected1Broad = np.broadcast_to(expected1[:, None, None], tPred.shape)
    np.testing.assert_allclose(tPred, expected1Broad, rtol=1e-4, atol=1e-2)
    # nPointsUsed should be N1 + N2 everywhere.
    assert (correction.diagnostics.nPointsUsed == N1 + N2).all()


def test_fitEmptyRampListRaises():
    import pytest
    with pytest.raises(ValueError):
        fit([])


# --- New behaviour: saturation-knee + HIGH_FIT_RESIDUAL flag ---------------

def _buildSaturatingRamp():
    """Build a 30-read 4x5 ramp where pixel (2, 2) saturates partway through.

    All other pixels ramp linearly at rate 1000 DN/read. Pixel (2, 2) ramps
    at the same rate through read 5, then its deltas drop to 1% of rate
    (200 → 10 DN) from read 6 onward — the classic 'faint saturating' case
    that the saturation-knee detector should catch from delta 0.
    """
    H, W, N = 4, 5, 30
    rate = 1000.0
    reads = (
        rate
        * np.arange(N, dtype=np.float32)[:, None, None]
        * np.ones((H, W), dtype=np.float32)
    )
    # Pixel (2, 2): keep first 5 deltas at rate, then collapse to 1% of rate.
    pr, pc = 2, 2
    for n in range(6, N):
        reads[n, pr, pc] = reads[n - 1, pr, pc] + 0.01 * rate
    return Ramp(reads=reads)


def test_saturationKneeFlagsSaturatingPixelByShrinkingFitWindow():
    ramp = _buildSaturatingRamp()
    # Knee 0.5 must catch pixel (2,2): once its delta drops to 1% of rate
    # (much less than 0.5 × refDelta), the read is masked.
    correction = fit(
        [ramp],
        model=PolynomialModel(order=2),
        saturationKnee=0.5,
        badLinearityMedianMultiplier=None,
    )
    nUsedSat = int(correction.diagnostics.nPointsUsed[2, 2])
    nUsedHealthy = int(correction.diagnostics.nPointsUsed[0, 0])
    # Saturator should drop the late-saturated reads; healthy pixel keeps everything.
    assert nUsedSat <= 7, f"expected <= 7 reads kept on saturator, got {nUsedSat}"
    assert nUsedHealthy == 30
    # Without the knee, the same ramp keeps every read on every pixel
    # (the existing deviationLimit-only path is disabled by default).
    correctionNoKnee = fit(
        [ramp],
        model=PolynomialModel(order=2),
        saturationKnee=None,
        badLinearityMedianMultiplier=None,
    )
    assert int(correctionNoKnee.diagnostics.nPointsUsed[2, 2]) == 30


def test_badLinearityMultiplierFlagsHighResidualPixel():
    from lsst.obs.pfs.h4Linearity.types import HIGH_FIT_RESIDUAL

    ramp = _buildSaturatingRamp()
    # Disable the saturation-knee so the saturating pixel survives into the
    # fit with a bad residual — that gives the residual-based flag something
    # to catch.
    correction = fit(
        [ramp],
        model=PolynomialModel(order=2),
        saturationKnee=None,
        badLinearityMedianMultiplier=5.0,
    )
    assert correction.badPixelMask[2, 2] & HIGH_FIT_RESIDUAL, (
        f"expected HIGH_FIT_RESIDUAL at (2,2); got mask = "
        f"{correction.badPixelMask[2, 2]:#06x}"
    )
    # Other pixels are clean.
    for r in range(4):
        for c in range(5):
            if (r, c) == (2, 2):
                continue
            assert (
                correction.badPixelMask[r, c] & HIGH_FIT_RESIDUAL
            ) == 0, f"unexpected HIGH_FIT_RESIDUAL at ({r}, {c})"
    # Summary captures the new fields.
    summary = correction.diagnostics.summary
    assert "badPixelFraction_highFitResidual" in summary
    assert summary["badPixelFraction_highFitResidual"] > 0
    assert "highFitResidualThresholdDN" in summary


def test_badLinearityMultiplierDisabledByPassingNone():
    from lsst.obs.pfs.h4Linearity.types import HIGH_FIT_RESIDUAL

    ramp = _buildSaturatingRamp()
    correction = fit(
        [ramp],
        model=PolynomialModel(order=2),
        saturationKnee=None,
        badLinearityMedianMultiplier=None,
    )
    assert ((correction.badPixelMask & HIGH_FIT_RESIDUAL) == 0).all()
