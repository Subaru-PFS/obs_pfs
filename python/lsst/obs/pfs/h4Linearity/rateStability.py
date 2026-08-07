"""Per-pixel rate-stability rejection for linearized H4 ramps.

After linearization a pixel's accumulation rate should be constant up
the ramp. This module splits the linearized delta cube into two halves,
computes each half's plain mean of un-flagged deltas, and flags pixels
whose half rates disagree by more than a fractional threshold:

    diff      = r1 - r2
    noise_var = sem(r1)**2 + sem(r2)**2   # per-pixel; sem from un-flagged deltas
    adjusted  = sqrt(max(0, diff**2 - noise_var))
    fraction  = adjusted / max(|r1|, |r2|, rateFloorADU)
    reject    = fraction > threshold

The metric is noise-aware via per-pixel quadrature subtraction:
pure-noise pixels collapse to ~0 in expectation (their half-vs-half
``diff ** 2`` equals ``noise_var`` on average), so the threshold's
meaning is the fractional disagreement of the *signal*, not of the raw
difference. The estimator is delta-space throughout; the cumulative
cube is never reconstructed by this module.

Note: the per-half plain mean used here is NOT the pipeline science
rate. The pipeline computes a UTR-weighted rate over the full ramp
(:meth:`PfsIsrTask.calcUTRrateFromDeltas`); the plain half-mean is the
simpler half-vs-half consistency check this gate runs on top.
Diagnostic plotters that visualise pixel ramps should call the
production rate explicitly rather than reuse ``segmentRates``.

The pure-numpy algorithm lives here; the wrapper that integrates it
with ``PfsIsrTask.makeNirExposure`` lives in ``isrTask``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

__all__ = ["RateStabilityResult", "detectRateInstability"]


@dataclass
class RateStabilityResult:
    """Outcome of :func:`detectRateInstability`.

    Attributes
    ----------
    rejectMask : np.ndarray
        ``(H, W)`` bool; True where the pixel's half rates are
        inconsistent and the pixel should be masked.
    fraction : np.ndarray
        ``(H, W)`` float32; the noise-adjusted fractional disagreement
        ``sqrt(max(0, (r1-r2)**2 - sem(r1)**2 - sem(r2)**2)) /
        max(|r1|, |r2|, rateFloorADU)``. The numerator quadrature-
        subtracts the per-pixel standard error of the per-half mean
        rate so pure-noise pixels yield ~0. ``NaN`` where the pixel was
        untestable or excluded by ``goodPixelMask``.
    nTestable : np.ndarray
        ``(H, W)`` int16; number of segments with enough un-flagged
        deltas to contribute to the test (0, 1, or 2 with N=2).
    segmentRates : np.ndarray
        ``(2, H, W)`` float32; per-half plain mean of un-flagged
        deltas. ``NaN`` for empty segments. NOT the pipeline rate
        (UTR-weighted) — used internally by the rejection metric and
        retained for diagnostic introspection of the per-half means.
    goodPixelMask : np.ndarray
        ``(H, W)`` bool; the input mask that selected which pixels the
        test was allowed to consider. Echoed back so a downstream
        diagnostic can rerun the metric at a different threshold /
        rate floor / segment-size knob against the exact same pixel
        set the pipeline used.
    nRejected : int
        Count of rejected pixels.
    nUntestable : int
        Count of good pixels that could not be tested (< 2 testable
        segments).
    """

    rejectMask: np.ndarray
    fraction: np.ndarray
    nTestable: np.ndarray
    segmentRates: np.ndarray
    goodPixelMask: np.ndarray
    nRejected: int
    nUntestable: int


def detectRateInstability(
    deltas: np.ndarray,
    flagMask: np.ndarray,
    *,
    goodPixelMask: Optional[np.ndarray] = None,
    threshold: float = 0.20,
    rateFloorADU: float = 5.0,
    minDeltasPerSegment: int = 3,
) -> RateStabilityResult:
    """Flag pixels whose linearized ramp rate is not constant.

    Splits the ``(H, W, N-1)`` delta cube into two contiguous halves
    along the time axis; each half's ``segmentRates`` entry is the
    plain mean of its un-flagged deltas. A pixel is rejected when
    ``sqrt(max(0, (r1-r2)**2 - sem(r1)**2 - sem(r2)**2)) / max(|r1|,
    |r2|, rateFloorADU)`` exceeds ``threshold``, where ``sem(r_k)`` is
    the standard error of that per-half mean. The pipeline science
    rate is a UTR-weighted rate over the full ramp; this metric uses
    plain half-means because the test is a half-vs-half consistency
    check, not the science rate.

    Parameters
    ----------
    deltas : np.ndarray
        ``(H, W, N-1)`` linearized delta cube, post CR/glitch repair.
        Time axis last to match the H4 ISR convention.
    flagMask : np.ndarray
        ``(H, W, N-1)`` bool; True at deltas flagged as CR or glitch.
        Flagged deltas are excluded from each half's mean.
    goodPixelMask : np.ndarray, optional
        ``(H, W)`` bool; pixels where this is False are never tested or
        flagged. Defaults to all-True.
    threshold : float
        Rejection threshold for the fractional disagreement.
    rateFloorADU : float
        Lower bound on the denominator ``max(|r1|, |r2|, rateFloorADU)``;
        protects dark / near-zero pixels from blowing up the fraction
        from noise.
    minDeltasPerSegment : int
        Minimum un-flagged deltas for a segment to count as testable.

    Returns
    -------
    RateStabilityResult
    """
    deltas = np.asarray(deltas)
    if deltas.ndim != 3:
        raise ValueError(f"deltas must be 3-D (H, W, N-1); got {deltas.shape}")
    H, W, nDeltas = deltas.shape
    if nDeltas < 2 * minDeltasPerSegment:
        raise ValueError(
            f"deltas has only {nDeltas} frames; need at least "
            f"{2 * minDeltasPerSegment} to form two testable halves"
        )

    flagMask = np.asarray(flagMask, dtype=bool)
    if flagMask.shape != deltas.shape:
        raise ValueError(
            f"flagMask shape {flagMask.shape} != deltas shape {deltas.shape}"
        )

    if goodPixelMask is None:
        goodPixelMask = np.ones((H, W), dtype=bool)
    else:
        goodPixelMask = np.asarray(goodPixelMask, dtype=bool)
        if goodPixelMask.shape != (H, W):
            raise ValueError(
                f"goodPixelMask shape {goodPixelMask.shape} != (H, W) "
                f"({H}, {W})"
            )

    # Two contiguous halves of the delta axis, as equal as possible.
    bounds = np.linspace(0, nDeltas, 3).astype(int)

    # Per-segment sum of un-flagged deltas and un-flagged count. One
    # segment at a time keeps the accumulation loop's where-temp to
    # 1/2 of the cube rather than a whole-cube copy.
    segSum = np.zeros((2, H, W), dtype=np.float64)
    segSumSq = np.zeros((2, H, W), dtype=np.float64)
    segCount = np.zeros((2, H, W), dtype=np.int32)
    for k in range(2):
        lo, hi = int(bounds[k]), int(bounds[k + 1])
        valid = ~flagMask[..., lo:hi]
        contrib = np.where(valid, deltas[..., lo:hi], np.float32(0.0))
        segSum[k] = contrib.sum(axis=-1, dtype=np.float64)
        # Sum of squares over un-flagged deltas only (squaring of zero
        # contributions from flagged entries is fine -- they add nothing).
        segSumSq[k] = (contrib.astype(np.float64) ** 2).sum(axis=-1)
        segCount[k] = valid.sum(axis=-1)

    testable = segCount >= minDeltasPerSegment            # (2, H, W)
    nTestable = testable.sum(axis=0).astype(np.int16)     # (H, W)

    with np.errstate(invalid="ignore", divide="ignore"):
        segRate = (segSum / segCount).astype(np.float32)  # NaN where count 0

    tested = goodPixelMask & (nTestable >= 2)

    # Per-half sample variance: Var = (sum(X^2) - sum(X)^2/n) / (n-1).
    # Numerically OK because we operate at float64 in the accumulators.
    # nDof is the ddof=1 divisor; <=0 where there are not enough samples.
    nDof = (segCount - 1).astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        sampleVar = (segSumSq - segSum ** 2 / segCount) / nDof
    # Clip floating-point noise that can push a near-zero variance
    # slightly negative on noiseless inputs. NaN entries (segCount<2)
    # are harmless: those pixels are excluded from `tested` below.
    sampleVar = np.maximum(sampleVar, 0.0)
    # Standard error of the per-half mean rate. NaN-safe for untestable
    # halves; those pixels are filtered out below via `tested`.
    with np.errstate(invalid="ignore", divide="ignore"):
        sem = np.sqrt(sampleVar / segCount)              # (2, H, W) float64

    r1, r2 = segRate[0], segRate[1]
    with np.errstate(invalid="ignore"):
        diff = (r1.astype(np.float64) - r2.astype(np.float64))
        noiseVar = sem[0] ** 2 + sem[1] ** 2             # var(r1 - r2)
        adjustedSq = np.maximum(diff ** 2 - noiseVar, 0.0)
        adjusted = np.sqrt(adjustedSq).astype(np.float32)
        denom = np.maximum(
            np.maximum(np.abs(r1), np.abs(r2)),
            np.float32(rateFloorADU),
        )
        fraction = adjusted / denom

    fractionOut = np.where(tested, fraction, np.nan).astype(np.float32)
    rejectMask = tested & (fraction > threshold)

    nRejected = int(rejectMask.sum())
    nUntestable = int((goodPixelMask & (nTestable < 2)).sum())

    return RateStabilityResult(
        rejectMask=rejectMask,
        fraction=fractionOut,
        nTestable=nTestable,
        segmentRates=segRate,
        goodPixelMask=goodPixelMask,
        nRejected=nRejected,
        nUntestable=nUntestable,
    )
