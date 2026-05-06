"""Apply a fitted LinearityCorrection to a new ramp or single cumulative frame."""

from __future__ import annotations

import numpy as np

from .types import (
    ABOVE_VALID_RANGE,
    BELOW_VALID_RANGE,
    MASKED_BY_INPUT,
    LinearityCorrection,
    LinearizedRamp,
    Ramp,
)


def apply(correction: LinearityCorrection, ramp: Ramp) -> LinearizedRamp:
    """Linearize a full ramp."""
    if ramp.reads.ndim != 3:
        raise ValueError(
            f"ramp.reads must be 3-D (N, H, W); got {ramp.reads.shape}"
        )
    if ramp.reads.shape[0] == 0:
        raise ValueError("ramp has zero reads")
    if ramp.reads.shape[1:] != correction.coefficients.shape[1:]:
        raise ValueError(
            f"ramp H,W = {ramp.reads.shape[1:]} does not match "
            f"correction H,W = {correction.coefficients.shape[1:]}"
        )

    m = ramp.reads.astype(np.float32)  # (N, H, W)

    # Map m → x ∈ [-1, 1] for Chebyshev evaluation
    denom = correction.fitMax - correction.fitMin
    denom = np.where(denom > 0, denom, 1.0)
    x = 2.0 * (m - correction.fitMin[None]) / denom[None] - 1.0

    t = correction.model.evaluate(correction.coefficients, x)
    below = m < correction.fitMin[None]
    above = m > correction.fitMax[None]

    # Merge the ramp's own validMask into the bad-pixel mask.
    bpm = correction.badPixelMask.copy()
    if ramp.validMask is not None:
        bpm[ramp.validMask != 0] |= MASKED_BY_INPUT

    # Bad-pixel pass-through: copy input m for any pixel with badPixelMask != 0.
    bad = bpm != 0
    if bad.any():
        t = np.where(bad[None], m, t)

    bpm[below.any(axis=0)] |= BELOW_VALID_RANGE
    bpm[above.any(axis=0)] |= ABOVE_VALID_RANGE

    return LinearizedRamp(
        cumulativeLinear=t.astype(np.float32),
        badPixelMask=bpm,
    )


def applyFrame(
    correction: LinearityCorrection, m: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Linearize a single already-cumulated frame."""
    m = np.asarray(m, dtype=np.float32)
    if m.shape != correction.coefficients.shape[1:]:
        raise ValueError(
            f"m shape {m.shape} does not match correction "
            f"H,W = {correction.coefficients.shape[1:]}"
        )

    # Map m → x ∈ [-1, 1] for Chebyshev evaluation
    denom = correction.fitMax - correction.fitMin
    denom = np.where(denom > 0, denom, 1.0)
    x = 2.0 * (m - correction.fitMin) / denom - 1.0

    t = correction.model.evaluate(correction.coefficients, x)
    oor = (m < correction.fitMin) | (m > correction.fitMax)

    bad = correction.badPixelMask != 0
    if bad.any():
        t = np.where(bad, m, t)

    return t.astype(np.float32), oor
