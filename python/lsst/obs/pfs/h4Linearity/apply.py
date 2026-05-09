"""Apply a fitted LinearityCorrection to a new ramp or single cumulative frame.

This module is the **production pipeline entry point** for nonlinearity
correction on H4 ramps. ``apply(correction, ramp)`` is the single call
PfsIsrTask makes per visit/detector after dark subtraction and before
CR/glitch detection; ``applyFrame`` is the per-frame variant used when
re-anchoring a partial range.
"""

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
    """Linearize a full ramp using a fitted per-pixel correction.

    Maps each cumulative-ADU value ``m`` through the per-pixel model
    ``t = model.evaluate(coefficients, x)`` where ``x`` is ``m`` rescaled
    to ``[-1, 1]`` over the per-pixel fit range ``[fitMin, fitMax]``
    (Chebyshev domain). The cube is processed read-by-read with vectorized
    numpy ops; no Python-level pixel loop.

    Pixel-by-pixel handling:

    - **Valid pixels** (``correction.badPixelMask == 0``, ``ramp.validMask``
      not flagging it): the linearized estimate ``t`` is returned.
    - **Bad pixels** (any bit set in ``correction.badPixelMask`` OR
      ``ramp.validMask``): the input value ``m`` is passed through
      unchanged. This is by design — the pipeline still wants a usable
      ADU value, just one not corrected for nonlinearity. The pixel
      stays flagged via the merged ``badPixelMask`` returned with the
      result.
    - **Out-of-range** (any read below ``fitMin`` or above ``fitMax``):
      the linearization formula is still evaluated (Chebyshev
      extrapolation), but the pixel's flag in the returned
      ``badPixelMask`` is OR'd with ``BELOW_VALID_RANGE`` /
      ``ABOVE_VALID_RANGE`` so callers can decide whether to use the
      extrapolated value. The current PfsIsrTask treats
      ``ABOVE_VALID_RANGE`` as ``SAT`` and ``BELOW_VALID_RANGE`` as
      informational (the rate is still well-defined there; CR detection
      still considers those pixels).

    The returned ramp is float32 cumulative-ADU values, same shape as
    ``ramp.reads``. ``correction.badPixelMask`` (the fit-time flags) is
    merged with ``ramp.validMask`` (a caller-supplied mask, e.g. defects)
    and with the freshly computed range flags, then returned as
    ``LinearizedRamp.badPixelMask``.

    Parameters
    ----------
    correction : LinearityCorrection
        Fitted correction from :func:`fit` or :func:`loadFits`. Must
        share ``(H, W)`` with ``ramp.reads``.
    ramp : Ramp
        Input ramp. ``reads`` must be 3-D ``(N, H, W)`` cumulative ADU.
        ``validMask`` is optional ``(H, W)`` — 0 means valid; any
        nonzero bit marks a pixel that should be passed through.

    Returns
    -------
    LinearizedRamp
        ``cumulativeLinear`` ``(N, H, W)`` float32 and ``badPixelMask``
        ``(H, W)`` uint8 (the merged fit-time + input + range bits).

    Raises
    ------
    ValueError
        If ``ramp.reads`` is not 3-D, has zero reads, or its ``(H, W)``
        does not match ``correction.coefficients[1:]``.
    """
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

    # copy=False: when the input is already float32 (the production
    # case after dark-subtraction), skip the wasted ~2.4 GB allocation
    # of an identical-dtype copy.
    m = ramp.reads.astype(np.float32, copy=False)  # (N, H, W)

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
    # In-place write avoids the ~2.4 GB extra cube allocation that
    # `t = np.where(bad[None], m, t)` would cost.
    bad = bpm != 0
    if bad.any():
        t[:, bad] = m[:, bad]

    bpm[below.any(axis=0)] |= BELOW_VALID_RANGE
    bpm[above.any(axis=0)] |= ABOVE_VALID_RANGE

    return LinearizedRamp(
        cumulativeLinear=t.astype(np.float32),
        badPixelMask=bpm,
    )


def applyFrame(
    correction: LinearityCorrection, m: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Linearize a single cumulative frame.

    Single-read variant of :func:`apply`. Used by PfsIsrTask when a
    partial ramp needs to be re-anchored at ``firstRead`` (the cumulative
    ADU at the first kept read is linearized once and subtracted from
    every read in the cube). Same model evaluation + bad-pixel
    pass-through as :func:`apply`, but does **not** return a merged
    bad-pixel mask — only the per-pixel out-of-range mask. Callers
    typically OR this into a larger mask themselves.

    Parameters
    ----------
    correction : LinearityCorrection
    m : np.ndarray
        Single cumulative-ADU frame, shape ``(H, W)``. Cast to float32
        internally.

    Returns
    -------
    t : np.ndarray
        ``(H, W)`` float32 linearized frame. Bad pixels (those with any
        bit in ``correction.badPixelMask``) are passed through unchanged.
    oor : np.ndarray
        ``(H, W)`` bool — True where ``m`` is below ``fitMin`` or above
        ``fitMax`` and therefore extrapolated.

    Raises
    ------
    ValueError
        If ``m.shape != correction.coefficients.shape[1:]``.
    """
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
        t[bad] = m[bad]

    return t.astype(np.float32, copy=False), oor
