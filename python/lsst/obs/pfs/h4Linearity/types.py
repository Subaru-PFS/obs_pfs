"""Data types and bad-pixel flag constants for h4Linearity.

The four public dataclasses below form the I/O surface of the package:

- :class:`Ramp` — the input to :func:`apply`. A cumulative cube plus
  an optional per-pixel mask. The caller (typically PfsIsrTask) builds
  one of these from raw H4 reads.
- :class:`LinearizedRamp` — the output of :func:`apply`. Linearized
  cube + merged bad-pixel mask. The pipeline consumes the cumulative
  values and uses the mask to drive ``exposure.mask`` plane stamping.
- :class:`Diagnostics` — per-pixel fit-quality arrays plus a
  global summary dict. Round-trips through FITS via :func:`saveFits` /
  :func:`loadFits` and surfaces in the persisted calibration product.
- :class:`LinearityCorrection` — the fitted calibration object. The
  product butler stores it on disk; :func:`apply` consumes it.

The bad-pixel bit-flag constants are stored in the uint8
``badPixelMask`` field of both :class:`Ramp.validMask` (input) and
:class:`LinearizedRamp.badPixelMask` / :class:`LinearityCorrection.badPixelMask`
(output). Bits at or below ``0x10`` are set at fit time; bits above are
populated when :func:`apply` discovers per-read excursions outside the
fitted range.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

# Bad-pixel bit-flag constants. See spec section 3.
#
# Fit-time flags (set by :func:`fit` and ride through :func:`saveFits` /
# :func:`loadFits`):
MASKED_BY_INPUT: int = 0x01      #: input data masked by caller (defects, etc.)
INSUFFICIENT_POINTS: int = 0x02  #: too few usable reads to fit
FIT_FAILED: int = 0x04           #: normal-equations condition number > limit
NON_MONOTONIC: int = 0x08        #: fitted polynomial is non-monotonic on its range
BORDER_PIX: int = 0x10           #: detector-edge reference pixel (skipped by fit)

# Apply-time flags (set by :func:`apply` per-ramp from out-of-range reads):
BELOW_VALID_RANGE: int = 0x20    #: at least one read was below the per-pixel fitMin
ABOVE_VALID_RANGE: int = 0x40    #: at least one read was above the per-pixel fitMax


@dataclass(frozen=True)
class Ramp:
    """A single ramp's cumulative flux at each non-destructive read.

    The input form to :func:`apply`.

    Attributes
    ----------
    reads : np.ndarray
        Shape ``(N, H, W)``, float32. Cumulative signal, already
        photodiode-corrected and (typically) dark-subtracted. ``reads[n]``
        is the total accumulated ADU through read *n*; the input to
        :func:`fit` and :func:`apply` is the raw cumulative cube, not
        per-read deltas.
    validMask : np.ndarray or None
        Shape ``(H, W)``, integer-valued. 0 means the pixel is valid;
        any nonzero value flags it as bad input (e.g. caller-supplied
        defects). :func:`apply` merges this into the output's
        ``badPixelMask`` by OR'ing ``MASKED_BY_INPUT`` for any nonzero
        entry. ``None`` (default) means "all pixels valid".
    """

    reads: np.ndarray
    validMask: np.ndarray | None = None


@dataclass(frozen=True)
class LinearizedRamp:
    """Output of :func:`apply` on a :class:`Ramp`.

    Attributes
    ----------
    cumulativeLinear : np.ndarray
        Shape ``(N, H, W)``, float32. The linearized cumulative cube —
        same shape and read ordering as the input ``Ramp.reads``, but
        each valid pixel has been mapped through its per-pixel model.
        Bad pixels (any bit set in the merged mask) are passed through
        unchanged from input ``reads``.
    badPixelMask : np.ndarray
        Shape ``(H, W)``, uint8. Bitwise OR of: the fit-time mask from
        :class:`LinearityCorrection.badPixelMask`, the caller-supplied
        ``Ramp.validMask`` (OR'd as ``MASKED_BY_INPUT``), and the
        per-pixel out-of-range bits ``BELOW_VALID_RANGE`` /
        ``ABOVE_VALID_RANGE`` set when any read of the ramp falls
        outside that pixel's fitted ``[fitMin, fitMax]`` interval.
    """

    cumulativeLinear: np.ndarray   # (N, H, W) float32
    badPixelMask: np.ndarray       # (H, W) uint8


@dataclass(frozen=True)
class Diagnostics:
    """Per-pixel fit-quality diagnostics plus a dataset-wide summary.

    Carried inside :class:`LinearityCorrection` and round-tripped
    through FITS by :func:`saveFits` / :func:`loadFits`.

    Attributes
    ----------
    residualRms : np.ndarray
        ``(H, W)`` float32. Per-pixel RMS of the fit residual ``t - model(m)``
        over the reads used in the fit, in true-signal units.
    maxAbsResidual : np.ndarray
        ``(H, W)`` float32. Per-pixel maximum absolute residual.
    nPointsUsed : np.ndarray
        ``(H, W)`` int32. Number of reads actually used in the per-pixel
        fit (after exclusion of masked/saturated/below-threshold reads).
    monotonic : np.ndarray
        ``(H, W)`` bool. True when the fitted polynomial is monotonic on
        the per-pixel ``[fitMin, fitMax]`` interval. Non-monotonic pixels
        also get the ``NON_MONOTONIC`` bit in ``badPixelMask``.
    conditionNumber : np.ndarray
        ``(H, W)`` float32. Per-pixel condition number of the
        normal-equations matrix. Pixels above ``conditionNumberLimit``
        in :func:`fit` get ``FIT_FAILED``.
    summary : dict
        Dataset-wide scalar summary (counts of flagged pixels, fit
        parameters used, etc.). Persisted as FITS PRIMARY-header cards
        and recovered on :func:`loadFits`. Long Python keys round-trip
        via HIERARCH cards (case-preserved); short keys via the card
        comment.
    """

    residualRms: np.ndarray        # (H, W) float32
    maxAbsResidual: np.ndarray     # (H, W) float32
    nPointsUsed: np.ndarray        # (H, W) int32
    monotonic: np.ndarray          # (H, W) bool
    conditionNumber: np.ndarray    # (H, W) float32
    summary: dict[str, Any]


@dataclass(frozen=True)
class LinearityCorrection:
    """A fitted per-pixel nonlinearity correction.

    Produced by :func:`fit`, persisted by :func:`saveFits` / :func:`loadFits`,
    and consumed by :func:`apply`.

    Attributes
    ----------
    model : Model
        The model object (an instance of a class registered in
        ``MODEL_REGISTRY`` — currently :class:`PolynomialModel`). Carries
        the form-specific evaluation code; this dataclass just holds the
        fitted coefficients alongside it. Typed as ``Any`` here to avoid
        an import cycle with ``models`` — duck-typed against the
        :class:`Model` protocol at runtime.
    coefficients : np.ndarray
        Model-specific shape. For :class:`PolynomialModel(order=N)` this
        is ``(N+1, H, W)`` float32 — one Chebyshev coefficient per pixel
        per order.
    fitMin, fitMax : np.ndarray
        ``(H, W)`` float32. The per-pixel data range used in the fit (in
        cumulative-ADU units of the input ramp). :func:`apply` rescales
        ``m`` to Chebyshev's ``[-1, 1]`` domain over this interval and
        flags out-of-range reads via ``BELOW/ABOVE_VALID_RANGE``.
    badPixelMask : np.ndarray
        ``(H, W)`` uint8. Fit-time bad-pixel flags
        (``MASKED_BY_INPUT`` / ``INSUFFICIENT_POINTS`` / ``FIT_FAILED`` /
        ``NON_MONOTONIC`` / ``BORDER_PIX``). :func:`apply` ORs this with
        its own merged mask to produce the result's mask.
    diagnostics : Diagnostics
        Per-pixel fit-quality arrays + the global summary dict; see
        :class:`Diagnostics`.
    """

    model: Any                     # Model protocol; runtime-checked to avoid a types <-> models cycle
    coefficients: np.ndarray       # shape depends on model; polynomial: (order+1, H, W) float32
    fitMin: np.ndarray             # (H, W) float32
    fitMax: np.ndarray             # (H, W) float32
    badPixelMask: np.ndarray       # (H, W) uint8
    diagnostics: Diagnostics
