"""h4Linearity — per-pixel nonlinearity correction for H4 IR-detector ramps.

This package fits and applies per-pixel polynomial nonlinearity
corrections to up-the-ramp H4 reads. It is consumed by ``PfsIsrTask``
on every NIR visit/detector.

Pipeline path (the single API call the pipeline makes):

    correction = loadFits(<calib-product>)          # or fit(...) offline
    linearized = apply(correction, Ramp(reads=cube, validMask=...))

``correction`` is a :class:`LinearityCorrection` (model coefficients,
per-pixel valid range, fit-time bad-pixel mask, fit diagnostics).
``linearized`` is a :class:`LinearizedRamp` (cumulative-linear cube
plus a merged bad-pixel mask).

For partial-ramp re-anchoring the pipeline also uses :func:`applyFrame`
on a single cumulative read. CR / ASIC-glitch detection on the
linearized cube lives in the sibling :mod:`cr` module.

Public entry points:

- :func:`apply`, :func:`applyFrame` — linearize a ramp / a frame.
- :func:`fit` — fit a correction from training ramps.
- :func:`loadFits`, :func:`saveFits`, :func:`isH4LinearityFile` —
  read / write / sniff the on-disk FITS calibration product.
- :class:`Ramp`, :class:`LinearizedRamp`, :class:`LinearityCorrection`,
  :class:`Diagnostics` — I/O dataclasses.
- :class:`Model`, :class:`PolynomialModel` — model protocol + the
  built-in Chebyshev implementation.
- Bad-pixel bit-flag constants (``MASKED_BY_INPUT``, ``BELOW_VALID_RANGE``
  etc.) for decoding the masks returned by :func:`apply`.
"""

from __future__ import annotations

from . import cr
from . import rateStability
from .apply import apply, applyFrame
from .fit import fit
from .io import isH4LinearityFile, loadFits, saveFits
from .models import Model, PolynomialModel
from .types import (
    ABOVE_VALID_RANGE,
    ASIC_GLITCH,
    BELOW_VALID_RANGE,
    BORDER_PIX,
    DEAD,
    FIT_FAILED,
    GLITCH_MASKED,
    HIGH_FIT_RESIDUAL,
    INSUFFICIENT_POINTS,
    MASKED_BY_INPUT,
    NON_MONOTONIC,
    RATE_UNSTABLE,
    UNCLASSIFIED,
    UNSTABLE,
    Diagnostics,
    LinearityCorrection,
    LinearizedRamp,
    Ramp,
)

__all__ = [
    "Ramp",
    "LinearizedRamp",
    "Diagnostics",
    "LinearityCorrection",
    "Model",
    "PolynomialModel",
    "fit",
    "apply",
    "applyFrame",
    "saveFits",
    "loadFits",
    "isH4LinearityFile",
    "MASKED_BY_INPUT",
    "INSUFFICIENT_POINTS",
    "FIT_FAILED",
    "NON_MONOTONIC",
    "BORDER_PIX",
    "ABOVE_VALID_RANGE",
    "BELOW_VALID_RANGE",
    "UNCLASSIFIED",
    "UNSTABLE",
    "ASIC_GLITCH",
    "GLITCH_MASKED",
    "HIGH_FIT_RESIDUAL",
    "RATE_UNSTABLE",
    "DEAD",
]
