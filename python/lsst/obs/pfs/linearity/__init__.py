"""nirLinearity — per-pixel nonlinearity correction for IR detector ramps."""

from __future__ import annotations

from .apply import apply, applyFrame
from .fit import fit
from .io import loadFits, saveFits
from .models import Model, PolynomialModel
from .types import (
    BORDER_PIX,
    FIT_FAILED,
    INSUFFICIENT_POINTS,
    MASKED_BY_INPUT,
    NON_MONOTONIC,
    ABOVE_VALID_RANGE,
    BELOW_VALID_RANGE,
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
    "MASKED_BY_INPUT",
    "INSUFFICIENT_POINTS",
    "FIT_FAILED",
    "NON_MONOTONIC",
    "BORDER_PIX",
    "ABOVE_VALID_RANGE",
    "BELOW_VALID_RANGE",
]
