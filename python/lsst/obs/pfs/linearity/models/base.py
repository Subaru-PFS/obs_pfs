"""The Model protocol and BlockFitResult — shared across all concrete model forms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol, runtime_checkable

import numpy as np


@dataclass(frozen=True)
class BlockFitResult:
    """Everything a model's ``fitBlock`` must return for a single tile."""

    coefficients: np.ndarray       # model-specific shape; polynomial: (order+1, hTile, wTile) float32
    fitMin: np.ndarray             # (hTile, wTile) float32 — min m used per pixel
    fitMax: np.ndarray             # (hTile, wTile) float32 — max m used per pixel
    residualRms: np.ndarray        # (hTile, wTile) float32
    maxAbsResidual: np.ndarray     # (hTile, wTile) float32
    nPointsUsed: np.ndarray        # (hTile, wTile) int32
    conditionNumber: np.ndarray    # (hTile, wTile) float32
    monotonic: np.ndarray          # (hTile, wTile) bool
    badPixelMask: np.ndarray       # (hTile, wTile) uint8 — fit-time flags only
                                   # (INSUFFICIENT_POINTS, FIT_FAILED, NON_MONOTONIC)


@runtime_checkable
class Model(Protocol):
    """Protocol implemented by every concrete fit form (polynomial, spline, ...)."""

    modelName: str  # e.g. "CHEBYSHEV"; written into the FITS PRIMARY header

    def fitBlock(
        self,
        m: np.ndarray,                  # (nPoints, hTile, wTile) float32
        t: np.ndarray,                  # (nPoints,) float32
        valid: np.ndarray,              # (nPoints, hTile, wTile) bool
        conditionNumberLimit: float,
    ) -> BlockFitResult: ...

    def evaluate(
        self, coefficients: np.ndarray, m: np.ndarray
    ) -> np.ndarray: ...

    def isMonotonic(
        self, coefficients: np.ndarray, mMin: np.ndarray, mMax: np.ndarray
    ) -> np.ndarray: ...

    def toFitsHdus(self, correction: Any) -> list: ...  # list[astropy.io.fits.HDU]

    @classmethod
    def fromFitsHdus(cls, hdus: list) -> tuple["Model", np.ndarray]: ...
