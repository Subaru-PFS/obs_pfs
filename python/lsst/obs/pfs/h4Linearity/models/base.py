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
    """Protocol implemented by every concrete fit form (polynomial, spline, ...).

    Concrete models live in this subpackage (currently
    :class:`PolynomialModel`) and register themselves via
    :func:`registerModel`. The protocol decouples :func:`fit` and
    :func:`apply` from any particular functional form, and lets the
    FITS layer round-trip arbitrary model variants by name (via the
    ``MODEL`` PRIMARY-header card).
    """

    modelName: str  #: e.g. "CHEBYSHEV"; written into the FITS PRIMARY header.

    def fitBlock(
        self,
        m: np.ndarray,                  # (nPoints, hTile, wTile) float32
        t: np.ndarray,                  # (nPoints,) float32
        valid: np.ndarray,              # (nPoints, hTile, wTile) bool
        conditionNumberLimit: float,
    ) -> BlockFitResult:
        """Fit the model independently on every pixel of a tile.

        Parameters
        ----------
        m : np.ndarray
            ``(nPoints, hTile, wTile)`` float32. Per-pixel cumulative-ADU
            values to use as the model's *measured* coordinate.
        t : np.ndarray
            ``(nPoints,)`` float32. Reference (true) signal corresponding
            to each input point — shared across all pixels of the tile.
        valid : np.ndarray
            ``(nPoints, hTile, wTile)`` bool. Per-pixel-per-read mask;
            False entries are excluded from the fit.
        conditionNumberLimit : float
            Pixels whose normal-equations matrix exceeds this condition
            number are marked ``FIT_FAILED`` and left with zero
            coefficients.
        """

    def evaluate(
        self, coefficients: np.ndarray, m: np.ndarray
    ) -> np.ndarray:
        """Map measured cumulative ADU ``m`` to linearized signal ``t``.

        Per-pixel vectorized; ``coefficients`` and ``m`` broadcast over
        the same trailing ``(H, W)`` axes.
        """

    def isMonotonic(
        self, coefficients: np.ndarray, mMin: np.ndarray, mMax: np.ndarray
    ) -> np.ndarray:
        """Return a per-pixel bool mask: True where the model is monotonic on ``[mMin, mMax]``.

        Used during :func:`fit` to flag ``NON_MONOTONIC`` pixels.
        """

    def toFitsHdus(self, correction: Any) -> list:  # list[astropy.io.fits.HDU]
        """Serialize the model's coefficients + form metadata to FITS HDUs.

        The list returned here is inserted after the PRIMARY HDU and
        before the standard image HDUs by :func:`saveFits`. ``correction``
        is the surrounding :class:`LinearityCorrection`, so the model
        can pull whatever per-pixel arrays it owns.
        """

    @classmethod
    def fromFitsHdus(cls, hdus: list) -> tuple["Model", np.ndarray]:
        """Rebuild ``(model, coefficients)`` from a list of HDUs.

        Inverse of :meth:`toFitsHdus`. The full HDU list (including the
        PRIMARY and the standard non-model HDUs) is passed in; the
        model picks out the entries it owns.
        """
