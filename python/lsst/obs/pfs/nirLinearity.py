from __future__ import annotations

from typing import TYPE_CHECKING  # noqa: F401

import numpy as np
import astropy.io.fits

from lsst.daf.base import PropertyList
from astro_metadata_translator import fix_header

from pfs.datamodel.utils import astropyHeaderFromDict

from .translator import PfsTranslator

__all__ = ("NirLinearity",)


class NirLinearity:
    """NIR Linearity correction data

    Consists of the following:
    - Primary HDU
    - Image HDU named ``LIMITS`` contains the per-pixel "full-well", as a
      detector-sized float array.
    - Image HDU named ``COEFFS`` contains the per-pixel polynomial coefficients
        for the linearity correction, as a cube of shape
        ``(NCOEFF, height, width)``.

    Parameters
    ----------
    limits : `np.ndarray`, shape ``(height, width)``
        The per-pixel "full-well".
    coeffs : `np.ndarray`, shape ``(NCOEFF, height, width)`
        The per-pixel polynomial coefficients for the linearity correction.
    metadata : `lsst.daf.base.PropertyList`
        Metadata (FITS header) for the images.
    method : `str`
        The method used to derive the linearity correction.
    """

    DEFAULT_METHOD = "np.polynomial.chebyshev"

    def __init__(
        self,
        limits: np.ndarray,
        coeffs: np.ndarray,
        metadata: "PropertyList",
        method: str | None = None,
    ) -> None:
        self.height, self.width = limits.shape
        if coeffs.ndim != 3 or coeffs.shape[1:] != limits.shape:
            raise ValueError(
                f"coeffs must have shape (NCOEFF, {self.height}, {self.width}), not {coeffs.shape}"
            )
        self.numCoeff = coeffs.shape[0]

        self.limits = limits
        self.coeffs = coeffs
        self.metadata = metadata
        self.method = method if method is not None else self.DEFAULT_METHOD

    @classmethod
    def empty(
        cls, shape: tuple, numCoeffs: int, metadata: "PropertyList", method: str | None = None
    ) -> "NirLinearity":
        """Construct an empty linearity correction

        Parameters
        ----------
        shape : `tuple`
            The shape of the detector, ``(height, width)``.
        numCoeffs : `int`
            The number of coefficients for the linearity correction.
        metadata : `lsst.daf.base.PropertyList`
            Metadata (FITS header) for the images.
        method : `str`, optional
            The method used to derive the linearity correction.

        Returns
        -------
        self : `NirLinearity`
            An empty linearity correction.
        """
        limits = np.zeros(shape, dtype=np.float32)
        coeffs = np.zeros((numCoeffs, shape[0], shape[1]), dtype=np.float32)
        return cls(limits, coeffs, metadata, method)

    @classmethod
    def fromFits(cls, fits: astropy.io.fits.HDUList) -> "NirLinearity":
        """Construct from a FITS file

        Parameters
        ----------
        fits : `astropy.io.fits.HDUList`
            The FITS file.

        Returns
        -------
        self : `NirLinearity`
            The linearity correction.
        """
        metadata = PropertyList.from_mapping(fits[0].header)
        fix_header(metadata, translator_class=PfsTranslator)
        limits = fits["LIMITS"].data.astype(np.float32, copy=False)
        coeffs = fits["COEFFS"].data.astype(np.float32, copy=False)
        method = fits[0].header.get("METHOD", cls.DEFAULT_METHOD)
        return cls(limits, coeffs, metadata, method)

    def toFits(self) -> astropy.io.fits.HDUList:
        """Write to a FITS file

        Returns
        -------
        fits : `astropy.io.fits.HDUList`
            The FITS file.
        """
        fits = astropy.io.fits.HDUList()

        header = astropyHeaderFromDict(self.metadata)
        header["METHOD"] = self.method
        fits.append(astropy.io.fits.PrimaryHDU(data=None, header=header))
        fits.append(astropy.io.fits.ImageHDU(self.limits.astype(np.float32), name="LIMITS"))
        fits.append(astropy.io.fits.ImageHDU(self.coeffs.astype(np.float32), name="COEFFS"))
        return fits

    def writeFits(self, path: str) -> None:
        """Write the linearity correction to a FITS file

        Parameters
        ----------
        path : `str`
            Path to the output FITS file.
        """
        with open(path, "wb") as fd:
            self.toFits().writeto(fd)

    @classmethod
    def readFits(cls, path: str) -> "NirLinearity":
        """Read the linearity correction from a FITS file

        Parameters
        ----------
        path : `str`
            Path to the FITS file.

        Returns
        -------
        self : `NirLinearity`
            The linearity correction.
        """
        with astropy.io.fits.open(path) as fits:
            return cls.fromFits(fits)
