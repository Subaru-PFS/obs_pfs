from typing import TYPE_CHECKING

import astropy.io.fits

from astro_metadata_translator import fix_header
from lsst.afw.fits import readMetadata
from lsst.afw.image import ExposureF


from .imageCube import ImageCube
from .translator import PfsTranslator

if TYPE_CHECKING:
    from lsst.daf.base import PropertyList

__all__ = ("PfsDark",)


class PfsDark:
    """Proxy for PFS darks

    PFS darks come in two flavors: CCD (for the b,r,m arms) and NIR (for the n
    arm). The CCD dark is an afw Exposure. The NIR dark is a multi-extension
    FITS file with a series of image HDUs, named "IMAGE_<index>", where <index>
    is an integer. This class provides a common interface for both.

    Parameters
    ----------
    path : `str`
        Path to the raw data file.

    See Also
    --------
    `lsst.obs.pfs.raw.PfsRaw`
    """
    def __init__(
        self,
        path: str | None = None,
        ccdDark: ExposureF | None = None,
        nirDark: ImageCube | None = None,
        metadata: "PropertyList | None" = None,
    ) -> None:
        if path is None and ccdDark is None and nirDark is None:
            raise ValueError("Must provide path, ccdDark, or nirDark")
        if ccdDark is not None and nirDark is not None:
            raise ValueError("Cannot provide both ccdDark and nirDark")
        if path is not None and (ccdDark is not None or nirDark is not None):
            raise ValueError("Cannot provide both path and ccdDark/nirDark")

        self.path = path
        self._metadata = metadata
        self._ccdDark = ccdDark
        self._nirDark = nirDark

    @classmethod
    def fromFile(cls, filename: str) -> "PfsDark":
        """Construct from a file

        Parameters
        ----------
        filename : `str`
            Path to the raw data file.

        Returns
        -------
        dark : `PfsDark`
            The dark.
        """
        return cls(filename)

    @classmethod
    def fromExposure(cls, exposure: ExposureF) -> "PfsDark":
        """Construct from a CCD dark

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The CCD dark.

        Returns
        -------
        dark : `PfsDark`
            The dark.
        """
        return cls(ccdDark=exposure, metadata=exposure.getMetadata())

    @classmethod
    def fromImageCube(cls, cube: ImageCube) -> "PfsDark":
        """Construct from a NIR dark

        Parameters
        ----------
        cube : `ImageCube`
            The NIR dark.

        Returns
        -------
        dark : `PfsDark`
            The dark.
        """
        return cls(nirDark=cube, metadata=cube.metadata)

    @property
    def metadata(self) -> "PropertyList":
        """Return the metadata"""
        if self._metadata is None:
            if self.path is not None:
                self._metadata = readMetadata(self.path, 0)
                fix_header(self._metadata, translator_class=PfsTranslator, filename=self.path)
            elif self._ccdDark is not None:
                self._metadata = self._ccdDark.getMetadata()
            elif self._nirDark is not None:
                self._metadata = self._nirDark.metadata
            else:
                raise RuntimeError("No metadata available")
        assert self._metadata is not None
        return self._metadata

    def isNir(self) -> bool:
        """Return if this is NIR data"""
        return self.metadata.get("W_ARM", 0) == 3

    def isCcd(self) -> bool:
        """Return if this is CCD data"""
        return not self.isNir()

    def getCcdDark(self, force: bool = False) -> ExposureF:
        """Return the CCD dark

        Parameters
        ----------
        force : `bool`, optional
            If True, return the CCD dark even if this is not CCD data.
            This is used to support using CCD-like darks for NIR data.

        Returns
        -------
        exposure : `ExposureF`
            The CCD dark.
        """
        if not force and not self.isCcd():
            raise RuntimeError("This is not a CCD")
        if self._ccdDark is not None:
            return self._ccdDark
        if self.path is None:
            raise RuntimeError("No path to CCD dark")
        return ExposureF(self.path)

    def getNirDark(self) -> ImageCube:
        """Return the NIR dark

        Returns
        -------
        cube : `ImageCube`
            The NIR dark.
        """
        if not self.isNir():
            raise RuntimeError("This is not NIR")
        if self._nirDark is not None:
            return self._nirDark
        if self.path is None:
            raise RuntimeError("No path to NIR dark")
        fits = astropy.io.fits.open(self.path)
        return ImageCube(fits, self.metadata)

    @classmethod
    def readFits(cls, path: str) -> "PfsDark":
        """Read the dark from a FITS file

        This is an alias for ``fromFile``.

        Parameters
        ----------
        path : `str`
            Path to the input FITS file.

        Returns
        -------
        dark : `PfsDark`
            The dark.
        """
        return cls.fromFile(path)

    def writeFits(self, path: str) -> None:
        """Write the dark to a FITS file

        Parameters
        ----------
        path : `str`
            Path to the output FITS file.
        """
        if self._ccdDark is not None and self._nirDark is not None:
            raise RuntimeError("Both CCD and NIR darks are set")
        if self._ccdDark is None and self._nirDark is None:
            raise RuntimeError("No darks are set")
        if self._ccdDark is not None:
            self._ccdDark.writeFits(path)
        if self._nirDark is not None:
            self._nirDark.write(path)
