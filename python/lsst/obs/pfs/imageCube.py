from typing import TYPE_CHECKING

import numpy as np
import astropy.io.fits

from astro_metadata_translator import fix_header
from lsst.afw.fits import readMetadata
from lsst.afw.image import ImageF

from pfs.datamodel.utils import astropyHeaderFromDict

from .translator import PfsTranslator

if TYPE_CHECKING:
    from lsst.daf.base import PropertyList

__all__ = ("ImageCube",)


class ImageCube:
    """A cube of images

    The images are stored in a FITS file, with each image in a separate HDU
    (in addition to a header-only primary HDU). The HDUs are named
    ``"IMAGE_<index>"``, where ``<index>`` is an integer. There is no assumption
    that the indices are contiguous, nor that they start at 0. There is no
    enforcement of common dimensions or type for the images.

    Images are not read until they are requested, and are cached once read.

    A new image cube can be created with the ``empty`` method, and images can be
    added with the ``__setitem__`` method, before calling ``write`` to save the
    cube to disk.

    This class holds a reference to the FITS file, so be sure to call it within
    a context manager (``with`` statement) to ensure the file is closed. You
    can use an instance of this class as that context manager, or explicitly
    delete the instance to close the file (if ``closeOnDel`` is True).

    Parameters
    ----------
    fits : `astropy.io.fits.HDUList`
        FITS file containing the images.
    metadata : `lsst.daf.base.PropertyList`
        Metadata (FITS header) for the images.
    closeOnDel : `bool`, optional
        If True, close the FITS file when the instance is deleted.
    """
    def __init__(
            self, fits: astropy.io.fits.HDUList, metadata: "PropertyList", closeOnDel: bool = True
    ) -> None:
        self.fits = fits
        self.metadata = metadata
        self._images: dict[int, ImageF] = {}
        self._closeOnDel = closeOnDel

    def __enter__(self) -> "ImageCube":
        """Enter context"""
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit context"""
        self.fits.close()

    def __del__(self):
        """Delete object"""
        if self._closeOnDel:
            self.fits.close()

    @classmethod
    def empty(cls, metadata: "PropertyList") -> "ImageCube":
        """Construct an empty cube

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`
            Metadata (FITS header) for the images.

        Returns
        -------
        cube : `ImageCube`
            An empty cube.
        """
        return cls(astropy.io.fits.HDUList(), metadata)

    @classmethod
    def fromFile(cls, path: str) -> "ImageCube":
        """Construct from a file

        To ensure the file is closed, it is recommended to use the instance
        returned from this method as a context manager, e.g.: ::
            with ImageCube.fromFile(path) as cube:
                # do something with cube

        Parameters
        ----------
        path : `str`
            Path to the FITS file.

        Returns
        -------
        cube : `ImageCube`
            The image cube.
        """
        fits = astropy.io.fits.open(path)
        metadata = readMetadata(path, 0)
        fix_header(metadata, translator_class=PfsTranslator, filename=path)
        return cls(fits, metadata)

    def flush(self) -> None:
        """Flush the cache"""
        self._images.clear()

    @classmethod
    def _getHduName(cls, index: int) -> str:
        """Return the name of the HDU for the image of interest

        Parameters
        ----------
        index : `int`
            The index of the image.

        Returns
        -------
        hduName : `str`
            The name of the HDU.
        """
        return f"IMAGE_{index}"

    @classmethod
    def _getHduIndex(cls, hduName: str) -> int:
        """Return the index of the image of interest

        Parameters
        ----------
        hduName : `str`
            The name of the HDU.

        Returns
        -------
        index : `int`
            The index of the image.
        """
        return int(hduName.split("_")[-1])

    def __getitem__(self, index: int) -> ImageF:
        """Return the image for the given index"""
        if index in self._images:
            return self._images[index]
        image = ImageF(self.fits[self._getHduName(index)].data.astype(np.float32))
        self._images[index] = image
        return image

    def __setitem__(self, index: int, image: ImageF) -> None:
        """Add the image for the given index"""
        self._images[index] = image

    def readAll(self) -> None:
        """Read all images into cache"""
        for hdu in self.fits[1:]:
            self[self._getHduIndex(hdu.name)] = ImageF(hdu.data.astype(np.float32))

    def writeFits(self, path: str) -> None:
        """Write the images

        Note that we only write images that have been explicitly read or set.

        Parameters
        ----------
        path : `str`
            Path to the output FITS file.
        """
        fits = astropy.io.fits.HDUList()

        fits.append(astropy.io.fits.PrimaryHDU(header=astropyHeaderFromDict(self.metadata)))
        for index in sorted(self._images):
            fits.append(astropy.io.fits.ImageHDU(self._images[index].array, name=self._getHduName(index)))
        with open(path, "wb") as fd:
            fits.writeto(fd)

    @classmethod
    def readFits(cls, path: str) -> "ImageCube":
        """Read the FITS file

        Parameters
        ----------
        path : `str`
            Path to the FITS file.

        Returns
        -------
        cube : `ImageCube`
            The image cube.
        """
        return cls.fromFile(path)
