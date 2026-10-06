from typing import TYPE_CHECKING, Optional

import numpy as np
import fitsio

from astro_metadata_translator import fix_header
from lsst.afw.fits import readMetadata
from lsst.afw.image import ImageF

from .translator import PfsTranslator

if TYPE_CHECKING:
    from lsst.daf.base import PropertyList

__all__ = ("ImageCube",)


class ImageCube:
    """A cube of images

    The images are stored in a FITS file, with each image in a separate HDU
    (in addition to a header-only primary HDU). The HDUs are named
    ``"IMAGE_<index>"``, where ``<index>`` is a 1-based integer. The image ids outside
    the file are 0-based.
    There is no enforcement of common dimensions or type for the images.

    Images are not read until they are requested, and are cached once read.

    A new image cube can be created with the ``empty`` method, and images can be
    added with the ``__setitem__`` method, before calling ``write`` to save the
    cube to disk.

    The file is read and written with fitsio (cfitsio). Images are written
    RICE-compressed with a fixed quantization step of 0.1 (in pixel units,
    i.e. e- for the NIR cubes), subtractively dithered. cfitsio's default
    picks the step per tile from a noise estimate, which images with vertical
    spectra break; a fixed, dithered step bounds the error at half a step
    everywhere and leaves it unbiased.

    A cube read from a file holds the file open, so be sure to use it within a
    context manager (``with`` statement) to ensure the file is closed. You can
    use an instance of this class as that context manager, or explicitly delete
    the instance to close the file.

    Parameters
    ----------
    metadata : `lsst.daf.base.PropertyList`
        Metadata (FITS header) for the images.
    path : `str`, optional
        Path to the FITS file holding the images; `None` for a cube that exists
        only in memory.
    """

    compression = dict(compress="RICE", qlevel=-0.1, qmethod="SUBTRACTIVE_DITHER_2", dither_seed=-1)
    """Compression for the image HDUs: a negative ``qlevel`` is an absolute
    quantization step; ``SUBTRACTIVE_DITHER_2`` dithers and preserves zeros;
    ``dither_seed=-1`` seeds the dither from the data checksum, so a rewrite
    of the same data gives the same values."""

    def __init__(self, metadata: "PropertyList", path: Optional[str] = None) -> None:
        self.metadata = metadata
        self._images: dict[int, ImageF] = {}
        self._path = path
        self._reader: Optional[fitsio.FITS] = None

        self.nreads = 0
        if path is not None:
            reader = self.reader
            header = reader[0].read_header()
            if "NREADS" in header:
                self.nreads = header["NREADS"]
            else:
                indices = [self._getHduIndex(name) for name in self._hduNames()]
                self.nreads = max(indices) + 1 if indices else 0

    @property
    def reader(self) -> Optional[fitsio.FITS]:
        """The cfitsio reader for the backing file, or None if in memory."""
        if self._reader is None and self._path is not None:
            self._reader = fitsio.FITS(self._path)
        return self._reader

    def _hduNames(self) -> list[str]:
        """Return the names of the image HDUs in the file."""
        if self.reader is None:
            return []
        return [hdu.get_extname() for hdu in self.reader[1:] if hdu.get_extname().startswith("IMAGE_")]

    def _readHduArray(self, index: int) -> np.ndarray:
        """Read one image from the file, without caching it anywhere."""
        name = self._getHduName(index)
        reader = self.reader
        if reader is None:
            # In-memory cube (``empty``/``fromCube``): no file to read from.
            raise KeyError(name)
        try:
            hdu = reader[name]
        except Exception as exc:
            # cfitsio reports a missing extension its own way; callers rely on
            # a missing read raising KeyError, which is how a dark that is
            # shorter than the ramp gets caught rather than passing silently.
            raise KeyError(name) from exc
        return np.asarray(hdu.read(), dtype=np.float32)

    def _closeReader(self) -> None:
        if self._reader is not None:
            self._reader.close()
            self._reader = None

    def __enter__(self) -> "ImageCube":
        """Enter context"""
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit context"""
        self._closeReader()

    def __del__(self):
        """Delete object"""
        self._closeReader()

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
        return cls(metadata)

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
        metadata = readMetadata(path, 0)
        fix_header(metadata, translator_class=PfsTranslator, filename=path)
        return cls(metadata, path=path)

    @classmethod
    def fromCube(cls, data: np.ndarray, metadata: "PropertyList") -> "ImageCube":
        """Construct ourselves from a data cube

        Note: this does *not* copy the data

        Parameters
        ----------
        data : `np.ndarray`
            The data to load from.
        metadata : `lsst.daf.base.PropertyList`
            Metadata (FITS header) for the images.

        Returns
        -------
        cube : `ImageCube`
            The image cube.
        """

        self = cls.empty(metadata)
        for i in range(len(data)):
            self[i] = ImageF(data[i].astype(np.float32, copy=False))
        self.nreads = len(data)
        return self

    def flush(self) -> None:
        """Flush the cache"""
        self._images.clear()

    @classmethod
    def _getHduName(cls, index: int) -> str:
        """Return the name of the HDU for the image of interest

        The HDU name is ``"IMAGE_<N>"``, where ``<N>`` is a 1-based integer.

        Parameters
        ----------
        index : `int`
            The index of the image. 0-based.

        Returns
        -------
        hduName : `str`
            The name of the HDU.
        """
        return f"IMAGE_{index+1}"

    @classmethod
    def _getHduIndex(cls, hduName: str) -> int:
        """Return the index of the image of interest

        The HDU name is ``"IMAGE_<N>"``, where ``<N>`` is a 1-based integer.
        The read index is 0-based.

        Parameters
        ----------
        hduName : `str`
            The name of the HDU.

        Returns
        -------
        index : `int`
            The index of the image. 0-based
        """
        return int(hduName.split("_")[-1])-1

    def __getitem__(self, index: int) -> ImageF:
        """Return the image for the given index"""
        if index in self._images:
            return self._images[index]
        image = ImageF(self._readHduArray(index))
        self._images[index] = image
        return image

    def __setitem__(self, index: int, image: ImageF) -> None:
        """Add the image for the given index"""
        if not isinstance(index, int):
            raise ValueError(f"Index must be an integer, not {index!r}")
        self._images[index] = image

    def getNumReads(self) -> int:
        """Return the number of images in the cube.

        Sugar to match the method name in PfsRaw
        """
        return self.nreads

    def readAll(self) -> None:
        """Read all images into cache"""
        for name in self._hduNames():
            index = self._getHduIndex(name)
            self[index] = ImageF(self._readHduArray(index))

    def getReadArray(self, index: int) -> np.ndarray:
        """Return the image for the given index, but do *not* cache it if it is not already cached

        Parameters
        ----------
        index : `int`
            The index of the image. 0-based.

        Returns
        -------
        image : `np.ndarray`
            The image.
        """
        if index in self._images:
            return self._images[index].array
        return self._readHduArray(index)

    def getImageCube(self, nreads: Optional[int] = None) -> np.ndarray:
        """Return the image cube as a 3D numpy array, trying not to cache new reads.

        The images are stacked along the first axis, so the shape of the
        returned array is (nreads, height, width).

        Parameters
        ----------
        nreads : `int`
            The number of reads to return. If None, return all reads.

        Returns
        -------
        imageCube : `np.ndarray`
            The (nreads, height, width) ndarray.
        """

        if nreads is None:
            nreads = self.nreads
        if nreads <= 0:
            return np.empty((0, 0, 0), dtype=np.float32)
        first = self.getReadArray(0)
        ret = np.empty((nreads, *first.shape), dtype=np.float32)
        ret[0] = first
        for i in range(1, nreads):
            ret[i] = self.getReadArray(i)
        return ret

    def writeFits(self, path: str) -> None:
        """Write the images

        Note that we only write images that have been explicitly read or set.

        Parameters
        ----------
        path : `str`
            Path to the output FITS file.
        """
        with fitsio.FITS(path, "rw", clobber=True) as fits:
            fits.write(None, header=self._makeHeader(self.metadata, self.nreads))
            for index in sorted(self._images):
                fits.write(self._images[index].array, extname=self._getHduName(index), **self.compression)

    @staticmethod
    def _makeHeader(metadata: "PropertyList", nreads: int) -> list[dict]:
        """Make the primary-HDU header records

        ``COMMENT`` and ``HISTORY`` cards are not preserved. Keywords longer
        than 8 characters are written with the HIERARCH convention.

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList` or `dict`
            FITS header keywords and values.
        nreads : `int`
            Number of reads in the cube, written as ``NREADS``.

        Returns
        -------
        records : `list` [`dict`]
            Header records for `fitsio`.
        """
        records = [dict(name=key, value=value) for key, value in metadata.items()
                   if key not in ("HISTORY", "COMMENT", "NREADS")]
        records.append(dict(name="NREADS", value=nreads))
        return records

    @classmethod
    def writeCubeData(cls, path: str, data: np.ndarray, metadata: "PropertyList") -> None:
        """Directly write the data to a FITS file

        Each image is written as soon as it is converted, so the whole cube is
        never held twice.

        Parameters
        ----------
        path : `str`
            Path to the output FITS file.
        data : `np.ndarray`
            The data to write.
        metadata : `lsst.daf.base.PropertyList`
            Metadata (FITS header) for the images.
        """
        with fitsio.FITS(path, "rw", clobber=True) as fits:
            fits.write(None, header=cls._makeHeader(metadata, len(data)))
            for i in range(len(data)):
                image = np.asarray(data[i], dtype=np.float32)
                fits.write(image, extname=cls._getHduName(i), **cls.compression)

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
