"""Support for CCD data binned in the parallel (row) direction

When the PFS CCDs are binned by ``N`` in the parallel direction, the DAQ
still writes a full-sized image: each binned row (the sum of the charge in
``N`` physical rows, read once) is written ``N`` times. The binning is
arranged so that the first row of the data region (after the 48 rows that
precede it in the raw image) is the first row of a group of ``N`` identical
rows; the final group in the data region may be incomplete (its remaining
rows are in the parallel overscan, which contains no charge).

Because the raw image is full-sized, the ISR can run unchanged provided that
the calibrations are binned to match (see `RowBinning.replicate`). The
post-ISR image may then be collapsed to one row per group (see
`RowBinning.binMaskedImage`) for extraction.

The binning factor is recorded in the ``W_CDRBIN`` header keyword.
"""

from typing import Any, Mapping

import numpy as np

import lsst.geom
from lsst.afw.image import Exposure, ExposureF, MaskedImage, MaskedImageF

__all__ = ("getRowBinning", "RowBinning")


def getRowBinning(metadata: Mapping[str, Any] | None) -> int:
    """Return the parallel (row) binning factor from a header

    Parameters
    ----------
    metadata : `dict`-like or `None`
        Header metadata.

    Returns
    -------
    binning : `int`
        Number of physical rows summed into each binned row; 1 if not binned.
    """
    if metadata is None:
        return 1
    try:
        binning = int(metadata["W_CDRBIN"])
    except (KeyError, TypeError, ValueError):
        return 1
    return max(binning, 1)


class RowBinning:
    """Mapping between unbinned rows and rows binned by ``factor``

    Unbinned row ``y`` belongs to binned row ``y//factor``; the final binned
    row may contain fewer than ``factor`` unbinned rows.

    Coordinates: the centre of binned row ``j`` is at unbinned row
    ``factor*j + (factor - 1)/2`` (see `toUnbinned`, `toBinned`).

    Parameters
    ----------
    factor : `int`
        Binning factor (number of unbinned rows per binned row).
    height : `int`
        Number of unbinned rows.
    """

    def __init__(self, factor: int, height: int):
        if factor < 1:
            raise ValueError(f"Binning factor must be positive: {factor}")
        self.factor = int(factor)
        self.height = int(height)
        self.index = np.arange(self.height)//self.factor  # binned row for each unbinned row
        self.numBinned = int(self.index[-1]) + 1 if self.height > 0 else 0
        self.starts = np.arange(self.numBinned)*self.factor  # first unbinned row of each binned row
        self.counts = np.bincount(self.index, minlength=self.numBinned)  # unbinned rows per binned row
        # Mean unbinned row (of those that exist) in each binned row
        self.centers = np.bincount(self.index, weights=np.arange(self.height, dtype=float))/self.counts

    @classmethod
    def fromMetadata(cls, metadata: Mapping[str, Any] | None, height: int) -> "RowBinning":
        """Construct from header metadata

        Parameters
        ----------
        metadata : `dict`-like or `None`
            Header metadata (containing ``W_CDRBIN``).
        height : `int`
            Number of unbinned rows.
        """
        return cls(getRowBinning(metadata), height)

    def __bool__(self) -> bool:
        """Is the data actually binned?"""
        return self.factor > 1

    def __repr__(self) -> str:
        return f"{type(self).__name__}(factor={self.factor}, height={self.height})"

    def toUnbinned(self, yBinned):
        """Convert binned row coordinate(s) to unbinned row coordinate(s)"""
        return self.factor*np.asarray(yBinned) + 0.5*(self.factor - 1)

    def toBinned(self, yUnbinned):
        """Convert unbinned row coordinate(s) to binned row coordinate(s)"""
        return (np.asarray(yUnbinned) - 0.5*(self.factor - 1))/self.factor

    def _checkHeight(self, array: np.ndarray, axis: int):
        if array.shape[axis] != self.height:
            raise RuntimeError(f"Array has {array.shape[axis]} rows along axis {axis}; "
                               f"expected {self.height}")

    def sum(self, array: np.ndarray, axis: int = 0) -> np.ndarray:
        """Sum the unbinned rows into each binned row

        Parameters
        ----------
        array : `numpy.ndarray`
            Array to bin; ``array.shape[axis]`` must equal ``height``.
        axis : `int`
            Axis along which to bin.
        """
        self._checkHeight(array, axis)
        return np.add.reduceat(array, self.starts, axis=axis)

    def mean(self, array: np.ndarray, axis: int = 0) -> np.ndarray:
        """Average the unbinned rows in each binned row"""
        shape = [1]*array.ndim
        shape[axis] = self.numBinned
        return self.sum(np.asarray(array, dtype=float), axis)/self.counts.reshape(shape)

    def bitwiseOr(self, array: np.ndarray, axis: int = 0) -> np.ndarray:
        """OR together the unbinned rows in each binned row (for masks)"""
        self._checkHeight(array, axis)
        return np.bitwise_or.reduceat(array, self.starts, axis=axis)

    def bin(self, array: np.ndarray, how: str, axis: int = 0) -> np.ndarray:
        """Bin an array along ``axis``

        Parameters
        ----------
        array : `numpy.ndarray`
            Array to bin.
        how : `str`
            One of ``sum``, ``mean``, ``or``.
        axis : `int`
            Axis along which to bin.
        """
        if how == "sum":
            return self.sum(array, axis)
        if how == "mean":
            return self.mean(array, axis)
        if how == "or":
            return self.bitwiseOr(array, axis)
        raise ValueError(f"Unrecognised binning method: {how}")

    def expand(self, binned: np.ndarray, axis: int = 0) -> np.ndarray:
        """Replicate each binned row back to the unbinned rows it contains"""
        return np.take(binned, self.index, axis=axis)

    def replicate(self, array: np.ndarray, how: str, axis: int = 0) -> np.ndarray:
        """Bin and then expand, so that each row in a group has the binned value"""
        return self.expand(self.bin(array, how, axis), axis)

    def binMaskedImage(self, maskedImage: MaskedImage, imageHow: str = "mean") -> MaskedImageF:
        """Collapse a full-sized MaskedImage to one row per binned row

        This is appropriate for a post-ISR image in which each binned row has
        been written ``factor`` times: ``imageHow="mean"`` returns the value
        of the binned pixel (the sum of the charge in ``factor`` physical
        rows), with the variance of a single read.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Full-sized image; its height must equal ``height``.
        imageHow : `str`
            How to combine the image (and variance) rows: ``mean`` or ``sum``.

        Returns
        -------
        binned : `lsst.afw.image.MaskedImageF`
            Binned image, with ``xy0`` equal to the input's ``xy0``.
        """
        if maskedImage.getHeight() != self.height:
            raise RuntimeError(f"Image height {maskedImage.getHeight()} != {self.height}")
        bbox = lsst.geom.Box2I(maskedImage.getXY0(),
                               lsst.geom.Extent2I(maskedImage.getWidth(), self.numBinned))
        binned = MaskedImageF(bbox)
        binned.mask.conformMaskPlanes(maskedImage.mask.getMaskPlaneDict())
        binned.image.array[:] = self.bin(maskedImage.image.array, imageHow)
        binned.variance.array[:] = self.bin(maskedImage.variance.array, imageHow)
        binned.mask.array[:] = self.bitwiseOr(maskedImage.mask.array)
        return binned

    def binExposure(self, exposure: Exposure, imageHow: str = "mean") -> ExposureF:
        """Collapse a full-sized Exposure to one row per binned row

        The metadata, visitInfo, detector, etc. are copied from the input.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Full-sized exposure.
        imageHow : `str`
            How to combine the image (and variance) rows: ``mean`` or ``sum``.

        Returns
        -------
        binned : `lsst.afw.image.ExposureF`
            Binned exposure.
        """
        binned = ExposureF(self.binMaskedImage(exposure.maskedImage, imageHow))
        binned.setInfo(exposure.getInfo())
        binned.setMetadata(exposure.getMetadata().deepCopy())
        return binned

    def replicateMaskedImage(self, maskedImage: MaskedImage, imageHow: str):
        """Bin a full-sized MaskedImage in place, keeping its dimensions

        Each row in a group is replaced by the binned value, so the image is
        suitable for use as a calibration for full-sized binned raw data.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Image to modify in-place.
        imageHow : `str`
            How to combine the image rows: ``mean`` or ``sum``. The variance
            is combined consistently (``mean``: divided by the number of
            rows; ``sum``: summed).
        """
        if maskedImage.getHeight() != self.height:
            raise RuntimeError(f"Image height {maskedImage.getHeight()} != {self.height}")
        image = maskedImage.image.array
        variance = maskedImage.variance.array
        if imageHow == "mean":
            counts = self.counts[:, np.newaxis]
            newVariance = self.expand(self.sum(variance)/counts**2)
        elif imageHow == "sum":
            newVariance = self.replicate(variance, "sum")
        else:
            raise ValueError(f"Unrecognised binning method: {imageHow}")
        image[:] = self.replicate(image, imageHow)
        variance[:] = newVariance
        maskedImage.mask.array[:] = self.replicate(maskedImage.mask.array, "or")
