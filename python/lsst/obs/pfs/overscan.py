from functools import partial

import numpy as np

import lsst.afw.math as afwMath
from lsst.afw.image import MaskedImageF
from lsst.ip.isr.overscan import OverscanCorrectionTask
from lsst.pipe.base import Struct

from .rowBinning import getRowBinning
from .utils import isWindowed

__all__ = ("PfsOverscanCorrectionTask",)


class PfsOverscanCorrectionTask(OverscanCorrectionTask):
    """Task to correct overscan in PFS data.

    We override the base implementation in two cases:

    - windowed data (we didn't read all rows), in which case we can't do
      anything fancy with the overscan fitting.
    - data binned in the parallel direction (``W_CDRBIN > 1``), in which each
      binned row (including its serial overscan) is written ``W_CDRBIN`` times.
      We fit the serial overscan using one sample per binned row (i.e., one
      per read), and apply the same value to every row of the group, so that
      the rows of a group remain identical.
    """

    def correctOverscan(self, exposure, amp, imageBBox, overscanBBox, *args, **kwargs):
        """Trim the exposure, fit the overscan, subtract the fit, and
        calculate statistics.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure containing the data.
        amp : `lsst.afw.cameraGeom.Amplifier`
            The amplifier that is to be corrected.
        imageBBox: `lsst.geom.Box2I`
            Bounding box of the image data that will have the overscan
            subtracted.  If parallel overscan will be performed, that
            area is added to the image bounding box during serial
            overscan correction.
        overscanBBox: `lsst.geom.Box2I`
            Bounding box for the overscan data.
        isTransposed: `bool`
            If true, then the data will be transposed before fitting
            the overscan.
        leadingToSkip : `int`, optional
            Leading rows/columns to skip.
        trailingToSkip : `int`, optional
            Leading rows/columns to skip.
        overscanFraction : `float`, optional
            If the overscan region median is greater than overscanFraction
            and the imaging region median is greater than imageThreshold
            then overscan correction will be skipped.
        maxLevel : `float`, optional
            If the overscan region median is greater than overscanFraction
            and the imaging region median is greater than imageThreshold
            then overscan correction will be skipped.
        maskedRowColumnGrowSize : `int`, optional
            If a column (parallel overscan) or row (serial overscan) is
            completely masked, then grow the mask by this radius. If the
            value is <=0 then this will not be checked.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            ``ampOverscanModel``
                Overscan model broadcast to the full image size.
                (`lsst.afw.image.Exposure`)
            ``overscanOverscanModel``
                Overscan model broadcast to the full overscan image
                size. (`lsst.afw.image.Exposure`)
            ``overscanImage``
                Overscan image with the overscan fit subtracted.
                (`lsst.afw.image.Exposure`)
            ``overscanValue``
                Overscan model. (`float` or `np.array`)
            ``overscanMean``
                Mean value of the overscan fit. (`float`)
            ``overscanMedian``
                Median value of the overscan fit. (`float`)
            ``overscanSigma``
                Standard deviation of the overscan fit. (`float`)
            ``overscanMeanResidual``
                Mean value of the overscan region after overscan
                subtraction. (`float`)
            ``overscanMedianResidual``
                Median value of the overscan region after overscan
                subtraction. (`float`)
            ``overscanSigmaResidual``
                Standard deviation of the overscan region after
                overscan subtraction. (`float`)
        """
        args = (exposure, amp, imageBBox, overscanBBox) + args
        binning = getRowBinning(exposure.getMetadata())
        if binning > 1 and not isWindowed(exposure.getMetadata()):
            self.fitOverscanOriginal = self.fitOverscan
            try:
                self.fitOverscan = partial(
                    self.fitOverscanBinned,
                    binning=binning,
                    firstBinnedRow=amp.getRawDataBBox().getMinY(),
                    fitOverscan=self.fitOverscanOriginal,
                )
                return super().correctOverscan(*args, **kwargs)
            finally:
                self.fitOverscan = self.fitOverscanOriginal

        if not isWindowed(exposure.getMetadata()):
            return super().correctOverscan(*args, **kwargs)

        # Handling windowed data is the same, except for how we fit the overscan
        # The 'fitOverscan' method doesn't have access to the metadata, so we
        # need to trigger the change here, but we don't want to reproduce the
        # full method, so we temporarily swap out the method. It's a little
        # dirty, but it's cleaner in other respects.
        self.fitOverscanOriginal = self.fitOverscan
        try:
            self.fitOverscan = self.fitOverscanWindowed
            with np.errstate(invalid="ignore"):  # imageMedian=0, so overscanMedian/imageMedian=nan
                return super().correctOverscan(*args, **kwargs)
        finally:
            self.fitOverscan = self.fitOverscanOriginal

    def fitOverscanWindowed(self, overscanImage, isTransposed=False):
        value = afwMath.makeStatistics(overscanImage, afwMath.MEDIAN, self.statControl).getValue()
        return Struct(
            overscanValue=value,
            overscanMean=value,
            overscanMedian=value,
            overscanSigma=0.0,
        )

    def fitOverscanBinned(self, overscanImage, isTransposed=False, *, binning, firstBinnedRow, fitOverscan):
        """Fit the serial overscan of data binned in the parallel direction

        Rows before ``firstBinnedRow`` are unbinned; thereafter, each group of
        ``binning`` rows is a single read, repeated. We fit the overscan using
        a single row from each group (so each read gets equal weight, and
        the abscissa of the fit is the order of the reads), and then expand the
        fit back to all rows.

        Parameters
        ----------
        overscanImage : `lsst.afw.image.MaskedImage`
            Overscan image (in raw coordinates).
        isTransposed : `bool`
            Has the overscan been transposed (parallel overscan)?
        binning : `int`
            Number of rows binned together.
        firstBinnedRow : `int`
            Row (in raw coordinates) of the first binned group.
        fitOverscan : callable
            Original ``fitOverscan`` method.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result of ``fitOverscan``, expanded to all rows.
        """
        if isTransposed:
            # Parallel overscan: not affected by the parallel binning
            return fitOverscan(overscanImage, isTransposed)

        bbox = overscanImage.getBBox()
        rows = np.arange(bbox.getMinY(), bbox.getMaxY() + 1)
        group = np.where(rows < firstBinnedRow, rows,
                         firstBinnedRow + binning*((rows - firstBinnedRow)//binning))
        _, first, inverse = np.unique(group, return_index=True, return_inverse=True)

        compressed = MaskedImageF(overscanImage.getWidth(), len(first))
        compressed.mask.conformMaskPlanes(overscanImage.mask.getMaskPlaneDict())
        compressed.image.array[:] = overscanImage.image.array[first]
        compressed.mask.array[:] = overscanImage.mask.array[first]
        compressed.variance.array[:] = overscanImage.variance.array[first]

        result = fitOverscan(compressed, isTransposed)
        value = np.asarray(result.overscanValue)
        if value.ndim == 1:
            if value.size != len(first):
                raise RuntimeError(f"Unexpected overscan vector length {value.size} vs {len(first)}")
            result.overscanValue = value[inverse]
        return result
