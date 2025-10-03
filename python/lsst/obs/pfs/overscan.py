import numpy as np

import lsst.afw.math as afwMath
from lsst.ip.isr.overscan import OverscanCorrectionTask
from lsst.pipe.base import Struct

from .utils import isWindowed

__all__ = ("PfsOverscanCorrectionTask",)


class PfsOverscanCorrectionTask(OverscanCorrectionTask):
    """Task to correct overscan in PFS data.

    We override the base implementation only in the case we're dealing with
    windowed data (we didn't read all rows), in which case we can't do anything
    fancy with the overscan fitting.
    """

    def correctOverscan(self, exposure, *args, **kwargs):
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
        if not isWindowed(exposure.getMetadata()):
            return super().correctOverscan(exposure, *args, **kwargs)

        # Handling windowed data is the same, except for how we fit the overscan
        # The 'fitOverscan' method doesn't have access to the metadata, so we
        # need to trigger the change here, but we don't want to reproduce the
        # full method, so we temporarily swap out the method. It's a little
        # dirty, but it's cleaner in other respects.
        self.fitOverscanOriginal = self.fitOverscan
        try:
            self.fitOverscan = self.fitOverscanWindowed
            with np.errstate(invalid="ignore"):  # imageMedian=0, so overscanMedian/imageMedian=nan
                return super().correctOverscan(exposure, *args, **kwargs)
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
