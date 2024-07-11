from typing import Literal, Optional

from lsst.obs.base.formatters.fitsExposure import FitsExposureFormatter
from lsst.obs.base import FitsRawFormatterBase

from .loadCamera import loadCamera
from .filterDefinitions import pfsFilterDefinitions
from .translator import PfsTranslator

__all__ = (
    "PfsRawFormatter",
    "PfsSimulatorRawFormatter",
    "PfsDevelopmentRawFormatter",
    "PfsFitsExposureFormatter",
)


class PfsRawFormatter(FitsRawFormatterBase):
    """Gen3 Butler Formatters for PFS raw data."""

    translatorClass = PfsTranslator
    filterDefinitions = pfsFilterDefinitions
    pfsCategory: Optional[Literal["F", "L", "S"]] = None

    def getDetector(self, id: int):
        """Return the detector that acquired this raw exposure.

        Parameters
        ----------
        id : `int`
            The identifying number of the detector to get.

        Returns
        -------
        detector : `~lsst.afw.cameraGeom.Detector`
            The detector associated with that ``id``.
        """
        return loadCamera(self.pfsCategory)[id]

    def makeWcs(self, visitInfo, detector):
        """Create a SkyWcs from information about the exposure.

        If VisitInfo is not None, use it and the detector to create a SkyWcs,
        otherwise return the metadata-based SkyWcs (always created, so that
        the relevant metadata keywords are stripped).

        Since PFS doesn't have a simple WCS for raw images, we simply disable
        this.

        Parameters
        ----------
        visitInfo : `~lsst.afw.image.VisitInfo`
            The information about the telescope boresight and camera
            orientation angle for this exposure.
        detector : `~lsst.afw.cameraGeom.Detector`
            The detector used to acquire this exposure.

        Returns
        -------
        skyWcs : `~lsst.afw.geom.SkyWcs`
            Reversible mapping from pixel coordinates to sky coordinates.

        Raises
        ------
        InitialSkyWcsError
            Raised if there is an error generating the SkyWcs, chained from the
            lower-level exception if available.
        """
        return None

    def _shiftAmpPixels(self, image):
        """Shift pixels in early raw frames.

        Early ADC versions insert a spurious pixel at the beginning of the
        readout. This affects all rows, pushing the last overscan pixel of a
        given row into the first leadin pixel on the following row.

        We strip out the 0th pixel, and leave the last one untouched.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Raw image; modified.
        """
        imArr = image.getArray()

        for amp in self.getDetector(self.observationInfo.detector_num):
            ySlice, xSlice = amp.getRawBBox().getSlices()
            ampIm = imArr[ySlice, xSlice]
            ampShape = ampIm.shape

            # Don't bother being tricky: make a copy.
            ampPixels = ampIm.flatten()
            ampPixels[:-1] = ampPixels[1:]

            imArr[ySlice, xSlice] = ampPixels.reshape(ampShape)

    def readImage(self):
        """Read just the image component of the Exposure.

        For PFS, we check an early ADC bug. Since this is almost certainly an
        FPGA bug, the decision is based on the FPGA version number. As of
        2016-12-01 the keyword is misnamed, so we can fix the format if the
        keyword does not exist.

        Returns
        -------
        image : `~lsst.afw.image.Image`
            In-memory image component.
        """
        image = super().readImage()

        try:
            dataVersion = int(self.metadata.get("W_VERSIONS_FPGA"), 16)
        except Exception:
            dataVersion = None
        if dataVersion is not None and dataVersion <= 0x0070:
            self._shiftAmpPixels(image)

        return image


class PfsSimulatorRawFormatter(PfsRawFormatter):
    """Formatter used for PFS simulator data"""

    pfsCategory = "F"


class PfsDevelopmentRawFormatter(PfsRawFormatter):
    """Formatter used for PFS development data from LAM"""

    pfsCategory = "L"


class PfsFitsExposureFormatter(FitsExposureFormatter):
    """FITS Exposure formatter that doesn't warn about filters

    The standard `FitsExposureFormatter` checks that the filter defined for the
    FITS file and the filter in the dataId match. But for PFS, we have different
    filters for different detectors rather than a single filter per exposure,
    so we suppress this check and the subsequent warning.
    """

    def _fixFilterLabels(self, file_filter_label, should_be_standardized=None):
        """DON'T compare the filter label read from the file with the one in the
        data ID.

        Parameters
        ----------
        file_filter_label : `lsst.afw.image.FilterLabel` or `None`
            Filter label read from the file, if there was one.
        should_be_standardized : `bool`, optional
            If `True`, expect ``file_filter_label`` to be consistent with the
            data ID and warn only if it is not.  If `False`, expect it to be
            inconsistent and warn only if the data ID is incomplete and hence
            the `FilterLabel` cannot be fixed.  If `None` (default) guess
            whether the file should be standardized by looking at the
            serialization version number in file, which requires this method to
            have been run after `readFull` or `readComponent`.

        Returns
        -------
        filter_label : `lsst.afw.image.FilterLabel` or `None`
            The preferred filter label; may be the given one or one built from
            the data ID.  `None` is returned if there should never be any
            filters associated with this dataset type.
        """
        return file_filter_label
