from typing import Any, Literal, Optional, Set

from lsst.daf.butler import StorageClassDelegate
from lsst.obs.base.formatters.fitsExposure import FitsExposureFormatter, standardizeAmplifierParameters
from lsst.obs.base import FitsRawFormatterBase
from lsst.afw.image import makeExposure, makeMaskedImage, FilterLabel

from .loadCamera import loadCamera
from .filterDefinitions import pfsFilterDefinitions
from .raw import PfsRaw
from .translator import PfsTranslator

__all__ = (
    "PfsRawFormatter",
    "PfsSimulatorRawFormatter",
    "PfsDevelopmentRawFormatter",
    "PfsFitsExposureFormatter",
)


class PfsRawFormatter(FitsRawFormatterBase):
    translatorClass = PfsTranslator
    filterDefinitions = pfsFilterDefinitions
    pfsCategory: Optional[Literal["F", "L", "S"]] = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._raw = PfsRaw(self.fileDescriptor.location.path, self.pfsCategory)

    def read(self, component: Optional[str] = None):
        # Docstring inherited.
        if component == "ramp":
            if (readNum := self.checked_parameters.get("readNum")) is None:
                raise ValueError(f"No 'readNum' specified for ramp")
            return self._raw.getCorrectedNirRead(readNum)
        if component == "exposure":
            exposure = makeExposure(makeMaskedImage(self._raw.getImage()))
            exposure.setDetector(self._raw.detector)
            info = exposure.getInfo()
            info.setVisitInfo(self._raw.visitInfo)
            info.setId(self._raw.visitInfo.id)
            info.setMetadata(self._raw.metadata)
            info.setDetector(self._raw.detector)
            arm = self._raw.obsInfo.ext_arm
            info.setFilter(FilterLabel(arm, arm))
            return exposure
        if component == "image":
            return self._raw.getImage()
        if component == "numReads":
            return self._raw.getNumReads()
        if component == "metadata":
            return self._raw.metadata
        if component == "obsInfo":
            return self._raw.obsInfo
        if component == "visitInfo":
            return self._raw.visitInfo
        if component == "detector":
            return self._raw.detector
        if component == "id":
            return self._raw.visitInfo.id
        if component == "bbox":
            return self._raw.bbox
        if component == "dimensions":
            return self._raw.dimensions
        if component == "xy0":
            return self._raw.xy0
        return self._raw

    def write(self, inMemoryDatset: Any):
        # Docstring inherited.
        raise NotImplementedError("PfsRawFormatter does not support writing")

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
