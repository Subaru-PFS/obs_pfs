from typing import Any, Literal, Optional
from importlib import import_module

import numpy as np

from lsst.obs.base.formatters.fitsExposure import FitsExposureFormatter
from lsst.obs.base.formatters.fitsGeneric import FitsGenericFormatter
from lsst.obs.base import FitsRawFormatterBase
from lsst.afw.image import Exposure, makeExposure, makeMaskedImage, FilterLabel
from lsst.resources import ResourcePath

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
        self._raw = PfsRaw(self.file_descriptor.location.path, self.pfsCategory)

    def read_from_local_file(self, path: str, component: str | None = None, expected_size: int = -1) -> Any:
        # Docstring inherited.
        if component == "ramp":
            if (readNum := self.checked_parameters.get("readNum")) is None:
                raise ValueError("No 'readNum' specified for ramp")
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

    def write_local_file(self, in_memory_dataset: Any, uri: ResourcePath) -> None:
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

    def write_local_file(self, in_memory_dataset: Any, uri: ResourcePath) -> None:
        """Write a Dataset.

        Workaround for DM-47152: FITS header keywords longer than 8 characters
        must have "HIERARCH " prepended to them, or there's a danger that a
        long string value will be damaged.

        Parameters
        ----------
        exposure : `object`
            The Dataset to store.
        """
        assert isinstance(in_memory_dataset, Exposure)
        exposure = in_memory_dataset

        metadata = exposure.getMetadata()
        fixes = []
        for key in metadata:
            if len(key) > 8 and not key.startswith("HIERARCH "):
                fixes.append(key)
        for key in fixes:
            value = metadata.get(key)
            comment = metadata.getComment(key)
            metadata.remove(key)
            metadata.add("HIERARCH " + key, value, comment)
        return super().write_local_file(exposure, uri)


class DetectorMapFormatter(FitsGenericFormatter):
    """Formatter for DetectorMap

    Provides components for the DetectorMap.
    """
    def read_from_local_file(self, path: str, component: str | None = None, expected_size: int = -1) -> Any:
        # Docstring inherited.
        from pfs.drp.stella import DetectorMap
        path = self.file_descriptor.location.path
        detectorMap = DetectorMap.readFits(path)
        if component is None:
            return detectorMap
        # This is inefficient, but it makes it very easy without implementing I/O
        if component == "bbox":
            return detectorMap.getBBox()
        if component == "slitOffsets":
            return np.vstack([detectorMap.getSpatialOffsets(), detectorMap.getSpectralOffsets()])
        raise ValueError(f"Unknown component {component}")


class PfsTargetSpectraFormatter(FitsGenericFormatter):
    """Formatter for PfsTargetSpectra

    Provides the 'single' component, that extracts a single spectrum.
    """
    PfsTargetSpectraClass: str = ""  # Subclasses must override

    def read_from_local_file(self, path: str, component: str | None = None, expected_size: int = -1) -> Any:
        # Docstring inherited.
        module = import_module("pfs.drp.stella.datamodel.pfsTargetSpectra", )
        cls = getattr(module, self.PfsTargetSpectraClass)

        path = self.file_descriptor.location.path
        spectra = cls.readFits(path)
        if component is None:
            return spectra

        parameters = self.file_descriptor.parameters
        if parameters is None:
            parameters = {}

        # This is inefficient, but it makes it very easy without implementing I/O
        if component == "single":
            # Don't care about catId, because all spectra in this file should have the same catId
            if "objId" in parameters and "obj_id" in parameters:
                raise ValueError("objId and obj_id both specified")
            objId = parameters.get("objId", parameters.get("obj_id"))
            select = [target for target in spectra if target.objId == objId]
            if not select:
                raise ValueError(f"objId={objId} not found in spectra at {path}")
            if len(select) > 1:
                raise ValueError(f"objId={objId} found multiple times in spectra at {path}")
            return spectra[select.pop()]
        raise ValueError(f"Unknown component {component}")


class PfsCalibratedSpectraFormatter(PfsTargetSpectraFormatter):
    PfsTargetSpectraClass = "PfsCalibratedSpectra"


class PfsObjectSpectraFormatter(PfsTargetSpectraFormatter):
    PfsTargetSpectraClass = "PfsObjectSpectra"
