#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Gen3 Butler registry declarations for Prime Focus Spectrograph."""

__all__ = ("PrimeFocusSpectrograph", "PfsSimulator", "PfsDevelopment")

import os

from functools import lru_cache
from typing import List

from lsst.utils import getPackageDir
from lsst.daf.butler import Registry
from lsst.afw.cameraGeom import Camera
from lsst.obs.base import FitsRawFormatterBase
from lsst.utils.introspection import get_full_type_name
from lsst.obs.base import Instrument
from lsst.obs.base.gen2to3 import TranslatorFactory, PhysicalFilterToBandKeyHandler
from lsst.obs.base.yamlCamera import makeCamera
from lsst.daf.butler import DatasetType
from .pfsMapper import pfsFilterDefinitions
from .translator import PfsTranslator


@lru_cache
def loadCamera(category: str = None) -> Camera:
    """Retrieve the cameraGeom representation of an instrument

    Parameters
    ----------
    category : `str`, optional
        Category of PFS. This should be one of ``S`` (Subaru), ``L`` (LAM) or
        ``F`` (fake/simulator), or ``None``.

    Returns
    -------
    camera : `Camera`
        Description of instrument.
    """
    if category is not None:
        yamlFile = os.path.join(
            getPackageDir("obs_pfs"), "pfs-" + category, "camera", "camera.yaml"
        )
        if os.path.exists(yamlFile):
            return makeCamera(yamlFile)
    yamlFile = os.path.join(getPackageDir("obs_pfs"), "pfs", "camera", "camera.yaml")
    return makeCamera(yamlFile)


class PfsRawFormatter(FitsRawFormatterBase):
    """Gen3 Butler Formatters for PFS raw data."""

    translatorClass = PfsTranslator
    filterDefinitions = pfsFilterDefinitions
    pfsCategory = None

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


class PrimeFocusSpectrograph(Instrument):
    """Gen3 Butler specialization class for Subaru's Prime Focus Spectrograph."""

    policyName = "pfs"
    obsDataPackage = "drp_pfs_data"
    filterDefinitions = pfsFilterDefinitions
    pfsCategory = None
    standardCuratedDatasetTypes = frozenset(["defects"])

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        packageDir = getPackageDir("obs_pfs")
        self.configPaths = [os.path.join(packageDir, "config")]

    @classmethod
    def getName(cls) -> str:
        # Docstring inherited from Instrument.getName
        return "PFS" + (f"-{cls.pfsCategory}" if cls.pfsCategory is not None else "")

    def register(self, registry: Registry, update: bool = False) -> None:
        # Docstring inherited from Instrument.register
        camera = self.getCamera()
        instrument = self.getName()
        # The maximum values below make Gen3's ObservationDataIdPacker produce
        # outputs that match Gen2's ccdExposureId.
        obsMax = 21474800
        with registry.transaction():
            registry.syncDimensionData(
                "instrument",
                {
                    "name": instrument,
                    "detector_max": 20,
                    "visit_max": obsMax,
                    "exposure_max": obsMax,
                    "class_name": get_full_type_name(self),
                },
                update=update,
            )
            for spectrograph in (
                0,
                1,
                2,
                3,
                4,
            ):  # 0 means not really a spectrograph (e.g., guiders)
                registry.syncDimensionData(
                    "spectrograph",
                    dict(instrument=instrument, num=spectrograph),
                    update=update,
                )
            for arm in "brnmx":
                registry.syncDimensionData(
                    "arm", dict(instrument=instrument, name=arm), update=update
                )

            for detector in camera:
                name = detector.getName()
                isGuider = name.startswith("AG")
                if isGuider:
                    spectrograph = 0
                    arm = "x"
                else:
                    arm, spectrograph = name
                    spectrograph = int(spectrograph)

                registry.syncDimensionData(
                    "detector",
                    {
                        "instrument": instrument,
                        "spectrograph": spectrograph,
                        "arm": arm,
                        "id": detector.getId(),
                        "full_name": name,
                        "raft": "guider" if isGuider else name[-1],
                        "name_in_raft": name[-1] if isGuider else name[0],
                        "purpose": "guider" if isGuider else "science",
                    },
                    update=update,
                )
            self._registerFilters(registry, update=update)

        # Register types for datasets that will be ingested into the datastore.
        # Other types will be defined by the pipelines.
        self.registerDatasetType(
            registry,
            "pfsConfig",
            (
                "instrument",
                "exposure",
            ),
            "PfsConfig",
        )
        self.registerDatasetType(
            registry, "detectorMap_bootstrap", ("instrument", "detector"), "DetectorMap"
        )

        # self.registerDatasetType(registry, "fiberProfiles", ("instrument", "detector"), "FiberProfileSet",
        #                          True)
        # self.registerDatasetType(registry, "fiberProfiles_subset",
        #                          ("instrument", "detector", "pfs_design_id"), "FiberProfileSet")
        # self.registerDatasetType(registry, "fiberTraces", ("instrument", "exposure", "detector"),
        #                          "FiberTraceSet")
        # self.registerDatasetType(registry, "detectorMap", ("instrument", "detector"), "DetectorMap", True)
        # self.registerDatasetType(registry, "pfsArm", ("instrument", "exposure", "detector"), "PfsArm")
        # self.registerDatasetType(registry, "pfsMerged", ("instrument", "exposure"), "PfsMerged")
        # self.registerDatasetType(registry, "detectorMap_used", ("instrument", "exposure", "detector"),
        #                          "DetectorMap")
        # self.registerDatasetType(registry, "sky2d", ("instrument", "exposure", "detector"), "SkyModel")
        # self.registerDatasetType(registry, "sky1d", ("instrument", "exposure", "detector"),
        #                          "FocalPlaneFunction")
        # self.registerDatasetType(registry, "apCorr", ("instrument", "exposure", "detector"),
        #                          "FocalPlaneFunction")
        # self.registerDatasetType(registry, "pfsFluxReference", ("instrument", "exposure"), "PfsFluxReference")
        # self.registerDatasetType(registry, "fluxCal", ("instrument", "exposure"), "FocalPlaneFunction")
        # self.registerDatasetType(registry, "pfsCalibrated", ("instrument", "exposure"), "PfsTargetSpectra")
        # self.registerDatasetType(registry, "pfsCoadd", ("instrument", "skymap", "tract", "patch"),
        #                          "PfsTargetSpectra")
        # self.registerDatasetType(registry, "centroids", ("instrument", "exposure", "detector"), "ArcLineSet")
        # self.registerDatasetType(registry, "photometry", ("instrument", "exposure", "detector"), "ArcLineSet")
        # self.registerDatasetType(registry, "psf", ("instrument", "exposure", "detector"), "NevenPsf")
        # self.registerDatasetType(registry, "pfsArmLsf", ("instrument", "exposure", "detector"), "LsfDict")
        # self.registerDatasetType(registry, "pfsMergedLsf", ("instrument", "exposure"), "LsfDict")
        # self.registerDatasetType(registry, "pfsCalibratedLsf", ("instrument", "exposure"), "LsfDict")
        # self.registerDatasetType(registry, "pfsCoaddLsf", ("instrument", "skymap", "tract", "patch"),
        #                          "LsfDict")

    def registerDatasetType(
        self,
        registry: Registry,
        name: str,
        dimensions: List[str],
        storageClass: str,
        isCalibration: bool = False,
    ):
        """Register a datasetType

        Parameters
        ----------
        registry : `Registry`
            Butler datastore registry.
        name : `str`
            Name of the dataset type to register.
        dimensions : list of `str`
            Relevant dimensions for the dataset type.
        storageClass : `str`
            Name of the storage class for the dataset type.
        isCalibration : `bool`, optional
            Is this dataset type a calibration?

        Returns
        -------
        datasetType : `lsst.daf.butler.DatasetType`
            Dataset type.
        """
        datasetType = DatasetType(
            name,
            dimensions,
            storageClass,
            universe=registry.dimensions,
            isCalibration=isCalibration,
        )
        registry.registerDatasetType(datasetType)
        return datasetType

    def getRawFormatter(self, dataId):
        # Docstring inherited from Instrument.getRawFormatter
        return {
            None: PfsRawFormatter,
            "S": PfsRawFormatter,
            "F": PfsSimulatorRawFormatter,
            "L": PfsDevelopmentRawFormatter,
        }[self.pfsCategory]

    def getCamera(self):
        """Retrieve the cameraGeom representation

        This is a temporary API that should go away once obs_ packages have
        a standardized approach to writing versioned cameras to a Gen3 repo.
        """
        return loadCamera(self.pfsCategory)

    @classmethod
    def _getSpecificCuratedCalibrationPath(cls, datasetTypeName):
        """Return the path of the curated calibration directory.

        Parameters
        ----------
        datasetTypeName : `str`
            The name of the standard dataset type to find.

        Returns
        -------
        path : `str`
            The path to the standard curated data directory.  `None` if the
            dataset type is not found or the obs data package is not
            available.
        """
        return os.path.join(
            cls.getObsDataPackageDir(), "curated", cls.policyName, datasetTypeName
        )

    def makeDataIdTranslatorFactory(self) -> TranslatorFactory:
        # Docstring inherited from lsst.obs.base.Instrument.
        factory = TranslatorFactory()
        factory.addGenericInstrumentRules(self.getName())
        # Translate Gen2 `filter` to band if it hasn't been consumed
        # yet and gen2keys includes tract.
        factory.addRule(
            PhysicalFilterToBandKeyHandler(self.filterDefinitions),
            instrument=self.getName(),
            gen2keys=("filter", "tract"),
            consume=("filter",),
        )
        return factory


class PfsSimulator(PrimeFocusSpectrograph):
    """Instrument used for PFS simulator data"""

    pfsCategory = "F"


class PfsSimulatorRawFormatter(PfsRawFormatter):
    """Formatter used for PFS simulator data"""

    pfsCategory = "F"


class PfsDevelopment(PrimeFocusSpectrograph):
    """Instrument used for PFS development data from LAM"""

    pfsCategory = "L"


class PfsDevelopmentRawFormatter(PfsRawFormatter):
    """Formatter used for PFS development data from LAM"""

    pfsCategory = "L"
