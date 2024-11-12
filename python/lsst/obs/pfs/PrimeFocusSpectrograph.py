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

from typing import List, Literal, Optional

from lsst.utils import getPackageDir
from lsst.daf.butler import Registry
from lsst.utils.introspection import get_full_type_name
from lsst.obs.base import Instrument
from lsst.daf.butler import DatasetType
from .filterDefinitions import pfsFilterDefinitions
from .formatters import PfsRawFormatter, PfsSimulatorRawFormatter, PfsDevelopmentRawFormatter
from .loadCamera import loadCamera

import lsst.obs.base._instrument
lsst.obs.base._instrument.StandardCuratedCalibrationDatasetTypes = dict(
    defects=dict(
        dimensions=("instrument", "detector", "arm", "spectrograph"),
        storageClass="Defects",
    ),
)


class PrimeFocusSpectrograph(Instrument):
    """Gen3 Butler specialization class for Subaru's Prime Focus Spectrograph."""

    policyName = "pfs"
    obsDataPackage = "drp_pfs_data"
    filterDefinitions = pfsFilterDefinitions
    pfsCategory: Optional[Literal["F", "L", "S"]] = None
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
            registry, "raw", ["instrument", "visit", "spectrograph", "arm"], "PfsRaw"
        )
        self.registerDatasetType(
            registry,
            "pfsConfig",
            ["instrument", "visit"],
            "PfsConfig",
        )
        self.registerDatasetType(
            registry,
            "detectorMap_calib",
            ["instrument", "arm", "spectrograph"],
            "DetectorMap",
            True,
        )

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

    def makeDefaultRawIngestRunName(self) -> str:
        """Make the default instrument-specific run collection string for raw
        data ingest.

        Returns
        -------
        coll : `str`
            Run collection name to be used as the default for ingestion of
            raws.
        """
        return self.makeCollectionName("raw", "sps")

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


class PfsSimulator(PrimeFocusSpectrograph):
    """Instrument used for PFS simulator data"""

    pfsCategory = "F"


class PfsDevelopment(PrimeFocusSpectrograph):
    """Instrument used for PFS development data from LAM"""

    pfsCategory = "L"
