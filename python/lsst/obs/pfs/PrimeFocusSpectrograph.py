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

"""Gen3 Butler registry declarations for Prime Focus Spectrograph.

This is a bare-bones placeholder, copied from obs_subaru and pared down.
It is not yet fully functional.
"""

__all__ = ("PrimeFocusSpectrograph",)

import os

from functools import lru_cache

from lsst.utils import getPackageDir
from lsst.afw.cameraGeom import makeCameraFromPath, CameraConfig
from lsst.obs.base import FitsRawFormatterBase
from lsst.daf.butler.core.utils import getFullTypeName
from lsst.obs.base import Instrument
from lsst.obs.base.gen2to3 import TranslatorFactory, PhysicalFilterToBandKeyHandler
from .pfsMapper import pfsFilterDefinitions
from .translator import PfsTranslator


class PrimeFocusSpectrographRawFormatter(FitsRawFormatterBase):
    """Gen3 Butler Formatters for PFS raw data.
    """
    translatorClass = PfsTranslator
    filterDefinitions = pfsFilterDefinitions

    def getDetector(self, id):
        return PrimeFocusSpectrograph().getCamera()[id]


class PrimeFocusSpectrograph(Instrument):
    """Gen3 Butler specialization class for Subaru's Prime Focus Spectrograph.
    """

    policyName = "pfs"
    obsDataPackage = "drp_pfs_data"
    filterDefinitions = pfsFilterDefinitions

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        packageDir = getPackageDir("obs_pfs")
        self.configPaths = [os.path.join(packageDir, "config"),
                            os.path.join(packageDir, "config", self.policyName)]

    @classmethod
    def getName(cls):
        # Docstring inherited from Instrument.getName
        return "PFS"

    def register(self, registry, update=False):
        # Docstring inherited from Instrument.register
        camera = self.getCamera()
        # The maximum values below make Gen3's ObservationDataIdPacker produce
        # outputs that match Gen2's ccdExposureId.
        obsMax = 21474800
        with registry.transaction():
            registry.syncDimensionData(
                "instrument",
                {
                    "name": self.getName(),
                    "detector_max": 10,
                    "visit_max": obsMax,
                    "exposure_max": obsMax,
                    "class_name": getFullTypeName(self),
                },
                update=update
            )
            for detector in camera:
                registry.syncDimensionData(
                    "detector",
                    {
                        "instrument": self.getName(),
                        "id": detector.getId(),
                        "full_name": detector.getName(),
                    },
                    update=update
                )
            self._registerFilters(registry, update=update)

    def getRawFormatter(self, dataId):
        # Docstring inherited from Instrument.getRawFormatter
        raise NotImplementedError("This method has not been implemented")

    def getCamera(self):
        """Retrieve the cameraGeom representation of HSC.

        This is a temporary API that should go away once obs_ packages have
        a standardized approach to writing versioned cameras to a Gen3 repo.
        """
        path = os.path.join(getPackageDir("obs_pfs"), self.policyName, "camera")
        return self._getCameraFromPath(path)

    @staticmethod
    @lru_cache()
    def _getCameraFromPath(path):
        """Return the camera geometry given solely the path to the location
        of that definition."""
        config = CameraConfig()
        config.load(os.path.join(path, "camera.py"))
        return makeCameraFromPath(
            cameraConfig=config,
            ampInfoPath=path,
            shortNameFunc=lambda name: name.replace(" ", "_"),
        )

    def makeDataIdTranslatorFactory(self) -> TranslatorFactory:
        # Docstring inherited from lsst.obs.base.Instrument.
        factory = TranslatorFactory()
        factory.addGenericInstrumentRules(self.getName())
        # Translate Gen2 `filter` to band if it hasn't been consumed
        # yet and gen2keys includes tract.
        factory.addRule(PhysicalFilterToBandKeyHandler(self.filterDefinitions),
                        instrument=self.getName(), gen2keys=("filter", "tract"), consume=("filter",))
        return factory
