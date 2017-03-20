import math
import numpy as np

from lsst.pex.config import Field, ListField
from lsst.pipe.drivers.constructCalibs import CalibConfig, CalibTask
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as afwcg

class PfsFlatCombineConfig(CalibConfig):
    """Configuration for flat construction.

    No changes required compared to the base class, but
    subclassed for distinction.
    """
    pass

class PfsFlatCombineTask(CalibTask):
    ConfigClass = PfsFlatCombineConfig
    _DefaultName = "flat"
    calibName = "flat"

    @classmethod
    def applyOverrides(cls, config):
        """Overrides to apply for bias construction"""
        config.isr.doFlat = False
        config.isr.doFringe = False
        config.isr.doLinearize = False
