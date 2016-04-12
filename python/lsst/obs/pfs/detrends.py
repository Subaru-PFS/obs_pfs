import math
import numpy as np

from lsst.pex.config import Field, ListField
from lsst.pipe.drivers import CalibConfig, CalibTask
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
#    """Mask the vignetted area"""
    ConfigClass = PfsFlatCombineConfig
    _DefaultName = "flat"
    calibName = "flat"
#    FilterName = "NONE"

    @classmethod
    def applyOverrides(cls, config):
        """Overrides to apply for bias construction"""
        config.isr.doBias = False
        config.isr.doDark = False
        config.isr.doFlat = False
        config.isr.doFringe = False
#    
#    def __init__(self, *args, **kwargs):
#        super(PfsFlatCombineTask, self).__init__(*args, **kwargs)
#        
#    def run(self, sensorRefList, *args, **kwargs):
#        """Mask vignetted pixels after combining

#        This returns an Exposure instead of an Image, but upstream shouldn't
#        care, as it just dumps it out via the Butler.
#        """
#        combined = super(HscFlatCombineTask, self).run(sensorRefList, *args, **kwargs)
#        mi = afwImage.makeMaskedImage(combined.getImage())
#        mi.getMask().set(0)

#        # Retrieve the detector
#        # XXX It's unfortunate that we have to read an entire image to get the detector, but there's no
#        # public API in the butler to get the same.
#        image = sensorRefList[0].get("postISRCCD")
#        detector = image.getDetector()
#        del image

#        self.maskVignetting(mi.getMask(), detector)
#        self.maskBadAmps(mi.getMask(), detector)

#        return afwImage.makeExposure(mi)

