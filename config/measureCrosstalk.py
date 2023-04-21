import os.path
from lsst.utils import getPackageDir

config.reduceExposure.load(os.path.join(getPackageDir("obs_pfs"), "config", "reduceExposure.py"))
config.reduceExposure.isr.doCrosstalk = False
config.reduceExposure.isr.growSaturationFootprintSize = 0  # Saturation spillover provides good signal
config.reduceExposure.isr.doDefect = False
config.reduceExposure.isr.doWidenSaturationTrails = False
config.reduceExposure.isr.doSaturationInterpolation = False
config.reduceExposure.isr.maskListToInterpolate = []
