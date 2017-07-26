import os.path
from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs.py"))
config.isr.doBias=True
config.isr.doDark=True
config.isr.doFlat=True
config.isr.doFringe=False
#config.isr.removePcCards=False
config.isr.doDefect=True
#config.isr.qa.doThumbnailOss=False
config.isr.doLinearize=False
#config.isr.qa.doThumbnailFlattened=False
#config.isr.doGuider=False
#config.isr.doWriteVignettePolygon=False
#config.isr.doCrosstalk=False
config.isr.doWrite=True
