import os.path

from lsst.utils import getPackageDir

config.isr.load(os.path.join(getPackageDir("obs_pfs"), "config", "isr.py"))

config.isr.overscanFitType = "AKIMA_SPLINE"
config.isr.overscanOrder = 30
config.isr.doBias = True
config.isr.doCrosstalk = False
config.isr.doDark = True
config.isr.doFlat = False
config.isr.doFringe = False
config.isr.doGuider = False
config.isr.doLinearize = False
config.isr.doDefect = False
config.isr.doWrite = True

