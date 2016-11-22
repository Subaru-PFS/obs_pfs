# Configuration for PFS ISR
from lsst.obs.subaru.isr import SubaruIsrTask
config.isr.retarget(SubaruIsrTask)
from lsst.obs.subaru.crosstalk import CrosstalkTask
config.isr.crosstalk.retarget(CrosstalkTask)

config.isr.overscanFitType = "AKIMA_SPLINE"
config.isr.overscanOrder = 30
config.isr.doBias = True
config.isr.doDark = True
config.isr.doFringe = False
config.isr.doWrite = False
config.isr.doCrosstalk = False
config.isr.doGuider = False
