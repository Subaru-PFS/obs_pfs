#import os.path

#from lsst.utils import getPackageDir

#config.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs", "isr.py"))

#from lsst.pipe.drivers.constructCalibs import CalibTask
##config.combination.retarget(CalibTask)
#from lsst.obs.pfs.detrends import PfsFlatCombineTask
#config.combination.retarget(PfsFlatCombineTask)
##config.combination.load(os.path.join(os.environ['OBS__DIR'], 'config', 'pfs', 'flat.py'))

config.isr.doBias = True
config.isr.doDark = True
config.isr.doFlat = False
config.isr.doFringe = False
#config.isr.doLinearize = False
