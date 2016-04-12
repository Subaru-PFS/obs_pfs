import os.path

from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs", "isr.py"))

#from pfs.drp.stella.detrends import FlatCombineTask
#config.combination.retarget(FlatCombineTask)
#from pfs.drp.stella.detrends import PfsFlatCombineTask
#config.combination.retarget(PfsFlatCombineTask)
#config.combination.load(os.path.join(os.environ['OBS_SUBARU_DIR'], 'config', 'pfs', 'vignette.py'))
