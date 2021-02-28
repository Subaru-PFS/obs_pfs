import os.path
from lsst.utils import getPackageDir

from lsst.obs.pfs.isrTask import PfsIsrTask
config.isr.retarget(PfsIsrTask)

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "calib.py"))
config.isr.load(os.path.join(getPackageDir("obs_pfs"), "config", "isr.py"))
config.repair.load(os.path.join(getPackageDir("obs_pfs"), "config", "repair.py"))

config.profiles.mask.append("BAD_FLAT")
config.profiles.centerFit.order = 9
