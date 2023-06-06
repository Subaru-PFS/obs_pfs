import os.path
from lsst.utils import getPackageDir

from lsst.obs.pfs.isrTask import PfsIsrTask
config.isr.retarget(PfsIsrTask)

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "calib.py"))
config.isr.load(os.path.join(getPackageDir("obs_pfs"), "config", "isr.py"))
config.repair.load(os.path.join(getPackageDir("obs_pfs"), "config", "repair.py"))

config.isr.doFlatNir = False

config.psfFwhm = 2.5
config.psfSize = 21
config.crGrow = 2
