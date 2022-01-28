import os.path
from lsst.utils import getPackageDir

from lsst.obs.pfs.isrTask import PfsIsrTask
config.reduceExposure.isr.retarget(PfsIsrTask)

config.reduceExposure.load(os.path.join(getPackageDir("obs_pfs"), "config", "reduceExposure.py"))
