import os.path
from lsst.utils import getPackageDir

from lsst.obs.pfs.isrTask import PfsIsrTask
config.isr.retarget(PfsIsrTask)

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "calib.py"))
config.isr.load(os.path.join(getPackageDir("obs_pfs"), "config", "isr.py"))

# Don't extend into neighbour's space
# This is specifically for flat creation, where the quartz spectra are narrowly spaced.
config.profiles.profileRadius = 3
