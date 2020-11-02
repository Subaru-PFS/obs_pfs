import os.path
from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "calib.py"))
config.isr.load(os.path.join(getPackageDir("obs_pfs"), "config", "isr.py"))
config.repair.load(os.path.join(getPackageDir("obs_pfs"), "config", "repair.py"))

# Don't extend into neighbour's space
# This is specifically for flat creation, where the quartz spectra are narrowly spaced.
config.profiles.profileRadius = 3

# Disable cosmic-ray detection: it takes too long, and we weed them out through combining with rejection
config.doRepair = False
