import os.path
from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "calib.py"))
config.isr.load(os.path.join(getPackageDir("obs_pfs"), "config", "isr.py"))
config.repair.load(os.path.join(getPackageDir("obs_pfs"), "config", "repair.py"))

config.profiles.mask.append("BAD_FLAT")
