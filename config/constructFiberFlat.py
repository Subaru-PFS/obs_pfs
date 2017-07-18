import os.path
from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs.py"))
config.profile.load(os.path.join(getPackageDir("drp_stella"), "config", "createFlatFiberTraceProfile.py"))
config.trace.load(os.path.join(getPackageDir("drp_stella"), "config", "findAndTraceApertures.py"))

config.doRepair = True
config.psfFwhm = 2.5
config.psfSize = 21
config.crGrow = 2
config.minSNR = 50.

