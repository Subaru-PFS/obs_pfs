import os.path

from lsst.utils import getPackageDir

config.trace.load(os.path.join(getPackageDir("drp_stella"), "config", "findAndTraceApertures.py"))
config.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs.py"))
