import os.path
from lsst.utils import getPackageDir

config.reduceExposure.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs.py"))