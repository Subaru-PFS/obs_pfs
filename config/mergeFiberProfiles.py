import os.path
from lsst.utils import getPackageDir

config.normalize.reduceExposure.load(os.path.join(getPackageDir("obs_pfs"), "config", "reduceExposure.py"))
