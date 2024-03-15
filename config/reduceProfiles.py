import os.path
from lsst.utils import getPackageDir

reduceExposure = os.path.join(getPackageDir("obs_pfs"), "config", "reduceExposure.py")
config.reduceExposure.load(reduceExposure)
config.normalize.reduceExposure.load(reduceExposure)
