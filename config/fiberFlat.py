import os.path
from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "fiberTrace.py"))

config.trace.function.xLow = -5.5
config.trace.function.xHigh = 5.5
