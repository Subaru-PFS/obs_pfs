import os.path
from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "fiberTrace.py"))

# Don't extend into neighbour's space
# This is specifically for flat creation, where the quartz spectra are narrowly spaced.
config.trace.function.xLow = -3
config.trace.function.xHigh = 3
