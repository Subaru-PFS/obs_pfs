import os.path
from lsst.utils import getPackageDir

config.isr.load(os.path.join(getPackageDir("obs_pfs"), "config", "isr.py"))
config.isr.doFlat = False
config.findLines.threshold = 50.0
