import os.path
from lsst.utils

config.isr.doBias = True
config.isr.doDark = True
config.isr.doFlat = True
config.isr.doLinearize = False
config.doRepair = False
config.load(os.path.join(lsst.utils.getPackageDir("obs_pfs"), "config", "pfs.py"))
