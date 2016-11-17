import os.path

from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs", "isr.py"))

config.isr.doBias = True
config.isr.doDark = True
config.isr.doFlat = False
config.isr.doFringe = False
config.isr.doLinearize = False
