import os.path

from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs", "isr.py"))
