import os.path

from lsst.utils import getPackageDir

from lsst.obs.pfs.isrTask import PfsIsrTask
config.isr.retarget(PfsIsrTask)

config.isr.load(os.path.join(getPackageDir("obs_pfs"), "config", "isr.py"))
config.isr.doWrite = True
