import os.path

from lsst.utils import getPackageDir

config.load(os.path.join(getPackageDir("obs_pfs"), "config", "pfs.py"))
config.profile.load(os.path.join(getPackageDir("drp_stella"), "config", "createFlatFiberTraceProfile.py"))
config.trace.load(os.path.join(getPackageDir("drp_stella"), "config", "findAndTraceApertures.py"))

config.trace.xLow = -5.4
config.trace.xHigh = 5.4
config.xOffsetHdrKeyWord = "sim.slit.xoffset"#"W_FCA_DITHER"

config.isr.doBias = True
config.isr.doDark = True
config.isr.doFlat = False
config.isr.doLinearize = False

config.repair.cosmicray.nCrPixelMax = 5000000
config.repair.cosmicray.minSigma = 5.0
config.repair.cosmicray.min_DN = 50.0
