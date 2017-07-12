import os.path

from lsst.utils import getPackageDir

config.profile.load(os.path.join(getPackageDir("drp_stella"), "config", "createFlatFiberTraceProfile.py"))
config.trace.load(os.path.join(getPackageDir("drp_stella"), "config", "findAndTraceApertures.py"))

config.isr.doBias = True
config.isr.doDark = True
config.isr.doFlat = False
config.isr.doLinearize = False

config.repair.cosmicray.nCrPixelMax = 5000000
config.repair.cosmicray.minSigma = 5.0
config.repair.cosmicray.min_DN = 50.0

