config.xOffsetHdrKeyWord = 'sim.slit.xoffset'
config.doRepair = True
config.psfFwhm = 2.5
config.psfSize = 21
config.crGrow = 2
config.darkTime = "DARKTIME"
config.display = True
config.minSNR = 50.

config.isr.doBias = True
config.isr.doDark = True
config.isr.doFlat = False
config.isr.doLinearize = False

config.repair.cosmicray.nCrPixelMax = 5000000
config.repair.cosmicray.minSigma = 5.0
config.repair.cosmicray.min_DN = 50.0

config.trace.interpolation = "POLYNOMIAL"
config.trace.order = 5
config.trace.xLow = -5.5
config.trace.xHigh = 5.5
config.trace.apertureFWHM = 2.5
config.trace.signalThreshold = 120.
config.trace.nTermsGaussFit = 3
config.trace.saturationLevel = 65000.0
config.trace.minLength = 3000
config.trace.maxLength = 4096
config.trace.nLost = 10

config.profile.swathWidth = 300
config.profile.telluric = "NONE"
config.profile.overSample = 100
config.profile.maxIterSig = 2
config.profile.lowerSigma = 3.0
config.profile.upperSigma = 3.0

