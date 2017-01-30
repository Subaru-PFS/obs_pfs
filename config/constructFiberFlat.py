config.xOffsetHdrKeyWord = 'sim.slit.xoffset'
config.doRepair = True
config.psfFwhm = 3.0
config.psfSize = 21
config.crGrow = 2
config.darkTime = "DARKTIME"
config.display = True
config.minSNR = 10.

config.isr.doBias = True
config.isr.doDark = True
config.isr.doFlat = False
config.isr.doLinearize = False

config.repair.cosmicray.nCrPixelMax = 5000000
config.repair.cosmicray.minSigma = 5.0
config.repair.cosmicray.min_DN = 50.0

config.trace.interpolation = "POLYNOMIAL"
config.trace.order = 5
config.trace.xLow = -5.0
config.trace.xHigh = 5.0
config.trace.apertureFWHM = 2.5
config.trace.signalThreshold = 120.
config.trace.nTermsGaussFit = 3
config.trace.saturationLevel = 65000.0
config.trace.minLength = 3000
config.trace.maxLength = 4096
config.trace.nLost = 10

config.traceWide.interpolation = "POLYNOMIAL"
config.traceWide.order = 5
config.traceWide.xLow = -9.0
config.traceWide.xHigh = 9.0
config.traceWide.apertureFWHM = 4.5
config.traceWide.signalThreshold = 120.
config.traceWide.nTermsGaussFit = 3
config.traceWide.saturationLevel = 65000.0
config.traceWide.minLength = 3000
config.traceWide.maxLength = 4096
config.traceWide.nLost = 10

config.profile.swathWidth = 300
config.profile.telluric = "NONE"
config.profile.overSample = 10
config.profile.maxIterSig = 2
config.profile.lowerSigma = 3.0
config.profile.upperSigma = 3.0

