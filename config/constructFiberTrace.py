config.interpolation = "POLYNOMIAL"
config.order = 5
config.xLow = -5.
config.xHigh = 5.
config.apertureFWHM = 2.5
config.signalThreshold = 100
config.nTermsGaussFit = 3
config.saturationLevel = 65000.
config.minLength = 3000
config.maxLength = 4096
config.nLost = 10
config.swathWidth = 300
config.telluric = "NONE"
config.overSample = 10
config.maxIterSig = 2
config.darkTime = None

config.isr.doBias = True
config.isr.doDark = True
config.isr.doFlat = False
config.isr.doLinearize = False

config.repair.cosmicray.nCrPixelMax = 5000000
config.repair.cosmicray.minSigma = 5.0
config.repair.cosmicray.min_DN = 50.0
