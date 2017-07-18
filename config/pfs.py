# Configuration for PFS ISR

if False:
    from lsst.obs.subaru.crosstalk import CrosstalkTask
    config.isr.crosstalk.retarget(CrosstalkTask)
    config.isr.doCrosstalk = False

if hasattr(config, 'ccdKeys'):
    config.ccdKeys = ['arm', 'spectrograph']
    
config.isr.expectWcs = False            # our spectrographs don't write a WCS to the header

config.isr.doLinearize = False
config.isr.doFringe = False
config.isr.doWrite = False

config.isr.overscanFitType = "AKIMA_SPLINE"
config.isr.overscanOrder = 30

if hasattr(config, 'repair'):
    config.repair.cosmicray.nCrPixelMax = 5000000
    config.repair.cosmicray.minSigma = 5.0
    config.repair.cosmicray.min_DN = 50.0

if hasattr(config, 'xOffsetHdrKeyWord'):
    config.xOffsetHdrKeyWord = 'sim.slit.xoffset'
