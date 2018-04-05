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
config.isr.expectWcs=False
config.isr.doAddDistortionModel=False

config.isr.overscanFitType = "AKIMA_SPLINE"
config.isr.overscanOrder = 30

if hasattr(config, 'repair'):
    # This CR tuning comes from looking at simulated r1 images w/o scattering.
    #

    config.repair.interp.modelPsf.defaultFwhm = 1.75
    config.repair.cosmicray.cond3_fac = 4
    config.repair.cosmicray.cond3_fac2 = 1
    config.repair.cosmicray.nCrPixelMax = 60000
    config.repair.cosmicray.minSigma = 10.0
    config.repair.cosmicray.min_DN = 500.0

    # Interpolation needs more work. In the meanwhile do not hide what is there
    # when it is likely to matter or to be confusing.
    #
    config.repair.cosmicray.keepCRs = True

if hasattr(config, 'trace'):
    config.trace.xLow = -4.5
    config.trace.xHigh = 4.5

if hasattr(config, 'xOffsetHdrKeyWord'):
    config.xOffsetHdrKeyWord = 'W_FCA_DITHER'
