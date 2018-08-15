# PFS configuration for lsst.ip.isr.IsrTask
# 'config' is a lsst.ip.isr.IsrConfig

if False:  # Disabled until we have crosstalk coefficients
    from lsst.obs.subaru.crosstalk import CrosstalkTask
    config.crosstalk.retarget(CrosstalkTask)
    config.doCrosstalk = True

config.expectWcs = False            # our spectrographs don't write a WCS to the header
config.doLinearize = False
config.doFringe = False
config.doWrite = False
config.doAddDistortionModel = False
config.overscanFitType = "AKIMA_SPLINE"
config.overscanOrder = 30
