# PFS configuration for lsst.ip.isr.IsrTask
# 'config' is a lsst.ip.isr.IsrConfig

config.expectWcs = False            # our spectrographs don't write a WCS to the header
config.doLinearize = False
config.doFringe = False
config.doWrite = False
config.overscan.fitType = "AKIMA_SPLINE"
config.overscan.order = 30
