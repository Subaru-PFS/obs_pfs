# PFS configuration for lsst.pipe.tasks.RepairTask
# 'config' is a lsst.pipe.tasks.RepairConfig

# This CR tuning comes from looking at SuNSS b/r data
#
config.interp.modelPsf.defaultFwhm = 1.5
config.cosmicray.nCrPixelMax = 5000000

# Interpolation needs more work. In the meanwhile do not hide what is there
# when it is likely to matter or to be confusing.
#
config.cosmicray.keepCRs = True
