# PFS configuration for lsst.pipe.tasks.RepairTask
# 'config' is a lsst.pipe.tasks.RepairConfig

# This CR tuning comes from looking at simulated r1 images w/o scattering.
#
config.interp.modelPsf.defaultFwhm = 1.75
config.cosmicray.cond3_fac = 4
config.cosmicray.cond3_fac2 = 1
config.cosmicray.nCrPixelMax = 5000000
config.cosmicray.minSigma = 10.0
config.cosmicray.min_DN = 500.0

# Interpolation needs more work. In the meanwhile do not hide what is there
# when it is likely to matter or to be confusing.
#
config.cosmicray.keepCRs = True
