"""
PFS-specific overrides for ProcessCcdTask
(applied after Subaru overrides in ../processCcd.py).
"""
import os.path

from lsst.utils import getPackageDir

config.load(os.path.join(os.path.join(getPackageDir("obs_pfs"), "config"), 'isr.py'))
#config.calibrate.photocal.colorterms.load(os.path.join(suprimecamConfigDir, 'colorterms.py'))

#config.measurement.algorithms["jacobian"].pixelScale = 0.2

