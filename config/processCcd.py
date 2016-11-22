"""
PFS-specific overrides for ProcessCcdTask
(applied after Subaru overrides in ../processCcd.py).
"""
import os.path

from lsst.utils import getPackageDir

pfsConfigDir = os.path.join(getPackageDir("obs_pfs"), "config", "pfs")
config.load(os.path.join(pfsConfigDir, 'isr.py'))
#config.calibrate.photocal.colorterms.load(os.path.join(suprimecamConfigDir, 'colorterms.py'))

#config.measurement.algorithms["jacobian"].pixelScale = 0.2

