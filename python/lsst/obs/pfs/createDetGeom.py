import os

import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
from lsst.pipe.base import Task
import lsst.utils

class CreateDetGeomConfig(pexConfig.Config):
    """Configuration for creating the PFS detector geometry files

    No changes required compared to the base class, but
    subclassed for distinction.
    """
    pass

class CreateDetGeomTask(Task):
    """Task to create the detector geometry files for PFS"""
    ConfigClass = CreateDetGeomConfig
    _DefaultName = "createDetGeomTask"

    def run(self):
        # Gain and RdNoise from Jim Gunn's talk in Marseille Dec 2015
        gain = (1.24, 1.24, 1.27, 1.18, 1.26, 1.20, 1.24, 1.26)
        rdnoise = (3.61, 3.78, 3.18, 2.95, 3.19, 3.80, 4.51, 3.18)

        # According to RHL's data model (https://github.com/Subaru-PFS/datamodel/blob/master/datamodel.txt):
        # The raw CCD images will have:
        #    NAMP=8
        #    AMPCOLS=520
        #    CCDROWS=4224
        #
        #    OVERROWS = 76 # as of 2016-02-01
        #    OVERCOLS = 32 # as of 2016-02-01
        #
        #    LEADINROWS = 48 # necked rows
        #    LEADINCOLS  = 8 # real leadin pixels
        #
        # So the raw CCD data will have
        #    NROW = 4300 #  CCDROWS + OVERROWS
        #    NCOL = 4416 # NAMP*(AMPCOLS + OVERCOLS)
        #
        # And the ISR-extracted CCD images will probably be close to:
        #    NROW = 4176  # CCDROWS - LEADINROWS
        #    NCOL = 4096  # NAMP*(AMPCOLS - LEADINCOLS)
        # This value of NROW is applicable to the pfsArm files.

        #                        according to Craig Loomis:
        nColsAmp = 512          #512
        nColsHPrescan = 8       #8 + 1 for first amp, then 8
        nColsHOverscan = 32     #32 - 1 for last amp, otherwise 32
        nRowsVPrescan = 48      #48-51
        nRowsVOverscan = 76     #76
        nRowsAmp = 4176         #actually 4224, ignore first 48 leadin rows
        saturationLevel = 60000

        outDir = os.path.join(lsst.utils.getPackageDir('obs_pfs'), 'pfs/camera')

        amp = 0
        for iArm in range(3):
            for iCCD in range(4):
                detector = afwTable.AmpInfoCatalog(afwTable.AmpInfoTable.makeMinimalSchema())

                for iAmp in range(8):

                    x0 = iAmp * (nColsAmp + nColsHPrescan + nColsHOverscan)
                    x1 = ((iAmp + 1) * (nColsAmp + nColsHPrescan + nColsHOverscan)) - 1
                    y0 = 0
                    y1 = nRowsAmp + nRowsVPrescan + nRowsVOverscan - 1
                    AmpA = detector.addNew()
                    AmpA.setName(str(amp))

                    bb1 = afwGeom.Point2I(iAmp * nColsAmp, 0)
                    bb2 = afwGeom.Point2I((nColsAmp*(iAmp+1))-1, nRowsAmp-1)
                    AmpA.setBBox(afwGeom.Box2I(bb1, bb2))

                    AmpA.setGain(gain[iAmp])
                    AmpA.setReadNoise(rdnoise[iAmp])
                    AmpA.setSaturation(saturationLevel)

                    rbb1 = afwGeom.Point2I(x0, y0)
                    rbb2 = afwGeom.Point2I(x1, y1)
                    AmpA.setRawBBox(afwGeom.Box2I(rbb1, rbb2))

                    rdbb1 = afwGeom.Point2I(x0 + nColsHPrescan, nRowsVPrescan)
                    rdbb2 = afwGeom.Point2I(x1 - nColsHOverscan, y1 - nRowsVOverscan)
                    AmpA.setRawDataBBox(afwGeom.Box2I(rdbb1, rdbb2))

                    if iAmp in [1,3,5,7]:
                        AmpA.setRawFlipX(True)
                    else:
                        AmpA.setRawFlipX(False)
                    AmpA.setRawFlipY(False)
                    AmpA.setRawXYOffset(afwGeom.Extent2I(0, 0))

                    rhobb1 = afwGeom.Point2I(x0 + nColsHPrescan + nColsAmp, nRowsVPrescan)
                    rhobb2 = afwGeom.Point2I(x1, y1 - nRowsVOverscan)
                    AmpA.setRawHorizontalOverscanBBox(afwGeom.Box2I(rhobb1, rhobb2))

                    rvobb1 = afwGeom.Point2I(x0 + nColsHPrescan, y1 - nRowsVOverscan + 1)
                    rvobb2 = afwGeom.Point2I(x1 - nColsHOverscan, y1)
                    AmpA.setRawVerticalOverscanBBox(afwGeom.Box2I(rvobb1, rvobb2))

                    rpbb1 = afwGeom.Point2I(x0 + nColsHPrescan, 0)
                    rpbb2 = afwGeom.Point2I(x1 - nColsHOverscan, nRowsVPrescan - 1)
                    AmpA.setRawPrescanBBox(afwGeom.Box2I(rpbb1, rpbb2))

                    AmpA.setHasRawInfo(True)
                    AmpA.setReadoutCorner(afwTable.LL)

                    amp = amp+1

                sArm = "brnm"[iArm]
                fitsName = os.path.join(outDir, sArm+"_"+str(iCCD+1)+".fits")
                detector.writeFits(fitsName)
