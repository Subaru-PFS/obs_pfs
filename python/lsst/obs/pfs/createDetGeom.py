# example: createDetGeom.py
import os
import lsst.pex.config                  as pexConfig
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
from lsst.pipe.base import Task
import lsst.utils

class CreateDetGeomConfig(pexConfig.Config):
    """Configuration for creating the PFS detector geometry files"""

class CreateDetGeomTask(Task):
    """Task to create the detector geometry files for PFS"""
    ConfigClass = CreateDetGeomConfig
    _DefaultName = "createDetGeomTask"

    def __init__(self, *args, **kwargs):
        Task.__init__(self, *args, **kwargs)

    def run(self):
        
        """ Gain and RdNoise from Jim Gunn's talk in Marseille Dec 2015"""
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
        #    NROW = 4172  # CCDROWS - LEADINROWS
        #    NCOL = 4096  # NAMP*(AMPCOLS - LEADINCOLS)
        # This value of NROW is applicable to the pfsArm files.

        #                        according to Craig Loomis:
        nColsAmp = 512#          512
        nColsHPrescan = 8#9#     8 + 1 for first, then 8
        nColsHOverscan = 32#31#  32 - 1 for last, otherwise 32
        nRowsVPrescan = 48#51#   48-51
        nRowsVOverscan = 76#75#  76
        nRowsAmp = 4176#4174#    actual 4224, ignore first ~50
        saturationLevel = 60000
        
        outDir = os.path.join(lsst.utils.getPackageDir('obs_pfs'),'pfs/camera')
        
        amp = 0
        for iArm in range(3):
            for iCCD in range(4):
                detector = afwTable.AmpInfoCatalog(afwTable.AmpInfoTable.makeMinimalSchema())

                for iAmp in range(8):

                    x0 = iAmp * (nColsAmp + nColsHPrescan + nColsHOverscan)
                    x1 = ((iAmp + 1) * (nColsAmp + nColsHPrescan + nColsHOverscan)) - 1
                    print 'iAmp = ',iAmp,': x0 = ',x0,', x1 = ',x1,', extend = ',x1-x0+1
                    y0 = 0
                    y1 = nRowsAmp + nRowsVPrescan + nRowsVOverscan - 1
                    print 'iAmp = ',iAmp,': y0 = ',y0,', y1 = ',y1,', extend = ',y1-y0+1
                    AmpA = detector.addNew()
                    AmpA.setName(str(amp))

                    bb1 = afwGeom.Point2I(iAmp * nColsAmp, 0)
                    bb2 = afwGeom.Point2I((nColsAmp*(iAmp+1))-1, nRowsAmp-1)
                    print 'iAmp = ',iAmp,': bb1 = ',bb1
                    print 'iAmp = ',iAmp,': bb2 = ',bb2
                    print 'iAmp = ',iAmp,': bb: extend_x = ',bb2[0]-bb1[0]+1,', extend_y = ',bb2[1] - bb1[1]+1
                    AmpA.setBBox(afwGeom.Box2I(bb1, bb2))

                    AmpA.setGain(gain[iAmp])
                    print 'gain[',iAmp,'] = ',gain[iAmp]
                    print 'rdnoise[',iAmp,'] = ',rdnoise[iAmp]
                    AmpA.setReadNoise(rdnoise[iAmp])
                    AmpA.setSaturation(saturationLevel)

                    rbb1 = afwGeom.Point2I(x0, y0)
                    rbb2 = afwGeom.Point2I(x1, y1)
                    print 'iAmp = ',iAmp,': rbb1 = ',rbb1
                    print 'iAmp = ',iAmp,': rbb2 = ',rbb2
                    print 'iAmp = ',iAmp,': rbb: extend_x = ',rbb2[0]-rbb1[0]+1,', extend_y = ',rbb2[1] - rbb1[1]+1
                    AmpA.setRawBBox(afwGeom.Box2I(rbb1, rbb2))

#                    rdbb1 = afwGeom.Point2I(x0+nColsHPrescan, 49)
#                    rdbb2 = afwGeom.Point2I(x1-nColsVOverscan, y1-77)
                    rdbb1 = afwGeom.Point2I(x0 + nColsHPrescan, nRowsVPrescan)
                    rdbb2 = afwGeom.Point2I(x1 - nColsHOverscan, y1 - nRowsVOverscan)
                    print 'iAmp = ',iAmp,': rdbb1 = ',rdbb1
                    print 'iAmp = ',iAmp,': rdbb2 = ',rdbb2
                    print 'iAmp = ',iAmp,': rdbb: extend_x = ',rdbb2[0]-rdbb1[0]+1,', extend_y = ',rdbb2[1] - rdbb1[1]+1
                    AmpA.setRawDataBBox(afwGeom.Box2I(rdbb1, rdbb2))

                    if iAmp in [1,3,5,7]:
                        AmpA.setRawFlipX(True)
                    else:
                        AmpA.setRawFlipX(False)
                    AmpA.setRawFlipY(False)
                    AmpA.setRawXYOffset(afwGeom.Extent2I(0, 0))

#                    rhobb1 = afwGeom.Point2I(x0+522, 49)
#                    rhobb2 = afwGeom.Point2I(x1-1, y1-77)
                    rhobb1 = afwGeom.Point2I(x0 + nColsHPrescan + nColsAmp, nRowsVPrescan)
                    rhobb2 = afwGeom.Point2I(x1, y1 - nRowsVOverscan)
                    print 'iAmp = ',iAmp,': rhobb1 = ',rhobb1
                    print 'iAmp = ',iAmp,': rhobb2 = ',rhobb2
                    print 'iAmp = ',iAmp,': rhobb: extend_x = ',rhobb2[0]-rhobb1[0]+1,', extend_y = ',rhobb2[1] - rhobb1[1]+1
                    AmpA.setRawHorizontalOverscanBBox(afwGeom.Box2I(rhobb1, rhobb2))

#                    rvobb1 = afwGeom.Point2I(x0+9, y1-76)
#                    rvobb2 = afwGeom.Point2I(x1-31, y1)
                    rvobb1 = afwGeom.Point2I(x0 + nColsHPrescan, y1 - nRowsVOverscan + 1)
                    rvobb2 = afwGeom.Point2I(x1 - nColsHOverscan, y1)
                    print 'iAmp = ',iAmp,': rvobb1 = ',rvobb1
                    print 'iAmp = ',iAmp,': rvobb2 = ',rvobb2
                    print 'iAmp = ',iAmp,': rvobb: extend_x = ',rvobb2[0]-rvobb1[0]+1,', extend_y = ',rvobb2[1] - rvobb1[1]+1
                    AmpA.setRawVerticalOverscanBBox(afwGeom.Box2I(rvobb1, rvobb2))

                    rpbb1 = afwGeom.Point2I(x0 + nColsHPrescan, 0)
                    rpbb2 = afwGeom.Point2I(x1 - nColsHOverscan, nRowsVPrescan - 1)
                    print 'iAmp = ',iAmp,': rpbb1 = ',rpbb1
                    print 'iAmp = ',iAmp,': rpbb2 = ',rpbb2
                    print 'iAmp = ',iAmp,': rpbb: extend_x = ',rpbb2[0]-rpbb1[0]+1,', extend_y = ',rpbb2[1] - rpbb1[1]+1
                    AmpA.setRawPrescanBBox(afwGeom.Box2I(rpbb1, rpbb2))

                    AmpA.setHasRawInfo(True)
                    AmpA.setReadoutCorner(afwTable.LL)

                    amp = amp+1

                sArm = "brnm"[iArm]
                fitsName = os.path.join(outDir,sArm+"_"+str(iCCD+1)+".fits")
                print "fitsName = <",fitsName,">"
                detector.writeFits(fitsName)

