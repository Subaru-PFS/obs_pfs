#!/usr/bin/env python
import math
import lsst.afw.image as afwImage
import numpy as np
import pfs.drp.stella as drpStella
import pfs.drp.stella.findAndTraceAperturesTask as fataTask
import pfs.drp.stella.createFlatFiberTraceProfileTask as cfftpTask
#from lsst.utils import getPackageDir
from lsst.pipe.drivers.constructCalibs import CalibConfig, CalibTask
from lsst.pex.config import Field, ConfigurableField
from lsst.pipe.tasks.repair import RepairTask
import lsst.meas.algorithms as measAlg
import lsst.afw.detection as afwDet
from lsst.pipe.drivers.utils import getDataRef
from lsst.ctrl.pool.pool import NODE
try:
    import lsst.afw.display as afwDisplay
except ImportError:
    afwDisplay = None

if afwDisplay:
    try:
        afwDisplay.setDefaultBackend("ds9" if True else "virtualDevice")
    except RuntimeError as e:
        print e

class ConstructNormFlatConfig(CalibConfig):
    """Configuration for flat construction"""
    doRepair = Field(dtype=bool, default=True, doc="Repair artifacts?")
    psfFwhm = Field(dtype=float, default=3.0, doc="Repair PSF FWHM (pixels)")
    psfSize = Field(dtype=int, default=21, doc="Repair PSF size (pixels)")
    crGrow = Field(dtype=int, default=2, doc="Grow radius for CR (pixels)")
    repair = ConfigurableField(target=RepairTask, doc="Task to repair artifacts")
    darkTime = Field(dtype=str, default="DARKTIME", doc="Header keyword for time since last CCD wipe, or None",
                     optional=True)
    display = Field(dtype=bool, default=True, doc="Display outputs?")
    profileInterpolation = Field(dtype=str, default='SPLINE3', doc="Method for determining the spatial profile, [PISKUNOV, SPLINE3]")
    ccdReadOutNoise = Field(dtype=float, default=3.19, doc="CCD readout noise")
    maxIterSF = Field(dtype=int, default=15, doc="profileInterpolation==PISKUNOV: Maximum number of iterations for the determination of the spatial profile")
    overSample = Field(dtype=int, default=15, doc="Oversampling factor for the determination of the spatial profile")
    swathWidth = Field(dtype=int, default=250, doc="Size of individual extraction swaths, set to 0 to calculate automatically")
    wingSmoothFactor = Field(dtype=float, default=0., doc="profileInterpolation==PISKUNOV: Lambda smoothing factor to remove possible oscillation of the wings of the spatial profile")

#class ConstructNormFlatConfig(Config):
#    """Configuration for reducing arc images"""
#    function = Field( doc = "Function for fitting the dispersion", dtype=str, default="POLYNOMIAL" );
#    order = Field( doc = "Fitting function order", dtype=int, default = 5 );
#    searchRadius = Field( doc = "Radius in pixels relative to line list to search for emission line peak", dtype = int, default = 2 );
#    fwhm = Field( doc = "FWHM of emission lines", dtype=float, default = 2.6 );
#    nRowsPrescan = Field( doc = "Number of prescan rows in raw CCD image", dtype=int, default = 49 );
#    wavelengthFile = Field( doc = "reference pixel-wavelength file including path", dtype = str, default=os.path.join(getPackageDir("obs_pfs"), "pfs/RedFiberPixels.fits.gz"));
#    lineList = Field( doc = "reference line list including path", dtype = str, default=os.path.join(getPackageDir("obs_pfs"), "pfs/lineLists/CdHgKrNeXe_red.fits"));

#class ConstructNormFlatTaskRunner(TaskRunner):
#    """Get parsed values into the ConstructNormFlatTask.run"""
#    @staticmethod
#    def getTargetList(parsedCmd, **kwargs):
#        print 'ConstructNormFlatTask.getTargetList: kwargs = ',kwargs
#        return [dict(expRefList=parsedCmd.id.refList, butler=parsedCmd.butler)]#, wLenFile=parsedCmd.wLenFile, lineList=parsedCmd.lineList)]
#
#    def __call__(self, args):
#        task = self.TaskClass(config=self.config, log=self.log)
#        if self.doRaise:
#            self.log.info('ConstructNormFlatTask.__call__: args = %s' % args)
#            result = task.run(**args)
#        else:
#            try:
#                result = task.run(**args)
#            except Exception, e:
#                task.log.fatal("Failed: %s" % e)
##                traceback.print_exc(file=sys.stderr)
#
#        if self.doReturnResults:
#            return Struct(
#                args = args,
#                metadata = task.metadata,
#                result = result,
#            )

class ConstructNormFlatTask(CalibTask):
    """Task to construct the normalized Flat"""
    ConfigClass = ConstructNormFlatConfig
    _DefaultName = "constructNormFlat"
    calibName = "flat"
#    filterName = "NONE"  # Sets this filter name in the output

    def __init__(self, *args, **kwargs):
        CalibTask.__init__(self, *args, **kwargs)
        self.makeSubtask("repair")

    @classmethod
    def applyOverrides(cls, config):
        """Overrides to apply for dark construction"""
        config.isr.doDark = False
        config.isr.doFlat = False
        config.isr.doFringe = False

    def processSingle(self, sensorRef):
        """Process a single CCD

        Besides the regular ISR, also masks cosmic-rays.
        """
        exposure = CalibTask.processSingle(self, sensorRef)

        if self.config.doRepair:
            psf = measAlg.DoubleGaussianPsf(self.config.psfSize, self.config.psfSize,
                                            self.config.psfFwhm/(2*math.sqrt(2*math.log(2))))
            exposure.setPsf(psf)
            self.repair.run(exposure, keepCRs=False)
            if self.config.crGrow > 0:
                mask = exposure.getMaskedImage().getMask().clone()
                mask &= mask.getPlaneBitMask("CR")
                fpSet = afwDet.FootprintSet(mask.convertU(), afwDet.Threshold(0.5))
                fpSet = afwDet.FootprintSet(fpSet, self.config.crGrow, True)
                fpSet.setMask(exposure.getMaskedImage().getMask(), "CR")
        return exposure

    def combine(self, cache, struct):
        """!Combine multiple exposures of a particular CCD and write the output

        Only the slave nodes execute this method.

        @param cache  Process pool cache
        @param struct  Parameters for the combination, which has the following components:
            * ccdIdList   List of data identifiers for combination
            * scales      Scales to apply (expScales are scalings for each exposure,
                               ccdScale is final scale for combined image)
            * outputId    Data identifier for combined image (fully qualified for this CCD)
        """
        dataRefList = [getDataRef(cache.butler, dataId) if dataId is not None else None for
                       dataId in struct.ccdIdList]
        self.log.info("Combining %s on %s" % (struct.outputId, NODE))

        self.log.info('len(dataRefList) = %d' % len(dataRefList))

#        flatVisits = []#29,41,42,44,45,46,47,48,49,51,53]
        
        myFindTask = fataTask.FindAndTraceAperturesTask()
        
        exposure = dataRefList[0].get('postISRCCD')
        
        sumFlats = exposure.getMaskedImage().getImage().getArray()
        sumVariances = exposure.getMaskedImage().getVariance().getArray()

        allFts = []
        xOffsets = []
        for expRef in dataRefList:
            exposure = expRef.get('postISRCCD')
            xOffsets.append(exposure.getMetadata().get('sim.slit.xoffset'))
            fts = myFindTask.run(exposure)
            allFts.append(fts)
            if expRef.dataId['visit'] != dataRefList[0].dataId['visit']:
                sumFlats += exposure.getMaskedImage().getImage().getArray()
                sumVariances += exposure.getMaskedImage().getVariance().getArray()

        self.log.info('=== xOffsets = '+str(xOffsets)+' ===')

        myProfileTask = cfftpTask.CreateFlatFiberTraceProfileTask()
        myProfileTask.config.profileInterpolation = self.config.profileInterpolation
        myProfileTask.config.ccdReadOutNoise = self.config.ccdReadOutNoise
        myProfileTask.config.maxIterSF = self.config.maxIterSF
        myProfileTask.config.overSample = self.config.overSample
        myProfileTask.config.swathWidth = self.config.swathWidth
        myProfileTask.config.lambdaSF = 1./float(self.config.overSample)
        myProfileTask.config.lambdaSP = 0.
        myProfileTask.config.wingSmoothFactor = self.config.wingSmoothFactor

        for fts in allFts:
            fts = myProfileTask.run(fts)

        sumRec = np.ndarray(shape=sumFlats.shape, dtype='float32')
        sumRec[:][:] = 0.
        sumRecIm = afwImage.ImageF(sumRec)
        rec = np.ndarray(shape=sumFlats.shape, dtype='float32')
        rec[:][:] = 0.
        recIm = afwImage.ImageF(rec)

        sumVar = np.ndarray(shape=sumFlats.shape, dtype='float32')
        sumVar[:][:] = 0.
        sumVarIm = afwImage.ImageF(sumVar)

        # Add all reconstructed FiberTraces of all dithered flats to one reconstructed image 'recIm'
        for fts in allFts:
            recIm.getArray()[:][:] = 0.
            for ft in fts.getTraces():
                spectrum = ft.extractFromProfile()
                recFt = ft.getReconstructed2DSpectrum(spectrum)
                recFtArr = recFt.getArray()
                imArr = ft.getImage().getArray()
                recFtArr = drpStella.where(imArr,'<=',0.,0.,recFtArr)
                drpStella.addFiberTraceToCcdImage(ft, recFt, sumRecIm)
                drpStella.addFiberTraceToCcdImage(ft, ft.getVariance(), sumVarIm)

        sumVariances = drpStella.where(sumVariances, '<', 0., 0., sumVariances)
        snrArr = sumFlats / np.sqrt(sumVariances)

        normalizedFlat = sumRecIm.getArray() / sumFlats
        normalizedFlat = drpStella.where(sumRecIm.getArray(), '<=', 0., 1., normalizedFlat)
        normalizedFlat = drpStella.where(snrArr, '<', 100., 1., normalizedFlat)
        normalizedFlat = drpStella.where(sumFlats, '<=', 0., 1., normalizedFlat)
        normalizedFlat = drpStella.where(sumVariances, '<=', 0., 1., normalizedFlat)
        
        normFlatOut = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(normalizedFlat)))
#        normFlatOut.getMaskedImage().getVariance().getArray()[:][:] = snrArr[:][:]
#        print 'dir(normFlatOut) = ',dir(normFlatOut)
#        normFlatOut.getMaskedImage().getMask().getArray()[:][:] = dataRefList[0].get('postISRCCD').getMaskedImage().getMask().getArray()[:][:]
#        print 'type(normFlatOut) = ',type(normFlatOut)

        self.recordCalibInputs(cache.butler, normFlatOut, struct.ccdIdList, struct.outputId)

        self.interpolateNans(normFlatOut)
        
        normFlatOutDec = afwImage.DecoratedImageF(normFlatOut.getMaskedImage().getImage())

        self.write(cache.butler, normFlatOutDec, struct.outputId)
        
        """Quality assessment"""
        snrMI = afwImage.makeMaskedImage( afwImage.ImageF( snrArr ) )
        sumFlatsMI = afwImage.makeMaskedImage( afwImage.ImageF( sumFlats ) )
        recFlatsMI = afwImage.makeMaskedImage( sumRecIm )
        
        ftsZeroOffset = drpStella.FiberTraceSetF()
        ftsSnr = drpStella.FiberTraceSetF()
        ftsSumFlats = drpStella.FiberTraceSetF()
        ftsRecFlats = drpStella.FiberTraceSetF()
        
        normArr = np.ndarray(shape=sumFlats.shape, dtype='float32')
        normArr[:][:] = 0.
        normIm = afwImage.ImageF(normArr)
        sumArr = np.ndarray(shape=sumFlats.shape, dtype='float32')
        sumArr[:][:] = 0.
        sumIm = afwImage.ImageF(sumArr)
        snrArr = np.ndarray(shape=sumFlats.shape, dtype='float32')
        snrArr[:][:] = 0.
        snrIm = afwImage.ImageF(snrArr)
        recArr = np.ndarray(shape=sumFlats.shape, dtype='float32')
        recArr[:][:] = 0.
        recIm = afwImage.ImageF(recArr)
        
        for i in range( len( xOffsets ) ):
            if xOffsets[ i ] == 0.:
                for iFt in range( allFts[ i ].size() ):
                    ftZero = drpStella.FiberTraceF( normFlatOut.getMaskedImage(), allFts[ i ].getFiberTrace( iFt ).getFiberTraceFunction() )
                    ftsZeroOffset.addFiberTrace( ftZero );
                    drpStella.addFiberTraceToCcdImage( ftZero, ftZero.getImage(), normIm )
                    
                    ftSnr = drpStella.FiberTraceF( snrMI, allFts[ i ].getFiberTrace( iFt ).getFiberTraceFunction() )
                    ftsSnr.addFiberTrace( ftSnr );
                    drpStella.addFiberTraceToCcdImage( ftSnr, ftSnr.getImage(), snrIm )
                    
                    ftSumFlats = drpStella.FiberTraceF( sumFlatsMI, allFts[ i ].getFiberTrace( iFt ).getFiberTraceFunction() )
                    ftsSumFlats.addFiberTrace( ftSumFlats );
                    drpStella.addFiberTraceToCcdImage( ftSumFlats, ftSumFlats.getImage(), sumIm )

                    ftRecFlats = drpStella.FiberTraceF( recFlatsMI, allFts[ i ].getFiberTrace( iFt ).getFiberTraceFunction() )
                    ftsRecFlats.addFiberTrace( ftRecFlats );
                    drpStella.addFiberTraceToCcdImage( ftRecFlats, ftRecFlats.getImage(), recIm )

        if ftsZeroOffset.size() != allFts[ 0 ].size():
            raise RunTimeError("constructNormFlatTask: ERROR: no Flat found with 0 xOffset")
        print 'ftsZeroOffset.size() = ',ftsZeroOffset.size()
        
        if self.config.display:
            if afwDisplay:
                afwDisplay.ds9.mtv(sumIm, title="sum", frame=0)
                afwDisplay.ds9.mtv(recIm, title="reconstructed", frame=1)
                afwDisplay.ds9.mtv(snrIm, title="SNR", frame=2)
                afwDisplay.ds9.mtv(normIm, title="normalizedFlat", frame=3)
        
#        for ft in ftsZeroOffset:
#            for x in range( ft.getWidth() ):
                
        
        if False:
            meanSDev = np.ndarray(shape=(normalizedFlat.shape[0],2), dtype=np.float)
            print 'meanSDev.shape = ',meanSDev.shape
            for x in range(normalizedFlat.shape[1]):
                ind = drpStella.getIndices(normalizedFlat[:][x])
                print 'x=',x,': ind = ',ind
                nInd = ind.shape[0]
                nPix = normalizedFlat.shape[0] - nInd
                print 'x=',x,': nPix = ',nPix
                pixArr = np.ndarray(shape=(nPix), dtype=np.float)
                index = 0
                for y in range(normalizedFlat.shape[0]):
                    if y not in ind:
                        print 'y = ',y,' not in ind = ',ind
                        pixArr[index] = normalizedFlat[y][x]
                        index += 1
                meanSDev[x][0] = np.mean(pixArr)
                meanSDev[x][1] = np.std(pixArr)
                print 'x=',x,': meanSDev[',x,'][:] = ',meanSDev[x][:]

        if False:
            allFtsSorted = []
            allReconstructedFtsSorted = []
            for iAp in range( allFts[ 0 ].size() ):
                ftsSorted = drpStella.FiberTraceSetF()
                reconstructedFtsSorted = drpStella.FiberTraceSetF()
                for iFts in range( len( allFts ) ):
                    ftsSorted.addFiberTrace( allFts[ iFts ].getFiberTrace( iAp ) )
                    ft = drpStella.FiberTraceF( allFts[ iFts ].getFiberTrace( iAp ) )
                    recIm = ft.getReconstructed2DSpectrum(ft.extractFromProfile())
                    ft.setImage(recIm)
                    reconstructedFtsSorted.addFiberTrace( ft )
                    print 'FiberTrace ',iAp,' from allFts[',iFts,'] added to ftsSorted -> ftsSorted.size() = ',ftsSorted.size()
                allFtsSorted.append( ftsSorted )
                allReconstructedFtsSorted.append( reconstructedFtsSorted )
                print 'ftsSorted added to allFtsSorted => len(allFtsSorted) = ',len( allFtsSorted )

            allSumsFtsSorted = []
            allSumsReconstructedFtsSorted = []
            normFlatTraces = []
            for iFts in range( len( allFtsSorted ) ):
                sumFtsSorted = drpStella.addFiberTraces( allFtsSorted[ iFts ] )
                allSumsFtsSorted.append( sumFtsSorted )
                print 'sumFtsSorted = ',sumFtsSorted.getImage().getArray().shape,': ',sumFtsSorted,' added to allSumsFtsSorted => len(allSumsFtsSorted) = ',len(allSumsFtsSorted)
                sumReconstructedFtsSorted = drpStella.addFiberTraces( allReconstructedFtsSorted[ iFts ] )
                allSumsReconstructedFtsSorted.append( sumReconstructedFtsSorted )

                snr = sumFtsSorted.getImage().getArray() / np.sqrt(sumFtsSorted.getVariance().getArray())
                normFlat = sumFtsSorted.getImage().getArray() / sumReconstructedFtsSorted.getImage().getArray()
                normFlat = drpStella.where(snr,'<',100.,1.,normFlat)

                normFlatTraces.append(normFlat)
#                drpStella.addArrayIntoArray( normFlat,
#                                 ndarray::Array< size_t, 2, 1 > const& minMax,
#                                 size_t const& yMin,
#                                 size_t const& yMax,
#                                 ndarray::Array< T, 2, J > & bigArr )
