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
    doRepair = Field(
        dtype=bool, 
        default=True, 
        doc="Repair artifacts?")
    psfFwhm = Field(
        dtype=float, 
        default=3.0, 
        doc="Repair PSF FWHM (pixels)",
        check = lambda x : x > 0)
    psfSize = Field(
        dtype=int, 
        default=21, 
        doc="Repair PSF size (pixels)",
        check = lambda x : x > 0)
    crGrow = Field(
        dtype=int, 
        default=2, 
        doc="Grow radius for CR (pixels)",
        check = lambda x : x >= 0)
    repair = ConfigurableField(
        target=RepairTask, 
        doc="Task to repair artifacts")
    darkTime = Field(
        dtype=str, 
        default="DARKTIME", 
        doc="Header keyword for time since last CCD wipe, or None",
        optional=True)
    display = Field(
        dtype=bool, 
        default=True, 
        doc="Display outputs?")
        
    """Parameters for the profile determination"""
    profileInterpolation = Field(
        dtype=str, 
        default='SPLINE3', 
        doc="Method for determining the spatial profile, [PISKUNOV, SPLINE3]")
    ccdReadOutNoise = Field(
        dtype=float, 
        default=3.19, 
        doc="CCD readout noise",
        check = lambda x : x > 0)
    maxIterSF = Field(
        dtype=int, 
        default=15, 
        doc="profileInterpolation==PISKUNOV: Maximum number of iterations for the determination of the spatial profile")
    maxIterSig = Field(
        doc = "Maximum number of iterations for masking bad pixels and CCD defects (default: 2)",
        dtype = int,
        default = 2,
        check = lambda x : x >= 0)
    overSample = Field(
        dtype=int, 
        default=15, 
        doc="Oversampling factor for the determination of the spatial profile",
        check = lambda x : x > 0)
    swathWidth = Field(
        dtype=int, 
        default=250, 
        doc="Size of individual extraction swaths, set to 0 to calculate automatically",
        check = lambda x : x > 0)
    wingSmoothFactor = Field(
        dtype=float, 
        default=0., 
        doc="profileInterpolation==PISKUNOV: Lambda smoothing factor to remove possible oscillation of the wings of the spatial profile",
        check = lambda x : x >= 0)
    lowerSigma = Field(
        dtype = float, 
        doc = "lower sigma rejection threshold if maxIterSig > 0 (default: 3.)", 
        default = 3.,
        check = lambda x : x >= 0 )
    upperSigma = Field(
        dtype = float,
        doc = "upper sigma rejection threshold if maxIterSig > 0 (default: 3.)",
        default = 3.,
        check = lambda x : x >= 0 )
        
    """parameters for tracing the apertures"""
    interpolation = Field(
        doc = "Interpolation schemes (CHEBYSHEV, LEGENDRE, CUBIC, LINEAR, POLYNOMIAL[only one implemented atm])",
        dtype = str,
        default = "POLYNOMIAL")
    order = Field(
        doc = "Polynomial order",
        dtype = int,
        default = 6,
        check = lambda x : x >= 0)
    xLow = Field(
        doc = "Lower (left) limit of aperture relative to center position of trace in x (< 0.)",
        dtype = float,
        default = -4.,
        check = lambda x : x < 0.)
    xHigh = Field(
        doc = "Upper (right) limit of aperture relative to center position of trace in x",
        dtype = float,
        default = 4.,
        check = lambda x : x > 0.)
    apertureFWHM = Field(
        doc = "FWHM of an assumed Gaussian spatial profile for tracing the spectra",
        dtype = float,
        default = 2.5,
        check = lambda x : x > 0.)
    signalThreshold = Field(
        doc = "Signal below this threshold is assumed zero for tracing the spectra",
        dtype = float,
        default = 120.,
        check = lambda x : x >= 0.)
    nTermsGaussFit = Field(
        doc = "1 to look for maximum only without GaussFit; 3 to fit Gaussian; 4 to fit Gaussian plus constant background, 5 to fit Gaussian plus linear term (sloped background)",
        dtype = int,
        default = 3,
        check = lambda x : x > 0)
    saturationLevel = Field(
        doc = "CCD saturation level",
        dtype = float,
        default = 65000.,
        check = lambda x : x > 0.)
    minLength = Field(
        doc = "Minimum aperture length to count as found FiberTrace",
        dtype = int,
        default = 3000,
        check = lambda x : x >= 0)
    maxLength = Field(
        doc = "Maximum aperture length to count as found FiberTrace",
        dtype = int,
        default = 4096,
        check = lambda x : x >= 0)
    nLost = Field(
        doc = "Number of consecutive times the trace is lost before aborting the trace",
        dtype = int,
        default = 10,
        check = lambda x : x >= 0)

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
        
        ignoreAps = [1,2,5,6]

#        flatVisits = []#29,41,42,44,45,46,47,48,49,51,53]
        
        myFindTask = fataTask.FindAndTraceAperturesTask()
        myFindTask.config.interpolation = self.config.interpolation
        myFindTask.config.order = self.config.order
        myFindTask.config.xLow = self.config.xLow
        myFindTask.config.xHigh = self.config.xHigh
        myFindTask.config.apertureFWHM = self.config.apertureFWHM
        myFindTask.config.signalThreshold = self.config.signalThreshold
        myFindTask.config.nTermsGaussFit = self.config.nTermsGaussFit
        myFindTask.config.saturationLevel = self.config.saturationLevel
        myFindTask.config.minLength = self.config.minLength
        myFindTask.config.maxLength = self.config.maxLength
        myFindTask.config.nLost = self.config.nLost
        
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
        myProfileTask.config.lowerSigma = self.config.lowerSigma
        myProfileTask.config.upperSigma = self.config.upperSigma
        myProfileTask.config.maxIterSig = self.config.maxIterSig

        for fts in allFts:
            """Remove next 2 lines when tested"""
#            for iFt in range(fts.size()):
#                fts.getFiberTrace(iFt).getImage().getArray()[:,:] += 10.
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
                if ft not in ignoreAps:
                    spectrum = ft.extractFromProfile()
                    recFt = ft.getReconstructed2DSpectrum(spectrum)
                    recFtArr = recFt.getArray()
                    imArr = ft.getImage().getArray()
                    recFtArr = drpStella.where(imArr,'<=',0.,0.,recFtArr)
                    drpStella.addFiberTraceToCcdImage(ft, recFt, sumRecIm)
                    drpStella.addFiberTraceToCcdImage(ft, ft.getVariance(), sumVarIm)

        sumVariances = drpStella.where(sumVariances, '<', 0., 0., sumVariances)
        snrArr = sumFlats / np.sqrt(sumVariances)

        normalizedFlat = sumFlats / sumRecIm.getArray()
        normalizedFlat = drpStella.where(sumRecIm.getArray(), '<=', 0., 1., normalizedFlat)
        #normalizedFlat = drpStella.where(snrArr, '<', 100., 1., normalizedFlat)
        normalizedFlat = drpStella.where(sumFlats, '<=', 0., 1., normalizedFlat)
        normalizedFlat = drpStella.where(sumVariances, '<=', 0., 1., normalizedFlat)
        
        normFlatOut = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(normalizedFlat)))
        if self.config.display:
            if afwDisplay:
                afwDisplay.ds9.mtv(normFlatOut, title="normalized Flat out", frame=1)
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
                meanArrs = []
                for iFt in range( allFts[ i ].size() ):
                    if iFt not in ignoreAps:
                        ftZero = drpStella.FiberTraceF( normFlatOut.getMaskedImage(), allFts[ i ].getFiberTrace( iFt ).getFiberTraceFunction() )
                        ftsZeroOffset.addFiberTrace( ftZero );
                        drpStella.addFiberTraceToCcdImage( ftZero, ftZero.getImage(), normIm )
                        meanArr = np.mean(ftZero.getImage().getArray(), 0)
                        stdDevArr = np.std(ftZero.getImage().getArray(), 0)
                        meanArrs.append( meanArr )
                        print 'FiberTrace ',iFt,': meanArr = ',meanArr
                        print 'FiberTrace ',iFt,': stdDevArr = ',stdDevArr

                        ftSnr = drpStella.FiberTraceF( snrMI, allFts[ i ].getFiberTrace( iFt ).getFiberTraceFunction() )
                        ftsSnr.addFiberTrace( ftSnr );
                        drpStella.addFiberTraceToCcdImage( ftSnr, ftSnr.getImage(), snrIm )

                        ftSumFlats = drpStella.FiberTraceF( sumFlatsMI, allFts[ i ].getFiberTrace( iFt ).getFiberTraceFunction() )
                        ftsSumFlats.addFiberTrace( ftSumFlats );
                        drpStella.addFiberTraceToCcdImage( ftSumFlats, ftSumFlats.getImage(), sumIm )

                        ftRecFlats = drpStella.FiberTraceF( recFlatsMI, allFts[ i ].getFiberTrace( iFt ).getFiberTraceFunction() )
                        ftsRecFlats.addFiberTrace( ftRecFlats );
                        drpStella.addFiberTraceToCcdImage( ftRecFlats, ftRecFlats.getImage(), recIm )
                meanVec = np.mean( meanArrs, 0 )
                print 'meanVec = ',meanVec

        if ftsZeroOffset.size() < 1:
            raise RunTimeError("constructNormFlatTask: ERROR: no Flat found with 0 xOffset")
        print 'ftsZeroOffset.size() = ',ftsZeroOffset.size()
        
        if self.config.display:
            if afwDisplay:
                afwDisplay.ds9.mtv(sumIm, title="sum", frame=2)
                afwDisplay.ds9.mtv(recIm, title="reconstructed", frame=3)
                sumMinusRec = afwImage.ImageF(sumIm.getArray() - recIm.getArray())
                afwDisplay.ds9.mtv(sumMinusRec, title="sum - reconstructed", frame=4)
                afwDisplay.ds9.mtv(snrIm, title="SNR", frame=5)
                afwDisplay.ds9.mtv(normIm, title="normalizedFlat", frame=6)
        
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
