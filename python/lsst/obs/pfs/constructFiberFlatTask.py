#!/usr/bin/env python
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
from lsst.ctrl.pool.pool import NODE
import lsst.meas.algorithms as measAlg
from lsst.pex.config import Field, ConfigurableField
from lsst.pipe.drivers.constructCalibs import CalibConfig, CalibTask
from lsst.pipe.drivers.utils import getDataRef
from lsst.pipe.tasks.repair import RepairTask
import math
import numpy as np
import pfs.drp.stella as drpStella
import pfs.drp.stella.createFlatFiberTraceProfileTask as profileTask
import pfs.drp.stella.findAndTraceAperturesTask as traceTask
from pfs.drp.stella.utils import makeFiberTraceSet

class ConstructFiberFlatConfig(CalibConfig):
    """Configuration for flat construction"""
    xOffsetHdrKeyWord = Field(dtype=str, default='sim.slit.xoffset', doc="Header keyword for fiber offset in input files")
    doRepair = Field(dtype=bool, default=True, doc="Repair artifacts?")
    psfFwhm = Field(dtype=float, default=3.0, doc="Repair PSF FWHM (pixels)")
    psfSize = Field(dtype=int, default=21, doc="Repair PSF size (pixels)")
    crGrow = Field(dtype=int, default=2, doc="Grow radius for CR (pixels)")
    repair = ConfigurableField(target=RepairTask, doc="Task to repair artifacts")
    darkTime = Field(dtype=str, default="DARKTIME", doc="Header keyword for time since last CCD wipe, or None",
                     optional=True)
    display = Field(dtype=bool, default=True, doc="Display images?")
    profileInterpolation = Field(
        doc = "Method for determining the spatial profile, [PISKUNOV, SPLINE3], default: PISKUNOV",
        dtype = str,
        default = "SPLINE3")
    swathWidth = Field(
        doc = "Size of individual extraction swaths",
        dtype = int,
        default = 500,
        check = lambda x : x > 10)
    telluric = Field(
        doc = "Method for determining the background (+sky in case of slit spectra, default: NONE)",
        dtype = str,
        default = "NONE")
    overSample = Field(
        doc = "Oversampling factor for the determination of the spatial profile (default: 10)",
        dtype = int,
        default = 10,
        check = lambda x : x > 0)
    maxIterSF = Field(
        doc = "Maximum number of iterations for the determination of the spatial profile (default: 8)",
        dtype = int,
        default = 8,
        check = lambda x : x > 0)
    maxIterSky = Field(
        doc = "Maximum number of iterations for the determination of the (constant) background/sky (default: 10)",
        dtype = int,
        default = 10,
        check = lambda x : x >= 0)
    maxIterSig = Field(
        doc = "Maximum number of iterations for masking bad pixels and CCD defects (default: 2)",
        dtype = int,
        default = 2,
        check = lambda x : x > 0)
    minSNR = Field(
        doc = "Minimum Signal-to-Noise Ratio for normalized Flat pixels",
        dtype = float,
        default = 100.,
        check = lambda x : x > 0.)
    lambdaSF = Field(
        doc = "Lambda smoothing factor for spatial profile (default: 1. / overSample)",
        dtype = float,
        default = 17000.,
        check = lambda x : x > 0.)
    lambdaSP = Field(
        doc = "Lambda smoothing factor for spectrum (default: 0)",
        dtype = float,
        default = 0.,
        check = lambda x : x >= 0)
    wingSmoothFactor = Field(
        doc = "Lambda smoothing factor to remove possible oscillation of the wings of the spatial profile (default: 0.)",
        dtype = float,
        default = 0.,
        check = lambda x : x >= 0)
    xLow = Field(
        doc = "Lower (left) limit of aperture relative to center position of trace in x (< 0.)",
        dtype = float,
        default = -5.)
    xHigh = Field(
        doc = "Upper (right) limit of aperture relative to center position of trace in x",
        dtype = float,
        default = 5.,
        check = lambda x : x > 0.)

class ConstructFiberFlatTask(CalibTask):
    """Task to construct the normalized Flat"""
    ConfigClass = ConstructFiberFlatConfig
    _DefaultName = "constructFiberFlat"
    calibName = "flat"

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

        myFindTask = traceTask.FindAndTraceAperturesTask()
        myFindTask.config.xLow = -5.
        myFindTask.config.xHigh = 5.

        exposure = dataRefList[0].get('postISRCCD')
        bias = dataRefList[0].get('bias')
        fiberTrace = dataRefList[0].get('fiberTrace')

        flatsShape = exposure.getMaskedImage().getImage().getArray().shape
        normFlatOut = np.zeros(shape = flatsShape, dtype = np.float32)
        normFlatOutIm = afwImage.ImageF(normFlatOut)

        allFts = []
        xOffsets = []
        for expRef in dataRefList:
            exposure = expRef.get('postISRCCD')
            try:
                xOffsets.append(exposure.getMetadata().get(self.config.xOffsetHdrKeyWord))
            except Exception:
                xOffsets.append(0.0)
            fts = myFindTask.run(exposure)
            allFts.append(fts)

        self.log.info('=== xOffsets = '+str(xOffsets)+' ===')

        myProfileTask = profileTask.CreateFlatFiberTraceProfileTask()
        myProfileTask.config.profileInterpolation = self.config.profileInterpolation
        myProfileTask.config.overSample = self.config.overSample
        myProfileTask.config.swathWidth = self.config.swathWidth
        myProfileTask.config.telluric = self.config.telluric
        myProfileTask.config.maxIterSF = self.config.maxIterSF
        myProfileTask.config.maxIterSky = self.config.maxIterSky
        myProfileTask.config.maxIterSig = self.config.maxIterSig
        myProfileTask.config.lambdaSF = self.config.lambdaSF
        myProfileTask.config.lambdaSP = self.config.lambdaSP
        myProfileTask.config.wingSmoothFactor = self.config.wingSmoothFactor

        sumOrig = np.zeros(shape=flatsShape, dtype=np.float32)
        normFlatTemp = np.zeros(shape=flatsShape, dtype=np.float32)
        normFlatTempIm = afwImage.ImageF(normFlatTemp)
        normFlatTempMIm = afwImage.makeMaskedImage(normFlatTempIm)
        sumVar = np.zeros(shape=flatsShape, dtype=np.float32)
        sumVarIm = afwImage.ImageF(sumVar)

        frame = 0
        # For each FiberTrace of all dithered flats add to one image 'sumOrigIm'
        for iFt in range(allFts[0].size()):
            xMinCenterRow = 100000.
            xMaxCenterRow = 0.;
            sumOrig[:][:] = 0.
            sumOrigIm = afwImage.ImageF(sumOrig)
            sumVar[:][:] = 0.
            sumVarIm = afwImage.ImageF(sumVar)
            for fts in allFts:
                ft = fts.getFiberTrace(iFt)
                drpStella.addFiberTraceToCcdImage(ft, ft.getImage(), sumOrigIm)
                drpStella.addFiberTraceToCcdImage(ft, ft.getVariance(), sumVarIm)

                xCenters = ft.getXCenters()
                xLeft = xCenters[xCenters.shape[0] / 2] + self.config.xLow
                xMinCenterRow = np.min([xMinCenterRow, xLeft])
                xRight = xCenters[xCenters.shape[0] / 2] + self.config.xHigh
                xMaxCenterRow = np.max([xMaxCenterRow, xRight])

            myFindTask = traceTask.FindAndTraceAperturesTask()
            myFindTask.config.xHigh = 0.5 * (xMaxCenterRow - xMinCenterRow)
            myFindTask.config.xLow = -1. * myFindTask.config.xHigh
            myFindTask.config.apertureFWHM = myFindTask.config.xHigh / 2.

            print 'mean(sumOrig) = ',np.mean(sumOrig)
            print 'mean(sumOrigIm.getArray()) = ',np.mean(sumOrigIm.getArray())
            sumOrig = np.where(sumOrig == 0.,
                               bias.getMaskedImage().getImage().getArray(),
                               sumOrig)
            maskedIm = afwImage.makeMaskedImage(sumOrigIm)
            sumVar = sumVarIm.getArray()
            maskedIm.getVariance().getArray()[:,:] = sumVar[:,:]
            exposure = afwImage.makeExposure(maskedIm)

            """find and trace FiberTrace in image containing sum of FiberTraces"""
            fts = myFindTask.run(exposure)

            """calculate profile for sum of FiberTraces"""
            myProfileTask.run(fts)
            ft = fts.getFiberTrace(0)

            if self.config.display:
#                ds9.setMaskTransparency(50)
#                drpStella.markFiberTraceInMask(ft,
#                                               exposure.getMaskedImage().getMask())
                frame += 1
                ds9.mtv(exposure, title='Flat iFT=%d' % iFt, frame=frame)

            imFt = ft.getImage()
            print 'mean(imFt) = ',np.mean(imFt.getArray())
            varFt = ft.getVariance()

            """extract spectrum and reconstruct FiberTrace"""
            specFt = ft.extractFromProfile()
            recFt = ft.getReconstructed2DSpectrum(specFt).getArray()
            print 'mean(recFt) = ',np.mean(recFt)
            recFt = np.where(recFt <= 0., 0.1, recFt)

            """calculate normalized Flat for FiberTrace"""
            normFt = imFt.getArray() / recFt
            normFt = np.where(recFt <= 0., 1., normFt)
            snrArr = imFt.getArray() / np.sqrt(varFt.getArray())
            normFt = np.where(snrArr <= self.config.minSNR, 1., normFt)
            print 'mean(normFt) = ',np.mean(normFt)

            """add normalized FiberTrace to temporary image, extract the"""
            """FiberTraceSet from the pfsFiberTrace"""
            """We do this so we don't add normalized FiberTraces on top of"""
            """each other"""
            """In order to avoid this step we need another function like"""
            """addFiberTraceToCcdImage but with setting the pixel values"""
            """instead of adding them"""
            drpStella.addFiberTraceToCcdImage(ft,
                                              afwImage.ImageF(normFt),
                                              normFlatTempIm)
            print 'after addNormFiberTrace: mean(normFlatTempIm) = ',np.mean(normFlatTempIm.getArray())
            print 'after addNormFiberTrace: mean(normFlatTempMIm) = ',np.mean(normFlatTempMIm.getImage().getArray())
            fts = makeFiberTraceSet(fiberTrace, normFlatTempMIm)
            print 'after makeFiberTraceSet: fts.size() = ',fts.size()
            ftNorm = fts.getFiberTrace(iFt)
            if self.config.display:
                normFlatTempMIm = afwImage.makeMaskedImage(normFlatTempIm)
                drpStella.markFiberTraceInMask(ft,
                                               normFlatTempMIm.getMask())
                frame += 1
                ds9.mtv(normFlatTempMIm, title='normalized Flat iFT=%d' % iFt, frame=frame)
            drpStella.addFiberTraceToCcdImage(ftNorm,
                                              ftNorm.getImage(),
                                              normFlatOutIm)
            print 'iFt = ',iFt,' finished: mean(normFt) = ',np.mean(normFt)

        normFlatOut = normFlatOutIm.getArray()
        normFlatOut = np.where(np.fabs(normFlatOut) < 0.01, 1., normFlatOut)
        normFlatOutIm = afwImage.ImageF(normFlatOut)
        if self.config.display:
            frame += 1
            ds9.mtv(normFlatOutIm, title='normalized Flat', frame=frame)

        """Write fiber flat"""
        normFlatExp = afwImage.makeExposure(afwImage.makeMaskedImage(normFlatOutIm))
        self.recordCalibInputs(cache.butler, normFlatExp, struct.ccdIdList, struct.outputId)
        self.interpolateNans(normFlatExp)
        normFlatOutDec = afwImage.DecoratedImageF(normFlatOutIm)
        self.write(cache.butler, normFlatOutDec, struct.outputId)
