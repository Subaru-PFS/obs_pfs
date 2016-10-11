#!/usr/bin/env python
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
import pfs.drp.stella.createFlatFiberTraceProfileTask as cfftpTask
from pfs.drp.stella.utils import makeFiberTraceSet

class ConstructFiberFlatConfig(CalibConfig):
    """Configuration for flat construction"""
    doRepair = Field(dtype=bool, default=True, doc="Repair artifacts?")
    psfFwhm = Field(dtype=float, default=3.0, doc="Repair PSF FWHM (pixels)")
    psfSize = Field(dtype=int, default=21, doc="Repair PSF size (pixels)")
    crGrow = Field(dtype=int, default=2, doc="Grow radius for CR (pixels)")
    repair = ConfigurableField(target=RepairTask, doc="Task to repair artifacts")
    darkTime = Field(dtype=str, default="DARKTIME", doc="Header keyword for time since last CCD wipe, or None",
                     optional=True)
    
    profileInterpolation = Field(
        doc = "Method for determining the spatial profile, [PISKUNOV, SPLINE3], default: PISKUNOV",
        dtype = str,
        default = "SPLINE3")
    ccdReadOutNoise = Field(
        doc = "CCD readout noise",
        dtype = float,
        default = 1.,
        check = lambda x : x > 0.)
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

#        flatVisits = []#29,41,42,44,45,46,47,48,49,51,53]
        
        exposure = dataRefList[0].get('postISRCCD')
        
        sumFlats = exposure.getMaskedImage().getImage().getArray()
        sumVariances = exposure.getMaskedImage().getVariance().getArray()

        allFts = []
        xOffsets = []
        for expRef in dataRefList:
            exposure = expRef.get('postISRCCD')
            xOffsets.append(exposure.getMetadata().get('sim.slit.xoffset'))
            fiberTrace = expRef.get('fiberTrace', immediate=True)
            fts = makeFiberTraceSet(fiberTrace, exposure.getMaskedImage())
            allFts.append(fts)
            if expRef.dataId['visit'] != dataRefList[0].dataId['visit']:
                sumFlats += exposure.getMaskedImage().getImage().getArray()
                sumVariances += exposure.getMaskedImage().getVariance().getArray()

        self.log.info('=== xOffsets = '+str(xOffsets)+' ===')

        myProfileTask = cfftpTask.CreateFlatFiberTraceProfileTask()
        myProfileTask.config.profileInterpolation = self.config.profileInterpolation
        myProfileTask.config.ccdReadOutNoise = self.config.ccdReadOutNoise
        myProfileTask.config.overSample = self.config.overSample
        myProfileTask.config.swathWidth = self.config.swathWidth
        myProfileTask.config.telluric = self.config.telluric
        myProfileTask.config.maxIterSF = self.config.maxIterSF
        myProfileTask.config.maxIterSky = self.config.maxIterSky
        myProfileTask.config.maxIterSig = self.config.maxIterSig
        myProfileTask.config.lambdaSF = self.config.lambdaSF
        myProfileTask.config.lambdaSP = self.config.lambdaSP
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

        sumVariances = drpStella.where(sumVariances, '<=', 0., 0.1, sumVariances)
        snrArr = sumFlats / np.sqrt(sumVariances)

        sumFlats = drpStella.where(sumFlats, '<=', 0., 0.1, sumFlats)
        normalizedFlat = sumRecIm.getArray() / sumFlats
        normalizedFlat = drpStella.where(sumRecIm.getArray(), '<=', 0., 1., normalizedFlat)
        normalizedFlat = drpStella.where(snrArr, '<', 100., 1., normalizedFlat)
        normalizedFlat = drpStella.where(sumFlats, '<=', 0., 1., normalizedFlat)
        normalizedFlat = drpStella.where(sumVariances, '<=', 0., 1., normalizedFlat)
        
        normFlatOut = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(normalizedFlat)))

        self.recordCalibInputs(cache.butler, normFlatOut, struct.ccdIdList, struct.outputId)

        self.interpolateNans(normFlatOut)
        
        normFlatOutDec = afwImage.DecoratedImageF(normFlatOut.getMaskedImage().getImage())

        self.write(cache.butler, normFlatOutDec, struct.outputId)
