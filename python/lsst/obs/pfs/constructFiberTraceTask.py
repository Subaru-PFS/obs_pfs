#!/usr/bin/env python
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
from lsst.ctrl.pool.pool import NODE
import lsst.meas.algorithms as measAlg
from lsst.pex.config import Field, ConfigurableField
from lsst.pipe.drivers.constructCalibs import CalibConfig, CalibTask
from lsst.pipe.tasks.repair import RepairTask
from lsst.pipe.drivers.utils import getDataRef
import math
import numpy as np
from pfs.datamodel.pfsFiberTrace import PfsFiberTrace
import pfs.drp.stella.createFlatFiberTraceProfileTask as cfftpTask
import pfs.drp.stella.findAndTraceAperturesTask as fataTask
from pfs.drp.stella.datamodelIO import PfsFiberTraceIO

class ConstructFiberTraceConfig(CalibConfig):
    """Configuration for FiberTrace construction"""
    xOffsetHdrKeyWord = Field(dtype=str, default='sim.slit.xoffset', doc="Header keyword for xOffset in input files")
    doRepair = Field(dtype=bool, default=True, doc="Repair artifacts?")
    psfFwhm = Field(dtype=float, default=3.0, doc="Repair PSF FWHM (pixels)")
    psfSize = Field(dtype=int, default=21, doc="Repair PSF size (pixels)")
    crGrow = Field(dtype=int, default=2, doc="Grow radius for CR (pixels)")
    repair = ConfigurableField(target=RepairTask, doc="Task to repair artifacts")
    darkTime = Field(dtype=str, default="DARKTIME", doc="Header keyword for time since last CCD wipe, or None",
                     optional=True)
    interpolation = Field(
          doc = "Interpolation schemes (CHEBYSHEV, LEGENDRE, CUBIC, LINEAR, POLYNOMIAL[only one implemented atm])",
          dtype = str,
          default = "POLYNOMIAL")
    order = Field(
        doc = "Polynomial order",
        dtype = int,
        default = 5,
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
        doc = "1 to look for maximum only without GaussFit; 3 to fit Gaussian; 4 to fit Gaussian plus constant background, 5 to fit Gaussian plus linear term (sloped backgfound)",
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

class ConstructFiberTraceTask(CalibTask):
    """Task to construct the normalized Flat"""
    ConfigClass = ConstructFiberTraceConfig
    _DefaultName = "constructFiberTrace"
    calibName = "fiberTrace"
#    filterName = "NONE"  # Sets this filter name in the output

    def __init__(self, *args, **kwargs):
        CalibTask.__init__(self, *args, **kwargs)
        self.makeSubtask("repair")

    @classmethod
    def applyOverrides(cls, config):
        """Overrides to apply for FiberTrace construction"""
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

    def run(self, expRefList, butler, calibId):
        """Only run for Flats with xOffset == 0.0"""
        newExpRefList = []
        for expRef in expRefList:
            exposure = expRef.get('raw')
            try:
                if exposure.getMetadata().get(self.config.xOffsetHdrKeyWord) == 0.:
                    newExpRefList.append(expRef)
            except:
                newExpRefList.append(expRef)
            
        return CalibTask.run(self, newExpRefList, butler, calibId)
    
    
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
        calib = self.combination.run(dataRefList, expScales=struct.scales.expScales,
                                     finalScale=struct.scales.ccdScale)

        self.recordCalibInputs(cache.butler, calib, struct.ccdIdList, struct.outputId)

        self.interpolateNans(calib)
        
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
        
        calExp = afwImage.makeExposure(afwImage.makeMaskedImage(calib.getImage()))
        fts = myFindTask.run(calExp)
        
        myProfileTask = cfftpTask.CreateFlatFiberTraceProfileTask()
        myProfileTask.config.profileInterpolation = self.config.profileInterpolation
        myProfileTask.config.maxIterSF = self.config.maxIterSF
        myProfileTask.config.overSample = self.config.overSample
        myProfileTask.config.swathWidth = self.config.swathWidth
        myProfileTask.config.lambdaSF = self.config.lambdaSF
        myProfileTask.config.lambdaSP = self.config.lambdaSP
        myProfileTask.config.wingSmoothFactor = self.config.wingSmoothFactor

        fts = myProfileTask.run(fts)
        
        dataId = struct.ccdIdList[0]
        pfsFT = PfsFiberTrace(struct.outputId['calibDate'], struct.outputId['spectrograph'], struct.outputId['arm'])
#        pfsFT = PfsFiberTrace(struct.outputId['visit'], struct.outputId['spectrograph'], struct.outputId['arm'])
        pfsFT.fwhm = self.config.apertureFWHM
        pfsFT.threshold = self.config.signalThreshold
        pfsFT.nTerms = self.config.nTermsGaussFit
        pfsFT.saturationLevel = self.config.saturationLevel
        pfsFT.minLength = self.config.minLength
        pfsFT.maxLength = self.config.maxLength
        pfsFT.nLost = self.config.nLost
        pfsFT.traceFunction = self.config.interpolation
        pfsFT.order = self.config.order
        pfsFT.xLow = self.config.xLow
        pfsFT.xHigh = self.config.xHigh
        pfsFT.nCutLeft = fts.getFiberTrace(0).getFiberTraceFunction().fiberTraceFunctionControl.nPixCutLeft
        pfsFT.nCutRight = fts.getFiberTrace(0).getFiberTraceFunction().fiberTraceFunctionControl.nPixCutRight
        pfsFT.interpol = self.config.profileInterpolation
        pfsFT.swathLength = self.config.swathWidth
        pfsFT.overSample = self.config.overSample
        pfsFT.maxIterSF = self.config.maxIterSF
        pfsFT.maxIterSig = fts.getFiberTrace(0).getFiberTraceProfileFittingControl().maxIterSig
        pfsFT.lambdaSF = self.config.lambdaSF
        pfsFT.lambdaSP = self.config.lambdaSP
        pfsFT.lambdaWing = self.config.wingSmoothFactor
        pfsFT.lSigma = fts.getFiberTrace(0).getFiberTraceProfileFittingControl().lowerSigma
        pfsFT.uSigma = fts.getFiberTrace(0).getFiberTraceProfileFittingControl().upperSigma
        
        pfsFT.fiberId = []
        pfsFT.xCenter = []
        pfsFT.yCenter = []
        pfsFT.yLow = []
        pfsFT.yHigh = []
        pfsFT.coeffs = []
        pfsFT.profiles = []
        
        nRows = calib.getImage().getHeight()

        for iFt in range(fts.size()):
            ft = fts.getFiberTrace(iFt)
            pfsFT.fiberId.append(ft.getITrace()+1)
            pfsFT.xCenter.append(ft.getFiberTraceFunction().xCenter)
            pfsFT.yCenter.append(ft.getFiberTraceFunction().yCenter)
            pfsFT.yLow.append(ft.getFiberTraceFunction().yLow)
            pfsFT.yHigh.append(ft.getFiberTraceFunction().yHigh)
            pfsFT.coeffs.append(ft.getFiberTraceFunction().coefficients)
            prof = ft.getProfile()
            profOut = np.ndarray(shape=[nRows,prof.getWidth()], dtype=np.float32)
            profOut[:,:] = 0.
            profOut[ft.getFiberTraceFunction().yCenter + ft.getFiberTraceFunction().yLow:ft.getFiberTraceFunction().yCenter + ft.getFiberTraceFunction().yHigh+1,:] = prof.getArray()[:,:]
            pfsFT.profiles.append(profOut)

        self.write(cache.butler, PfsFiberTraceIO(pfsFT), struct.outputId)
