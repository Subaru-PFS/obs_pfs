#!/usr/bin/env python
import lsst.afw.detection as afwDet
import lsst.afw.display as afwDisplay
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
import pfs.drp.stella.utils as dsUtils
from pfs.drp.stella.createFlatFiberTraceProfileTask import CreateFlatFiberTraceProfileTask
from pfs.drp.stella.findAndTraceAperturesTask import FindAndTraceAperturesTask
from pfs.drp.stella.datamodelIO import PfsFiberTraceIO

class ConstructFiberTraceConfig(CalibConfig):
    """Configuration for FiberTrace construction"""
    crGrow = Field(
        dtype=int,
        default=2,
        doc="Grow radius for CR (pixels)"
    )
    doRepair = Field(
        dtype=bool,
        default=True,
        doc="Repair artifacts?"
    )
    profile = ConfigurableField(
        target=CreateFlatFiberTraceProfileTask,
        doc="Task to calculate the spatial profile"
    )
    psfFwhm = Field(
        dtype=float,
        default=3.0,
        doc="Repair PSF FWHM (pixels)"
    )
    psfSize = Field(
        dtype=int,
        default=21,
        doc="Repair PSF size (pixels)"
    )
    repair = ConfigurableField(
        target=RepairTask,
        doc="Task to repair artifacts"
    )
    trace = ConfigurableField(
        target=FindAndTraceAperturesTask,
        doc="Task to trace apertures"
    )
    xOffsetHdrKeyWord = Field(
        dtype=str,
        default='sim.slit.xoffset',
        doc="Header keyword for fiber offset in input files"
    )

    def setDefaults(self):
        CalibConfig.setDefaults(self)

class ConstructFiberTraceTask(CalibTask):
    """Task to construct the normalized Flat"""
    ConfigClass = ConstructFiberTraceConfig
    _DefaultName = "constructFiberTrace"
    calibName = "fibertrace"

    def __init__(self, *args, **kwargs):
        CalibTask.__init__(self, *args, **kwargs)
        self.makeSubtask("repair")
        self.makeSubtask("profile")
        self.makeSubtask("trace")

        import lsstDebug
        self.debugInfo = lsstDebug.Info(__name__)

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
                fpSet = afwDet.FootprintSet(mask, afwDet.Threshold(0.5))
                fpSet = afwDet.FootprintSet(fpSet, self.config.crGrow, True)
                fpSet.setMask(exposure.getMaskedImage().getMask(), "CR")

        if self.debugInfo.display and self.debugInfo.display_inputs:
            disp = afwDisplay.Display(frame=self.debugInfo.inputs_frame)

            visit = sensorRef.dataId['visit']
            disp.mtv(exposure, "raw %d" % (visit))
                
        return exposure

    def run(self, expRefList, butler, calibId):
        #Only run for Flats with xOffset == 0.0
        newExpRefList = []
        for expRef in expRefList:
            exposure = expRef.get('raw')
            md = exposure.getMetadata()

            if self.config.xOffsetHdrKeyWord not in md.names():
                self.log.warn("Keyword %s not found in metadata; ignoring flat for %s" %
                              (self.config.xOffsetHdrKeyWord, expRef.dataId))
            else:
                self.log.info('offset = %f' % (md.get(self.config.xOffsetHdrKeyWord)))
                if md.get(self.config.xOffsetHdrKeyWord) == 0.:
                    newExpRefList.append(expRef)

        if len(newExpRefList) == 0:
            self.log.fatal("No flats were found with valid xOffset keyword %s" %
                           self.config.xOffsetHdrKeyWord)
            raise RuntimeError("Unable to find any valid flats")

        return CalibTask.run(self, newExpRefList, butler, calibId)


    def combine(self, cache, struct, outputId):
        """!Combine multiple exposures of a particular CCD and write the output

        Only the slave nodes execute this method.

        @param cache  Process pool cache
        @param struct  Parameters for the combination, which has the following components:
            * ccdName     Name tuple for CCD
            * ccdIdList   List of data identifiers for combination
            * scales      Scales to apply (expScales are scalings for each exposure,
                               ccdScale is final scale for combined image)
        @param outputId    Data identifier for combined image (exposure part only)
        """
        # Check if we need to look up any keys that aren't in the output dataId
        fullOutputId = {k: struct.ccdName[i] for i, k in enumerate(self.config.ccdKeys)}
        self.addMissingKeys(fullOutputId, cache.butler)
        fullOutputId.update(outputId)  # must be after the call to queryMetadata
        outputId = fullOutputId
        del fullOutputId

        dataRefList = [getDataRef(cache.butler, dataId) if dataId is not None else None for
                       dataId in struct.ccdIdList]

        self.log.info("Combining %s on %s" % (outputId, NODE))
        calib = self.combination.run(dataRefList, expScales=struct.scales.expScales,
                                     finalScale=struct.scales.ccdScale)
        calExp = afwImage.makeExposure(calib)

        self.recordCalibInputs(cache.butler, calExp, struct.ccdIdList, outputId)

        self.interpolateNans(calExp)

        if self.debugInfo.display:
            disp = afwDisplay.Display(frame=self.debugInfo.combined_frame)

            disp.mtv(calExp, "Combined")

        fts = self.trace.run(calExp)
        self.log.info('%d FiberTraces found on combined flat' % (fts.size()))

        self.profile.run(fts)

        if self.debugInfo.display:
            disp = afwDisplay.Display(frame=self.debugInfo.combined_frame)

            dsUtils.addFiberTraceSetToMask(calExp.getMaskedImage().getMask(), fts.getTraces(), disp)
            disp.setMaskTransparency(50)
            disp.mtv(calExp, "Traced")
        #
        # Package up fts (a FiberTraceSet) into a pfsFiberTrace for I/O;  these two classes
        # should be merged at some point
        #
        pfsFT = PfsFiberTrace(
            outputId['calibDate'],
            outputId['spectrograph'],
            outputId['arm']
        )
        pfsFT.fwhm = self.trace.config.apertureFWHM
        pfsFT.threshold = self.trace.config.signalThreshold
        pfsFT.nTerms = self.trace.config.nTermsGaussFit
        pfsFT.saturationLevel = self.trace.config.saturationLevel
        pfsFT.minLength = self.trace.config.minLength
        pfsFT.maxLength = self.trace.config.maxLength
        pfsFT.nLost = self.trace.config.nLost
        pfsFT.traceFunction = self.trace.config.interpolation
        pfsFT.order = self.trace.config.order
        pfsFT.xLow = self.trace.config.xLow
        pfsFT.xHigh = self.trace.config.xHigh
        pfsFT.nCutLeft = self.trace.config.nPixCutLeft
        pfsFT.nCutRight = self.trace.config.nPixCutRight
        pfsFT.interpol = self.profile.config.profileInterpolation
        pfsFT.swathLength = self.profile.config.swathWidth
        pfsFT.overSample = self.profile.config.overSample
        pfsFT.maxIterSF = self.profile.config.maxIterSF
        pfsFT.maxIterSig = self.profile.config.maxIterSig
        pfsFT.lambdaSF = self.profile.config.lambdaSF
        pfsFT.lambdaSP = self.profile.config.lambdaSP
        pfsFT.lambdaWing = self.profile.config.wingSmoothFactor
        pfsFT.lSigma = self.profile.config.lowerSigma
        pfsFT.uSigma = self.profile.config.upperSigma

        pfsFT.fiberId = []
        pfsFT.xCenter = []
        pfsFT.yCenter = []
        pfsFT.yLow = []
        pfsFT.yHigh = []
        pfsFT.coeffs = []
        pfsFT.profiles = None           # we'll set it when we know its dimension

        nRows = calib.getImage().getHeight()
        #
        # The code doesn't impose any requirements on the widths of the traces, while the
        # datamodel requires that they are all the same.
        #
        # For now, trim off extra columns
        #
        width = min([fts.getFiberTrace(i).getProfile().getWidth() for i in range(fts.size())])
        pfsFT.profiles = np.empty((fts.size(), nRows, width), dtype=np.float32)
        pfsFT.profiles[:] = np.nan
        
        for iFt in range(fts.size()):
            ft = fts.getFiberTrace(iFt)
            ftFunc = ft.getFiberTraceFunction()
            pfsFT.fiberId.append(ft.getITrace()+1)
            pfsFT.xCenter.append(ftFunc.xCenter)
            pfsFT.yCenter.append(ftFunc.yCenter)
            pfsFT.yLow.append(ftFunc.yLow)
            pfsFT.yHigh.append(ftFunc.yHigh)
            pfsFT.coeffs.append(ftFunc.coefficients)
            yStart = ft.getFiberTraceFunction().yCenter + ft.getFiberTraceFunction().yLow
            yEnd = ft.getFiberTraceFunction().yCenter + ft.getFiberTraceFunction().yHigh + 1

            pfsFT.profiles[iFt, yStart:yEnd, :] = ft.getProfile().getArray()[:,:width]

        pfsFiberTrace = PfsFiberTraceIO(pfsFT)
        self.write(cache.butler, pfsFiberTrace, outputId)
