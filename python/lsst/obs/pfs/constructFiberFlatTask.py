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
from pfs.drp.stella.createFlatFiberTraceProfileTask import CreateFlatFiberTraceProfileTask
from pfs.drp.stella.findAndTraceAperturesTask import FindAndTraceAperturesTask

class ConstructFiberFlatConfig(CalibConfig):
    """Configuration for flat construction"""
    xOffsetHdrKeyWord = Field(dtype=str,
                              default='sim.slit.xoffset',
                              doc="Header keyword for fiber offset in input files")
    doRepair = Field(dtype=bool,
        default=True,
        doc="Repair artifacts?")
    psfFwhm = Field(dtype=float,
        default=3.0,
        doc="Repair PSF FWHM (pixels)")
    psfSize = Field(dtype=int,
        default=21,
        doc="Repair PSF size (pixels)")
    crGrow = Field(dtype=int,
        default=2,
        doc="Grow radius for CR (pixels)")
    repair = ConfigurableField(target=RepairTask,
        doc="Task to repair artifacts")
    trace = ConfigurableField(target=FindAndTraceAperturesTask,
        doc="Task to trace apertures")
    profile = ConfigurableField(target=CreateFlatFiberTraceProfileTask,
        doc="Task to calculate the spatial profile")
    darkTime = Field(dtype=str,
        default="DARKTIME",
        doc="Header keyword for time since last CCD wipe, or None",
        optional=True)
    display = Field(dtype=bool,
        default=True,
        doc="Display images?")

class ConstructFiberFlatTask(CalibTask):
    """Task to construct the normalized Flat"""
    ConfigClass = ConstructFiberFlatConfig
    _DefaultName = "constructFiberFlat"
    calibName = "flat"

    def __init__(self, *args, **kwargs):
        CalibTask.__init__(self, *args, **kwargs)
        self.makeSubtask("profile")
        self.makeSubtask("repair")
        self.makeSubtask("trace")

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

        exposure = dataRefList[0].get('postISRCCD')
        sumFlats = exposure.getMaskedImage().getImage().getArray()
        sumVariances = exposure.getMaskedImage().getVariance().getArray()

        allFts = []
        xOffsets = []
        for expRef in dataRefList:
            exposure = expRef.get('postISRCCD')
            try:
                xOffsets.append(exposure.getMetadata().get(self.config.xOffsetHdrKeyWord))
            except Exception:
                xOffsets.append(0.0)
            fts = self.trace.run(exposure)
            allFts.append(fts)
            if expRef.dataId['visit'] != dataRefList[0].dataId['visit']:
                sumFlats += exposure.getMaskedImage().getImage().getArray()
                sumVariances += exposure.getMaskedImage().getVariance().getArray()

        self.log.info('=== xOffsets = '+str(xOffsets)+' ===')

        # Calculate spatial profiles for all FiberTraceSets
        for fts in allFts:
            fts = self.profile.run(fts)

        sumRec = np.zeros(shape=sumFlats.shape, dtype='float32')
        sumRecIm = afwImage.ImageF(sumRec)
        sumVar = np.zeros(shape=sumFlats.shape, dtype='float32')
        sumVarIm = afwImage.ImageF(sumVar)

        # Add all reconstructed FiberTraces of all dithered flats to one
        # reconstructed image 'sumRecIm'
        for fts in allFts:
            for ft in fts.getTraces():
                spectrum = ft.extractFromProfile()
                recFt = ft.getReconstructed2DSpectrum(spectrum)
                recFtArr = recFt.getArray()
                imArr = ft.getImage().getArray()
                recFtArr[imArr <= 0] = 0.0
                drpStella.addFiberTraceToCcdImage(ft, recFt, sumRecIm)
                drpStella.addFiberTraceToCcdImage(ft, ft.getVariance(), sumVarIm)

        sumVariances[sumVariances <= 0.0] = 0.1
        snrArr = sumFlats / np.sqrt(sumVariances)

        sumRecImArr = sumRecIm.getArray()
        #to avoid division by zero and remove non-physical negative values,
        #we set the reconstructed values <= 0.0 to 0.01. We later set the normalized
        #Flat to unity where sumRecImArr <= 0.01, as well as all other pixels
        #with a SNR < 100
        sumRecImArr[sumRecImArr <= 0.0] = 0.01
        normalizedFlat = sumFlats / sumRecImArr
        normalizedFlat[sumRecImArr <= 0.01] = 1.0
        normalizedFlat[snrArr < 100.0] = 1.0
        normalizedFlat[sumFlats <= 0.0] = 1.0
        normalizedFlat[sumVariances <= 0.0] = 1.0

        import lsstDebug
        try:
            import debug
            if lsstDebug.Info(__name__).display:
                import lsst.afw.display as afwDisplay
                display = afwDisplay.getDisplay(frame=1)
                display.mtv(afwImage.ImageF(sumFlats - sumRecImArr), title='sumFlats - sumRecIm')
                display.frame = 2
                display.mtv(afwImage.ImageF(normalizedFlat), title='normalized Flat')
        except ImportError:
            debug = None

        """Write fiber flat"""
        normFlatOut = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(normalizedFlat)))
        self.recordCalibInputs(cache.butler, normFlatOut, struct.ccdIdList, struct.outputId)
        self.interpolateNans(normFlatOut)
        normFlatOutDec = afwImage.DecoratedImageF(normFlatOut.getMaskedImage().getImage())
        self.write(cache.butler, normFlatOutDec, struct.outputId)
