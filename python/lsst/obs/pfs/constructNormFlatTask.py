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

class ConstructNormFlatConfig(CalibConfig):
    """Configuration for flat construction"""
    doRepair = Field(dtype=bool, default=True, doc="Repair artifacts?")
    psfFwhm = Field(dtype=float, default=3.0, doc="Repair PSF FWHM (pixels)")
    psfSize = Field(dtype=int, default=21, doc="Repair PSF size (pixels)")
    crGrow = Field(dtype=int, default=2, doc="Grow radius for CR (pixels)")
    repair = ConfigurableField(target=RepairTask, doc="Task to repair artifacts")
    darkTime = Field(dtype=str, default="DARKTIME", doc="Header keyword for time since last CCD wipe, or None",
                     optional=True)

class ConstructNormFlatTask(CalibTask):
    """Task to construct the normalized Flat"""
    ConfigClass = ConstructNormFlatConfig
    _DefaultName = "constructNormFlat"
    calibName = "fiberFlat"

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
        myProfileTask.config.profileInterpolation = 'PISKUNOV'
        myProfileTask.config.ccdReadOutNoise = 3.19
        myProfileTask.config.maxIterSF = 15
        myProfileTask.config.overSample = 15
        myProfileTask.config.swathWidth = 250
        myProfileTask.config.lambdaSF = 1./float(myProfileTask.config.overSample)
        myProfileTask.config.lambdaSP = 0.
        myProfileTask.config.wingSmoothFactor = 0.

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
