import os

import matplotlib.pyplot as plt
import numpy as np
try:
    import pyfits
except ImportError:
    pyfits = None

import lsst.afw.image as afwImage
import lsst.log as log
from lsst.obs.subaru.isr import SubaruIsrTask
from lsst.pex.config import Config, ConfigurableField, Field, ListField
from lsst.pipe.base import Struct, TaskRunner, ArgumentParser, CmdLineTask
from lsst.pipe.drivers.constructCalibs import CalibCombineTask
from lsst.utils import getPackageDir
import pfs.drp.stella as drpStella
from pfs.drp.stella.extractSpectraTask import ExtractSpectraTask
from pfs.drp.stella.utils import createLineListForFiberTrace
from pfs.drp.stella.utils import createLineListForLamps
from pfs.drp.stella.utils import getElements, getElementsString
from pfs.drp.stella.utils import makeFiberTraceSet, measureLinesInWavelengthSpace
from pfs.drp.stella.utils import plotWavelengthResiduals
from pfs.drp.stella.utils import readWavelengthFile
from pfs.drp.stella.utils import writePfsArm

# For meanings of flags please refer to drp_stella/include/pfs/drp/stella/Lines.h

fileNameLampsOut = '%s_%s.fits'# % (elements, lineListSuffix)

class ConstructLineListConfig(Config):
    """Configuration for reducing arc images"""
    aperturesToExtract = ListField(
        doc = "Indix of FiberTrace to use (-1 to use all FiberTraces)",
        dtype = int,
        default = [-1]
    )
    clobber = Field(
        doc = "Overwrite existing line lists?",
        dtype = bool,
        default = False
    )
    combine = ConfigurableField(
        target=CalibCombineTask,
        doc="Task for the combination of multiple images"
    )
    extract = ConfigurableField(
        target=ExtractSpectraTask,
        doc="Task for the Instrumental Signature Removal"
    )
    isr = ConfigurableField(
        target=SubaruIsrTask,
        doc="Task for the Instrumental Signature Removal"
    )
    lamps = ListField(
        doc = "Arc lamp(s) to create line list for, separated by ',' (e.g. Hg,Ar)",
        dtype = str,
        default = ["Ne", "Xe", "Hg,Ar"]
    )
    lampKeyword = Field(
        doc = "Fits header keywords for the lamps",
        dtype = str,
        default = "W_AIT_SRC_%s"
    )
    lineListSuffix = Field(
        doc = "Suffix of line list to read (vac or air)",
        dtype = str,
        default = 'vac'
    )
    function = Field(
        doc = "Function for fitting the dispersion",
        dtype = str,
        default = "POLYNOMIAL"
    )
    fwhm = Field(
        doc = "FWHM of emission lines in pixels",
        dtype = float,
        default = 2.2
    )
    maxDistance = Field(
        doc = "Reject emission lines which center is more than this value in pixels away from the predicted position",
        dtype = float,
        default = 1.5
    )
    maxDLambda = Field(
        doc = "Reject emission lines which center is more than this value in nm away from the predicted position",
        dtype = float,
        default = 0.02
    )
    minDistance = Field(
        doc="Minimum distance between lines for creation of line list in FWHM (see below)",
        dtype = float,
        default = 2.2
    )
    minDistanceLines = Field(
        doc="Minimum distance in pixels between 2 lines to be identified",
        dtype = float,
        default = 2.5
    )
    minErr = Field(
        doc="Minimum measure error for PolyFit",
        dtype = float,
        default = 0.00001
    )
    minPercentageOfLines = Field(
        doc = "Minimum percentage of lines to be identified for <identify> to pass",
        dtype = float,
        default = 66.6
    )
    minRatio = Field(
        doc = "Minimum strength ratio between adjacent lines",
        dtype = float,
        default = 100.
    )
    minStrength = Field(
        doc = "Minimum strength of lines",
        dtype = int,
        default = 50
    )
    nIterReject = Field(
        doc = "Number of sigma rejection iterations",
        dtype = int,
        default = 3
    )
    nRowsPrescan = Field(
        doc = "Number of prescan rows in raw CCD image",
        dtype = int,
        default = 48
    )
    order = Field(
        doc = "Fitting function order",
        dtype = int,
        default = 6
    )
    plot = ListField(
        doc = "FiberTraceId to plot (set to -1 for none)",
        dtype = int,
        default = [-1]
    )
    removeLines = Field(
        doc = "Remove lines from line list (e.g. HgII,NeII)",
        dtype = str,
        default = 'HgII,NeII'
    )
    searchRadius = Field(
        doc = "Radius in pixels relative to line list to search for emission line peak",
        dtype = int,
        default = 2
    )
    sigmaReject = Field(
        doc = "Sigma rejection threshold for polynomial fitting",
        dtype = float,
        default = 2.5,
    )
    wavelengthFile = Field(
        doc = "Reference pixel-wavelength file including path",
        dtype = str,
        default = os.path.join(getPackageDir("obs_pfs"), "pfs/RedFiberPixels.fits.gz")
    )
    xCorRadius = Field(
        doc = "Radius in pixels for cross correlating spectrum and line list",
        dtype = int,
        default = 15,
    )

class ConstructLineListTaskRunner(TaskRunner):
    """Get parsed values into the ConstructLineListTask.run"""
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [dict(expRefList=parsedCmd.id.refList,
                     butler=parsedCmd.butler,
                     wLenFile=parsedCmd.wLenFile,
                     )]

    def __call__(self, args):
        task = self.TaskClass(config=self.config, log=self.log)
        self.doRaise = True
        if self.doRaise:
            self.log.debug('ConstructLineListTask.__call__: args = %s' % args)
            result = task.run(**args)
        else:
            try:
                result = task.run(**args)
            except Exception, e:
                task.log.warn("Failed: %s" % e)

        if self.doReturnResults:
            return Struct(
                args = args,
                metadata = task.metadata,
                result = result,
            )

class ConstructLineListTask(CmdLineTask):
    """Task to construct the line list"""
    ConfigClass = ConstructLineListConfig
    RunnerClass = ConstructLineListTaskRunner
    _DefaultName = "constructLineList"

    def __init__(self, *args, **kwargs):
        CmdLineTask.__init__(self, *args, **kwargs)
        self.makeSubtask("combine")
        self.makeSubtask("extract")
        self.makeSubtask("isr")

        # create DispCorControl
        self.dispCorControl = drpStella.DispCorControl()
        self.dispCorControl.fittingFunction = self.config.function
        self.dispCorControl.fwhm = self.config.fwhm
        self.dispCorControl.maxDistance = self.config.maxDistance
        self.dispCorControl.minDistanceLines = self.config.minDistanceLines
        self.dispCorControl.minErr = self.config.minErr
        self.dispCorControl.minPercentageOfLines = self.config.minPercentageOfLines
        self.dispCorControl.nIterReject = self.config.nIterReject
        self.dispCorControl.order = self.config.order
        self.dispCorControl.percentageOfLinesForCheck = 0
        self.dispCorControl.searchRadius = self.config.searchRadius
        self.dispCorControl.sigmaReject = self.config.sigmaReject
        self.dispCorControl.verticalPrescanHeight = self.config.nRowsPrescan


    @classmethod
    def applyOverrides(cls, config):
        """Overrides to apply for line list construction"""
        config.isr.doLinearize = False
        config.isr.doGuider = False
        config.isr.doDefect = False
        config.isr.doFringe = False
        config.isr.doCrosstalk = False
        config.isr.doWrite = True

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", datasetType="raw",
                               help="input identifiers, e.g., --id visit=123 ccd=4")
        parser.add_argument("--wLenFile", help='directory and name of pixel vs. wavelength file')
        return parser

    def run(self,
            expRefList,
            butler,
            wLenFile=None,
            immediate=True,
           ):
        """
        The method for the line list creation is as follows:
        * combine the NIST line lists for the individual elements of the Arc lamps
        * extract the Arc FiberTraces
        * create an artificial spectrum from the NIST wavelengths and predicted
          strengths
        * cross-correlate the Arc spectrum with the artificial spectrum to find
          the offset
        * fit each line with a Gaussian to find the position and strength
        * remove lines which have a measured strength of less than minStrength,
          blends, and lines which could not be fit or are flagged as problematic

        Blends are identified as follows:
        For each line search for lines with a fitted position close to its own
        position (within minDistance * fwhm pixels)
        If a possible blend is identified, compare the predicted lines strengths.
        If one line is at least minRatio times stronger than the other lines in
        the blend, keep the strong line, otherwise discard all lines in the blend
        @param expRefList : reference list of Arc exposures
        @param butler : butler to use
        @param wLenFile : simulator output with the predicted wavelengths for
                          each pixel
        @param immediate : let butler read file immediately or only when needed?
        """
        # Silence verbose loggers
        log.setLevel("afw.ExposureFormatter", log.WARN)
        log.setLevel("afw.image.ExposureInfo", log.WARN)
        log.setLevel("afw.image.Mask", log.WARN)
        log.setLevel("afw.ImageFormatter", log.WARN)
        log.setLevel("CameraMapper", log.WARN)
        log.setLevel("createLineListForFiberTrace", log.WARN)
        log.setLevel("createLineListForLamps", log.WARN)
        log.setLevel("daf.persistence.butler", log.WARN)
        log.setLevel("daf.persistence.LogicalLocation", log.WARN)
        log.setLevel("extractSpectra", log.WARN)
        log.setLevel("gaussFit", log.WARN)
        log.setLevel("gaussFunc", log.FATAL)
        log.setLevel("getLinesInWavelengthRange", log.WARN)
        log.setLevel("makeArtificialSpectrum", log.WARN)
        log.setLevel("measureLines", log.WARN)
        log.setLevel("measureLinesInPixelSpace", log.WARN)
        log.setLevel("removeBadLines", log.WARN)
        log.setLevel("pfs.drp.stella.math.assignITrace", log.WARN)
        log.setLevel("pfs.drp.stella.math.CurfFitting.PolyFit", log.WARN)
        log.setLevel("pfs.drp.stella.math.findITrace", log.WARN)
        log.setLevel("pfs.drp.stella.Spectra.identify", log.TRACE)

        if wLenFile == None:
            wLenFile = self.config.wavelengthFile
        self.log.debug('len(expRefList) = %d' % len(expRefList))
        self.log.debug('wLenFile = %s' % wLenFile)

        # read wavelength file
        xCenters, lambdaPix, traceIds = readWavelengthFile(wLenFile)
        self.log.debug('unique traceIds = %s' % (np.array_str(np.unique(traceIds))))

        self.log.trace('dispCorControl.fittingFunction = %s' % self.dispCorControl.fittingFunction)
        self.log.trace('dispCorControl.order = %d' % self.dispCorControl.order)
        self.log.trace('dispCorControl.searchRadius = %d' % self.dispCorControl.searchRadius)
        self.log.trace('dispCorControl.fwhm = %g' % self.dispCorControl.fwhm)
        self.log.trace('dispCorControl.maxDistance = %g' % self.dispCorControl.maxDistance)

        self.combine.calibName = 'arc'

        oneLampOnly = []
        lineListPerLamp = []
        spectrumSetPerLamp = []

        # get lamps and lampKeywords
        lamps = []
        lampKeywords = []
        for elem in self.config.lamps:
            lamp = getElementsString(elem)
            lamps.append(lamp)
            lampKeywords.append(self.config.lampKeyword % (lamp))
            oneLampOnly.append([])
        self.log.trace('lamps = %s' % (lamps))
        self.log.trace('lampKeywords = %s' % (lampKeywords))

        # find possible combinations of lamps
        doubleCombos = self.getDoubleCombinations(lamps)
        nCombs = len(doubleCombos)
        self.log.trace('number of combinations with 2 lamps on = %d' % (nCombs))
        self.log.trace('doubleCombos = %s' % (doubleCombos))
        doubleComboList = dict()
        for combo in doubleCombos:
            doubleComboList[combo] = []
        self.log.trace('doubleComboList = %s' % (doubleComboList))

        # for each input exposure reference find Arc lamps which were switched on
        for expRef in expRefList:
            try:
                md = expRef.get('raw_md')
            except Exception as e:
                continue

            lampsOn = []
            for key in lampKeywords:
                try:
                    lampsOn.append(md.get(key))
                except Exception as e:
                    self.log.info('%s' % (e.message))
            self.log.trace('lampsOn = %s' % (lampsOn))

            lampsOnInd = self.getLamps(lampsOn)
            nLamps = len(lampsOnInd)
            if nLamps == 1:
                oneLampOnly[lampsOnInd[0]].append(expRef)
            elif nLamps == 2:
                elemStr = ''
                for i in range(len(lampsOn)):
                    if lampsOn[i]:
                        elemStr += lamps[i] + ','
                elemStr = elemStr[0:len(elemStr)-1]
                doubleComboList[elemStr].append(expRef)

        # One Arc lamp on only:
        # Run ISR Task and combine multiple images, extract spectra, create line list
        for iLamp in range(len(oneLampOnly)):
            lamp = oneLampOnly[iLamp]
            element = self.config.lamps[iLamp]
            self.log.trace('element = %s' % (element))
            self.log.trace('len(lamp) = %d' % (len(lamp)))
            goodLines = drpStella.getNistLineMeasVec()
            spectrumSet = drpStella.SpectrumSet()
            if len(lamp) > 0:
                arcRef = lamp[0]
                for expRef in lamp:
                    self.isr.runDataRef(expRef)
                    if np.isnan(expRef.get('postISRCCD').getMaskedImage().getImage().getArray()).any():
                        raise RuntimeError('after ISR task: NaNs detected')
                if len(lamp) > 1:
                    subset = []
                    for expRef in lamp:
                        subset.append(expRef)#.get('postISRCCD'))
                    combinedImage = self.combine.run(subset)
                    combinedExposure = afwImage.ExposureF(afwImage.MaskedImageF(combinedImage.getImage()))
                else:
                    combinedExposure = expRef.get('postISRCCD')
                if np.isnan(combinedExposure.getMaskedImage().getImage().getArray()).any():
                    raise RuntimeError('combinedExposure: NaNs detected')

                spectrumSet, goodLines = self.createLineList(
                    butler = butler,
                    arcRef = arcRef,
                    arcExp = combinedExposure,
                    elements = element,
                    lambdaPix = lambdaPix,
                    traceIds = traceIds,
                    xCenters = xCenters,
                    immediate = immediate
                )
                self.log.trace('len(goodLines) = %d' % (len(goodLines)))
            lineListPerLamp.append(goodLines)
            spectrumSetPerLamp.append(spectrumSet)

        # Find blends for each possible combination of Arc lamps
        for iCombo in range(len(doubleCombos)):
            combo = doubleCombos[iCombo]
            self.log.trace('combo = %s' % (combo))
            elements = getElements(combo)
            spectra = []
            lineLists = drpStella.getNistLineMeasVec()
            for element in elements:
                for iLamp in range(len(self.config.lamps)):
                    if element in getElementsString(self.config.lamps[iLamp]):
                        self.log.trace('element %s found in self.config.lamps[%d] = %s'
                                       % (element, iLamp, self.config.lamps[iLamp]))
                        if spectrumSetPerLamp[iLamp].size() > 0:
                            spectra.append(spectrumSetPerLamp[iLamp].getSpectrum(
                                           int(spectrumSetPerLamp[iLamp].size() / 2)))
                            self.log.trace('element <%s> found in self.config.lamps[%d] = %s'
                                           % (element, iLamp, self.config.lamps[iLamp]))
                            lineListTemp = lineListPerLamp[iLamp]
                            self.log.trace('len(lineLists) = %d, len(lineListTemp) = %d'
                                           % (len(lineLists), len(lineListTemp)))
                            lineLists.extend(lineListTemp)
                            self.log.trace('len(lineLists) = %d' % (len(lineLists)))
            self.log.trace('len(spectra) = %d' % (len(spectra)))
            if len(spectra) > 1:
                # Add spectra in wavelength space
                spec = spectra[0]
                specSum = spec.getSpectrum()
                for iSpec in np.arange(1,len(spectra)):
                    spectrum = spectra[iSpec]
                    specSum += np.interp(spec.getWavelength(),
                                         spectrum.getWavelength(),
                                         spectrum.getSpectrum())
                    if spec.getITrace() in self.config.plot:
                        plt.plot(spec.getWavelength(), spec.getSpectrum(), 'g-', label = 'spectrum %s' % (getElements(combo)[0]))
                        plt.plot(spectrum.getWavelength(), spectrum.getSpectrum(), 'b-', label = 'spectrum %s' % (getElements(combo)[1]))
                        plt.plot(spec.getWavelength(), specSum, 'r-', label = 'sum')
                        plt.legend()
                        plt.show()

                # fit Gaussians in wavelength space
                measureLinesInWavelengthSpace(lineLists,
                                              specSum,
                                              spec.getWavelength(),
                                              sigma=0.2,
                                              plot=(spec.getITrace() in self.config.plot))
                self.flagBlends(lineLists, combo)
                nGoodLines = 0
                for iLine in range(len(lineLists)):
                    if 'l' in lineLists[iLine].flags:
                        self.log.debug('line %d with ID %d marked as blend, lambda=%f'
                                       % (iLine,
                                          lineLists[iLine].nistLine.id,
                                          lineLists[iLine].nistLine.laboratoryWavelength))
                    else:
                        nGoodLines += 1
                self.log.info('nGoodLines = %d' % (nGoodLines))
        finalLineList = drpStella.getNistLineMeasVec()
        for lineList in lineListPerLamp:
            finalLineList.extend(lineList)
        self.log.info('Final line list contains %d lines' % (len(finalLineList)))

        # write output line list
        elements = ""
        for lamp in self.config.lamps:
            elements += getElementsString(lamp)
        outputFileName = os.path.join(getPackageDir("obs_pfs"),
                                      "pfs/lineLists/%s_%d%s.fits"
                                      % (elements,
                                         expRefList[0].dataId['spectrograph'],
                                         expRefList[0].dataId['arm']))
        if os.path.isfile(outputFileName) and not self.config.clobber:
            self.log.warn('existing line list found, set clobber to True to overwrite')
        else:
            self.writeLineList(finalLineList,
                               expRefList[0].dataId['spectrograph'],
                               expRefList[0].dataId['arm'])

    def getLamps(self, lampsOn):
        """
        Check which Arc lamps were on and return a tuple containing the indices of those lamps
        @param lampsOn : tuple of booleans, one per lamp
        @return : tuple of indices where lampsOn is True
        """
        lamps = []
        for iLamp in range(len(lampsOn)):
            if lampsOn[iLamp]:
                lamps.append(iLamp)
        return lamps

    def getDoubleCombinations(self, lamps):
        """
        Return all possible combinations of 2 Arc lamps
        @param lamps : tuple of individual Arc lamps
        @return : tuple containing the possible combinations, elements separated by ','
        """
        combos = []
        for first in range(len(lamps)-1):
            for second in np.arange(first+1, len(lamps)):
                combos.append('%s,%s' % (lamps[first], lamps[second]))
        return combos

    def createLineList(self,
                       butler,
                       arcRef,
                       arcExp,
                       elements,
                       lambdaPix,
                       traceIds,
                       xCenters,
                       immediate=False):
        """
        extract Arc spectra and create a list of good lines for the wavelength calibration
        The method for the line list creation is as follows:
        * combine the NIST line lists for the individual elements of the Arc lamps
        * extract the Arc FiberTraces
        * create an artificial spectrum from the NIST wavelengths and predicted
          strengths
        * cross-correlate the Arc spectrum with the artificial spectrum to find
          the offset
        * fit each line with a Gaussian to find the position and strength
        * remove lines which have a measured strength of less than minStrength,
          blends, and lines which could not be fit or are flagged as problematic

        Blends are identified as follows:
        For each line search for lines with a fitted position close to its own
        position (within minDistance * fwhm pixels)
        If a possible blend is identified, compare the predicted lines strengths.
        If one line is at least minRatio times stronger than the other lines in
        the blend, keep the strong line, otherwise discard all lines in the blend
        @param butler : Butler to use
        @param arcRef : Arc reference exposure
        @param arcExp : Arc exposure, possible combined from multiple exposures
        @param elements : Arc lamp used for the exposure (e.g. 'Hg,Ar', 'Ne')
        @param lambdaPix : predicted wavelength per pixel and fiberTrace from wavelength file
        @param traceIds : Trace IDs for each FiberTrace from wavelength file
        @param xCenters : center positions of each FiberTrace on CCD in pixels
        @param immediate : Load objects from butler immediately?
        @return : [spectrumSetFromProfile, goodLines]
        """
        elementsString = getElementsString(elements)
        self.log.trace('elementsString = %s' % (elementsString))

        lines = createLineListForLamps(elements,
                                       self.config.lineListSuffix,
                                       self.config.removeLines)
        self.log.info('raw line list contains %d lines' % (len(lines)))

        measuredLinesPerArc = []
        offsetPerArc = []

        self.log.debug('arcRef.dataId = %s' % arcRef.dataId)

        # construct fiberTraceSet from pfsFiberTrace
        try:
            fiberTrace = arcRef.get('fibertrace', immediate=immediate)
            flatFiberTraceSet = makeFiberTraceSet(fiberTrace)
        except Exception, e:
            raise RuntimeError("Unable to load fiberTrace for %s from %s: %s" %
                               (arcRef.dataId,
                                arcRef.get('fiberTrace_filename')[0], e))

        for iFT in range(flatFiberTraceSet.size()):
            ft = flatFiberTraceSet.getFiberTrace(iFT)
            self.log.trace(
                'flatFiberTraceSet[%d].getFiberTraceFunction().fiberTraceFunctionControl.xLow = %f'
                % (iFT, ft.getFiberTraceFunction().fiberTraceFunctionControl.xLow))
            self.log.trace(
                'flatFiberTraceSet[%d].getFiberTraceFunction().fiberTraceFunctionControl.xHigh = %f'
                % (iFT, ft.getFiberTraceFunction().fiberTraceFunctionControl.xHigh))
            self.log.trace(
                'flatFiberTraceSet[%d].getFiberTraceFunction().fiberTraceFunctionControl.nPixCutLeft = %f'
                % (iFT, ft.getFiberTraceFunction().fiberTraceFunctionControl.nPixCutLeft))
            self.log.trace(
                'flatFiberTraceSet[%d].getFiberTraceFunction().fiberTraceFunctionControl.nPixCutRight = %f'
                % (iFT, ft.getFiberTraceFunction().fiberTraceFunctionControl.nPixCutRight))
            if iFT == 4:
                self.log.trace('ft.getXCenters() = %s' % np.array_str(ft.getXCenters()))

        self.log.trace('arcExp = %s' % arcExp)
        self.log.trace('type(arcExp) = %s' % type(arcExp))

        self.log.debug('extracting arc spectra')

        # assign trace number to flatFiberTraceSet
        drpStella.assignITrace( flatFiberTraceSet, traceIds, xCenters )
        for i in range( flatFiberTraceSet.size() ):
            self.log.info('iTraces[%d] = %d' % (i, flatFiberTraceSet.getFiberTrace(i).getITrace()))

        # optimally extract arc spectra
        for iFT in range(flatFiberTraceSet.size()):
            ft = flatFiberTraceSet.getFiberTrace(iFT)
            assert(ft.getProfile().getHeight() == ft.getTrace().getHeight())
            assert(ft.getProfile().getWidth() == ft.getTrace().getWidth())
        spectrumSetFromProfile = self.extract.run(arcExp,
                                                  flatFiberTraceSet,
                                                  self.config.aperturesToExtract)
        self.log.info('spectrumSetFromProfile.size() = %d' % (spectrumSetFromProfile.size()))

        measuredLinesPerFiberTrace = []
        offsetPerFiberTrace = []
        idPerFiberTrace = []
        goodLines = drpStella.getNistLineMeasVec()
        if self.config.aperturesToExtract[0] == -1:
            apsToExtract = range(spectrumSetFromProfile.size())
            extractedSpectra = range(spectrumSetFromProfile.size())
        else:
            apsToExtract = self.config.aperturesToExtract
            extractedSpectra = [0]
        for i in range(len(apsToExtract)):
            spec = spectrumSetFromProfile.getSpectrum(extractedSpectra[i])
            self.log.trace('i = %d: spec.getITrace() = %d' % (i, spec.getITrace()))
            fluxPix = spec.getSpectrum()
            self.log.trace('fluxPix.shape = %d' % fluxPix.shape)
            self.log.trace('type(fluxPix) = %s: <%s>' % (type(fluxPix),type(fluxPix[0])))
            self.log.trace('type(spec) = %s: <%s>: <%s>'
                % (type(spec),type(spec.getSpectrum()),type(spec.getSpectrum()[0])))

            traceId = spec.getITrace()
            idPerFiberTrace.append(traceId)
            self.log.debug('traceId = %d' % traceId)

            # cut off both ends of wavelengths where is no signal
            yMin = (flatFiberTraceSet.getFiberTrace(apsToExtract[i]).getFiberTraceFunction().yCenter +
                    flatFiberTraceSet.getFiberTrace(apsToExtract[i]).getFiberTraceFunction().yLow)
            yMax = (flatFiberTraceSet.getFiberTrace(apsToExtract[i]).getFiberTraceFunction().yCenter +
                    flatFiberTraceSet.getFiberTrace(apsToExtract[i]).getFiberTraceFunction().yHigh)
            self.log.debug('fiberTrace %d: yMin = %d' % (apsToExtract[i], yMin))
            self.log.debug('fiberTrace %d: yMax = %d' % (apsToExtract[i], yMax))

            startIndex = drpStella.firstIndexWithValueGEFrom(traceIds,traceId)
            self.log.trace('startIndex = %d' % startIndex)
            self.log.trace('copying lambdaPix[%d:%d]'
                % (startIndex + self.config.nRowsPrescan + yMin,
                   startIndex + self.config.nRowsPrescan + yMax + 1))
            lambdaPixFiberTrace = lambdaPix[startIndex + self.config.nRowsPrescan + yMin:
                                            startIndex + self.config.nRowsPrescan + yMax + 1]
            if len(lambdaPixFiberTrace) != len(fluxPix):
                raise RuntimeError(
                    "constructLineListTask.py: ERROR: len(lambdaPixFiberTrace)(=%d) != len(fluxPix)(=%d)"
                    % (len(lambdaPixFiberTrace), len(fluxPix)))

            createLineListFiberTraceResult = createLineListForFiberTrace(
                lambdaPix = lambdaPixFiberTrace,
                fluxPix = spec.getSpectrum(),
                lines = lines,
                plot = (traceId in self.config.plot),
                fwhm = self.config.fwhm,
                xCorRadius = self.config.xCorRadius,
                minDistance = self.config.minDistance,
                maxDistance = self.config.maxDistance,
                minStrength = self.config.minStrength,
                minRatio = self.config.minRatio
            )
            self.log.trace('constructLineListTask: createLineListForFiberTrace finished')
            lineList = createLineListFiberTraceResult["lineList"]
            self.log.debug('type(lineList) = %s, len(lineList) = %d'
                           % (type(lineList), len(lineList)))
            measuredLines = createLineListFiberTraceResult["measuredLines"]
            self.log.debug('type(measuredLines) = %s, len(measuredLines) = %d'
                           % (type(measuredLines), len(measuredLines)))
            offset = createLineListFiberTraceResult["offset"]
            self.log.trace('offset = %s' % (np.array_str(offset)))

            self.log.debug('new line list created')
            self.log.trace('len(lineList) = %d' % (len(lineList)))

            # use line list to calibrate the Arc spectrum
            try:
                spec.identify(lineList, self.dispCorControl)
            except Exception, e:
                raise RuntimeError(
                    "constructLineListTask.py: %dth FiberTrace: traceId = %d: ERROR: %s"
                    % (i,traceId,e.message))
            spectrumSetFromProfile.setSpectrum(i, spec)
            if spec.getITrace() in self.config.plot:
                plotWavelengthResiduals(lineList, spec.getITrace(), spec.getDispRms())

            for line in lineList:
                if 'g' in line.flags and 'i' not in line.flags and 'f' not in line.flags:
                    if len(drpStella.getLinesWithID(goodLines, line.nistLine.id)) == 0:
                        goodLines.append(line)

            measuredLinesPerFiberTrace.append(measuredLines)
            offsetPerFiberTrace.append(offset)

            self.log.trace("FiberTrace %d: spec.getWavelength() = %s"
                % (i, np.array_str(spec.getWavelength())))
            self.log.trace("FiberTrace %d: spec.getDispCoeffs() = %s"
                % (i,np.array_str(spec.getDispCoeffs())))
            self.log.info("FiberTrace %d (ID=%d): spec.getDispRms() = %f"
                % (i, spec.getITrace(), spec.getDispRms()))
            self.log.info("FiberTrace %d (ID=%d): spec.getDispRmsCheck() = %f"
                % (i, spec.getITrace(), spec.getDispRmsCheck()))
            spectrumSetFromProfile.setSpectrum(i, spec)

        measuredLinesPerArc.append(measuredLinesPerFiberTrace)
        offsetPerArc.append(offsetPerFiberTrace)

        return [spectrumSetFromProfile, goodLines]

    def flagBlends(self, lineList, lamps):
        """
        Flag blended lines in combined spectra
        @param lineList : input line list containing NistLineMeas
        @param lamps : tuple containing the 2 lamps, e.g. ['HgAr', 'Ne']
        @return : None
        """
        elements = getElements(lamps)
        for line in lineList:
            if (np.fabs(line.gaussCoeffsLambda.mu - line.nistLine.laboratoryWavelength)
                > self.config.maxDLambda):
                if 'l' not in line.flags:
                    line.flags += 'l'
                if line.nistLine.element in elements[0]:
                    line.flags += elements[1]
                else:
                    line.flags += elements[0]

    def writeLineList(self, lineList, spectrograph, arm):
        """
        Write line list to output file. Note that the output file will only contain the NistLines
        from NistLineMeas in lineList
        @param lineList : list of NistLineMeas to write to output file
        """
        if not pyfits:
            raise RuntimeError("I failed to import pyfits, so cannot read from disk")

        elements = ""
        for lamp in self.config.lamps:
            elements += getElementsString(lamp)
        outputFileName = os.path.join(getPackageDir("obs_pfs"),
                                      "pfs/lineLists/%s_%d%s.fits"
                                      % (elements, spectrograph, arm))
        self.log.trace('writing output file <%s>' % outputFileName)

        hdus = pyfits.HDUList()

        hdr = pyfits.Header()
        hdr.update()
        hdus.append(pyfits.PrimaryHDU(header=hdr))

        colDefs = []
        # create data columns
        col1 = pyfits.Column(name='element',
                             format='2A',
                             array=[line.nistLine.element for line in lineList])
        colDefs.append(col1)

        col2 = pyfits.Column(name='ion',
                             format='2A',
                             array=[line.nistLine.ion for line in lineList])
        colDefs.append(col2)

        col3 = pyfits.Column(name='wavelength',
                             format='E',
                             array=[line.nistLine.laboratoryWavelength for line in lineList])
        colDefs.append(col3)

        col4 = pyfits.Column(name='strength',
                             format='E',
                             array=[line.nistLine.predictedStrength for line in lineList])
        colDefs.append(col4)

        col5 = pyfits.Column(name='sources',
                             format='128A',
                             array=[line.nistLine.sources for line in lineList])
        colDefs.append(col5)

        flags = []
        for line in lineList:
            if 'l' in line.flags:
                flags.append(line.flags[line.flags.index('l')+1:])
                self.log.trace(
                    'flag l found in line with ID %d, element=%s: line.flags = <%s> => adding flag <%s>'
                    % (line.nistLine.id,line.nistLine.element, line.flags, flags[len(flags)-1]))
            else:
                flags.append('')
        col6 = pyfits.Column(name='flags',
                             format='16A',
                             array=flags)
        colDefs.append(col6)

        # write fits file
        cols = pyfits.ColDefs(colDefs)
        tbhdu = pyfits.BinTableHDU.from_columns(cols)

        tbhdu.name = ("NistLines")
        hdus.append(tbhdu)
        with open(os.path.join(getPackageDir("obs_pfs"), outputFileName), "w") as fd:
            hdus.writeto(fd)
