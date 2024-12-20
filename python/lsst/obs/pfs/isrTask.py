# This file is part of obs_pfs
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import os
import numpy as np
import numpy.linalg
import lsst.log
import lsst.geom as geom                # noqa F401; used in eval(darkBBoxes)
import lsst.pex.config as pexConfig
from lsst.utils import getPackageDir

import lsst.afw.display as afwDisplay
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
from lsst.afw.image import DecoratedImageF
import lsst.ip.isr as ipIsr
from lsst.ip.isr.assembleCcdTask import AssembleCcdTask
from lsst.ip.isr import isrQa
from lsst.ip.isr.defects import Defects
from lsst.utils.timer import timeMethod
import lsst.pipe.base as pipeBase
from lsst.pipe.base.connectionTypes import Input as InputConnection
from lsst.pipe.base.connectionTypes import PrerequisiteInput as PrerequisiteConnection
from lsst.pipe.base.connectionTypes import Output as OutputConnection
from lsst.pipe.base import Struct, PipelineTaskConnections
from lsst.daf.butler import DimensionGraph

from .utils import getCalibPath
from pfs.drp.stella.crosstalk import PfsCrosstalkTask

___all__ = ["IsrTask", "IsrTaskConfig"]


class BrokenShutterConfig(pexConfig.Config):
    """Configuration parameters for working around a broken shutter
    """
    brokenShutterList = pexConfig.ListField(dtype=str, default=["r1"],
                                            doc="List of cameras with broken shutters (e.g. ['r1'])")

    checkParallelOverscan = pexConfig.Field(dtype=bool, default=False,
                                            doc="Use parallel overscan to guess if the shutter is broken?")
    maximumAllowedParallelOverscanFraction = \
        pexConfig.Field(dtype=float, default=0.01,
                        doc="Largest permitted fraction of flux seen "
                        "in parallel overscan if shutter's working")

    useAnalytic = pexConfig.Field(dtype=bool, default=False,
                                  doc="Use an analytic correction")
    window = pexConfig.Field(dtype=int, default=1,
                             doc="Number of frames +- to search for a bias")

    causalWindow = pexConfig.Field(dtype=bool, default=True,
                                   doc="Only search for later biases?")

    doCacheCorrections = pexConfig.Field(dtype=bool, default=False,
                                         doc="Use a cache of analytic correction matrices?")
    cacheSize = pexConfig.Field(dtype=int, default=3,
                                doc="Number of analytic correction matrices to cache")

    t_wipe = pexConfig.Field(dtype=float, default=5.23,
                             doc="Time taken to wipe the CCD (s)")
    t_read = pexConfig.Field(dtype=float, default=36.79,
                             doc="Time taken to read the CCD (s)")
    t_stare = pexConfig.Field(dtype=float, default=3.06,
                              doc="Time during readout when charge isn't being clocked CCD (s)")

    def validate(self):
        super().validate()
        if self.window > 1 and self.useAnalytic:
            raise ValueError(f"You may not specify both (window == {self.window}) > 1 and useAnalytic")


class PfsAssembleCcdTask(AssembleCcdTask):
    def assembleCcd(self, exposure):
        """Assemble CCD and mask pixels with value of zero

        Pixels with value of zero were intentionally not read out, and only
        exist in the data as padding. We have masked them as NO_DATA
        (PfsIsrTask.maskAmplifier), and here we adjust the value after
        overscan subtraction. The value at the end of ISR will still be slightly
        negative due to bias- and dark-subtraction, but at least it won't be
        negative thousands.

        We multiply each amplifier by the gain, so all data is in units of
        electrons. This is safe because CCD assembly is performed after
        saturation detection (where the units must be ADU), and we will reset
        the gain to unity so any subsequent operation will be satisfied.

        exposure : `lsst.afw.image.Exposure`
            Raw exposure containing all amplifiers and their overscans.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            Assembled exposure.
        """
        exposure = super().assembleCcd(exposure)
        noData = (exposure.mask.array & exposure.mask.getPlaneBitMask("NO_DATA")) != 0
        if np.any(noData):
            exposure.image.array[noData] = 0.0
            self.log.info("Flagging %d zero pixels", noData.sum())

        # Convert to electrons
        for amp in exposure.getDetector():
            ampImage = exposure.maskedImage[amp.getBBox()]
            ampImage *= amp.getGain()
        # Reset gain to unity
        detector = exposure.getDetector().rebuild()
        for amp in detector:
            amp.setGain(1.0)
        exposure.setDetector(detector.finish())

        return exposure


def lookupDefects(datasetType, registry, dataId, collections):
    """Look up defects

    We need to determine the appropriate 'detector' value from arm,spectrograph.

    Parameters
    ----------
    datasetType : `str`
        The dataset type to look up.
    registry : `lsst.daf.butler.Registry`
        The butler registry.
    dataId : `lsst.daf.butler.DataCoordinate`
        The data identifier.
    collections : `list` of `str`
        The collections to search.

    Returns
    -------
    refs : `list` of `lsst.daf.butler.Reference`
        The references to the bias or dark frame.
    """
    results = list(registry.queryDimensionRecords("detector", dataId=dataId))
    if len(results) != 1:
        raise RuntimeError(f"Unable to find detector for {dataId}: {results}")
    detector = results[0].id

    return [registry.findDataset(
        datasetType, collections=collections, dataId=dataId, timespan=dataId.timespan, detector=detector,
    )]


def lookupCrosstalkSources(datasetType, registry, dataId, collections):
    """Look up crosstalk sources

    We don't do inter-CCD crosstalk correction.

    Parameters
    ----------
    datasetType : `str`
        The dataset type to look up.
    registry : `lsst.daf.butler.Registry`
        The butler registry.
    dataId : `lsst.daf.butler.DataCoordinate`
        The data identifier.
    collections : `list` of `str`
        The collections to search.

    Returns
    -------
    refs : `list` of `lsst.daf.butler.Reference`
        The references to the bias or dark frame.
    """
    return []


def lookupBiasDark(datasetType, registry, dataId, collections):
    """Look up a bias or dark frame

    This is a lookup function for the IsrTaskConnections PrerequisiteConnection
    that finds a bias or dark frame for a given dataId. We use it to provide an
    r or m bias or dark frame for r/m arms.

    Parameters
    ----------
    datasetType : `str`
        The dataset type to look up.
    registry : `lsst.daf.butler.Registry`
        The butler registry.
    dataId : `lsst.daf.butler.DataCoordinate`
        The data identifier.
    collections : `list` of `str`
        The collections to search.

    Returns
    -------
    refs : `list` of `lsst.daf.butler.Reference`
        The references to the bias or dark frame.
    """
    timespan = dataId.timespan
    if dataId["arm"] not in "rm":
        return [registry.findDataset(datasetType, collections=collections, dataId=dataId, timespan=timespan)]

    reducedGraph = DimensionGraph(dataId.universe, names=("instrument", "spectrograph"))
    reducedId = dataId.subset(reducedGraph)
    arm = dataId["arm"]
    spectrograph = dataId["spectrograph"]

    ref1 = None
    ref2 = None
    kwargs = dict(
        collections=collections, dataId=reducedId, spectrograph=spectrograph, timespan=timespan
    )
    try:
        ref1 = registry.findDataset(datasetType, arm=arm, **kwargs)
    except Exception:
        pass
    if arm in "rm":
        otherArm = dict(r="m", m="r")[arm]
        try:
            ref2 = registry.findDataset(datasetType, arm=otherArm, **kwargs)
        except Exception:
            pass

    if ref1 is None and ref2 is None:
        raise RuntimeError(f"Unable to find {datasetType} for {dataId}")
    if ref1 is None:
        return [ref2]
    return [ref1]  # Same arm as original, so use this.


class PfsIsrConnections(PipelineTaskConnections, dimensions=("instrument", "visit", "arm", "spectrograph")):
    """Connections for IsrTask"""

    ccdExposure = InputConnection(
        name="raw.exposure",
        doc="Input exposure to process.",
        storageClass="Exposure",
        dimensions=["instrument", "visit", "arm", "spectrograph"],
    )
    camera = PrerequisiteConnection(
        name="camera",
        storageClass="Camera",
        doc="Input camera to construct complete exposures.",
        dimensions=["instrument"],
        isCalibration=True,
    )

    crosstalk = PrerequisiteConnection(
        name="crosstalk",
        doc="Input crosstalk object",
        storageClass="CrosstalkCalib",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        minimum=0,  # can fall back to cameraGeom
    )
    crosstalkSources = PrerequisiteConnection(
        name="isrOverscanCorrected",
        doc="Overscan corrected input images.",
        storageClass="Exposure",
        dimensions=["instrument", "visit", "arm", "spectrograph"],
        deferLoad=True,
        multiple=True,
        lookupFunction=lookupCrosstalkSources,
        minimum=0,  # not needed for all instruments, no config to control this
    )
    bias = PrerequisiteConnection(
        name="bias",
        doc="Input bias calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        minimum=0,
        lookupFunction=lookupBiasDark,
    )
    dark = PrerequisiteConnection(
        name="dark",
        doc="Input dark calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        lookupFunction=lookupBiasDark,
    )
    flat = PrerequisiteConnection(
        name="fiberFlat",
        doc="Combined flat",
        storageClass="ExposureF",
        dimensions=("instrument", "arm", "spectrograph"),
        isCalibration=True,
    )
    ptc = PrerequisiteConnection(
        name="ptc",
        doc="Input Photon Transfer Curve dataset",
        storageClass="PhotonTransferCurveDataset",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
    )
    fringes = PrerequisiteConnection(
        name="fringe",
        doc="Input fringe calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "physical_filter", "arm", "spectrograph"],
        isCalibration=True,
        minimum=0,  # only needed for some bands, even when enabled
    )
    strayLightData = PrerequisiteConnection(
        name='yBackground',
        doc="Input stray light calibration.",
        storageClass="StrayLightData",
        dimensions=["instrument", "physical_filter", "arm", "spectrograph"],
        deferLoad=True,
        isCalibration=True,
        minimum=0,  # only needed for some bands, even when enabled
    )
    bfKernel = PrerequisiteConnection(
        name='bfKernel',
        doc="Input brighter-fatter kernel.",
        storageClass="NumpyArray",
        dimensions=["instrument"],
        isCalibration=True,
        minimum=0,  # can use either bfKernel or newBFKernel
    )
    newBFKernel = PrerequisiteConnection(
        name='brighterFatterKernel',
        doc="Newer complete kernel + gain solutions.",
        storageClass="BrighterFatterKernel",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        minimum=0,  # can use either bfKernel or newBFKernel
    )
    defects = PrerequisiteConnection(
        name='defects',
        doc="Input defect tables.",
        storageClass="Defects",
        dimensions=["instrument", "detector", "arm", "spectrograph"],
        isCalibration=True,
        lookupFunction=lookupDefects,
    )
    linearizer = PrerequisiteConnection(
        name='linearizer',
        storageClass="Linearizer",
        doc="Linearity correction calibration.",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        minimum=0,  # can fall back to cameraGeom
    )
    opticsTransmission = PrerequisiteConnection(
        name="transmission_optics",
        storageClass="TransmissionCurve",
        doc="Transmission curve due to the optics.",
        dimensions=["instrument"],
        isCalibration=True,
    )
    filterTransmission = PrerequisiteConnection(
        name="transmission_filter",
        storageClass="TransmissionCurve",
        doc="Transmission curve due to the filter.",
        dimensions=["instrument", "physical_filter"],
        isCalibration=True,
    )
    sensorTransmission = PrerequisiteConnection(
        name="transmission_sensor",
        storageClass="TransmissionCurve",
        doc="Transmission curve due to the sensor.",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
    )
    atmosphereTransmission = PrerequisiteConnection(
        name="transmission_atmosphere",
        storageClass="TransmissionCurve",
        doc="Transmission curve due to the atmosphere.",
        dimensions=["instrument"],
        isCalibration=True,
    )
    illumMaskedImage = PrerequisiteConnection(
        name="illum",
        doc="Input illumination correction.",
        storageClass="MaskedImageF",
        dimensions=["instrument", "physical_filter", "arm", "spectrograph"],
        isCalibration=True,
    )
    deferredChargeCalib = PrerequisiteConnection(
        name="cpCtiCalib",
        doc="Deferred charge/CTI correction dataset.",
        storageClass="IsrCalib",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
    )

    outputExposure = OutputConnection(
        name='postISRCCD',
        doc="Output ISR processed exposure.",
        storageClass="Exposure",
        dimensions=["instrument", "visit", "arm", "spectrograph"],
    )
    preInterpExposure = OutputConnection(
        name='preInterpISRCCD',
        doc="Output ISR processed exposure, with pixels left uninterpolated.",
        storageClass="ExposureF",
        dimensions=["instrument", "visit", "arm", "spectrograph"],
    )
    outputOssThumbnail = OutputConnection(
        name="OssThumb",
        doc="Output Overscan-subtracted thumbnail image.",
        storageClass="Thumbnail",
        dimensions=["instrument", "visit", "arm", "spectrograph"],
    )
    outputFlattenedThumbnail = OutputConnection(
        name="FlattenedThumb",
        doc="Output flat-corrected thumbnail image.",
        storageClass="Thumbnail",
        dimensions=["instrument", "visit", "arm", "spectrograph"],
    )
    outputStatistics = OutputConnection(
        name="isrStatistics",
        doc="Output of additional statistics table.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "visit", "arm", "spectrograph"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if config.doBias is not True:
            self.prerequisiteInputs.remove("bias")
        if config.doLinearize is not True:
            self.prerequisiteInputs.remove("linearizer")
        if config.doCrosstalk is not True:
            self.prerequisiteInputs.remove("crosstalkSources")
            self.prerequisiteInputs.remove("crosstalk")
        if config.doBrighterFatter is not True:
            self.prerequisiteInputs.remove("bfKernel")
            self.prerequisiteInputs.remove("newBFKernel")
        if config.doDefect is not True:
            self.prerequisiteInputs.remove("defects")
        if config.doDark is not True:
            self.prerequisiteInputs.remove("dark")
        if config.doFlat is not True:
            self.prerequisiteInputs.remove("flat")
        if config.doFringe is not True:
            self.prerequisiteInputs.remove("fringes")
        if config.doStrayLight is not True:
            self.prerequisiteInputs.remove("strayLightData")
        if config.usePtcGains is not True and config.usePtcReadNoise is not True:
            self.prerequisiteInputs.remove("ptc")
        if config.doAttachTransmissionCurve is not True:
            self.prerequisiteInputs.remove("opticsTransmission")
            self.prerequisiteInputs.remove("filterTransmission")
            self.prerequisiteInputs.remove("sensorTransmission")
            self.prerequisiteInputs.remove("atmosphereTransmission")
        else:
            if config.doUseOpticsTransmission is not True:
                self.prerequisiteInputs.remove("opticsTransmission")
            if config.doUseFilterTransmission is not True:
                self.prerequisiteInputs.remove("filterTransmission")
            if config.doUseSensorTransmission is not True:
                self.prerequisiteInputs.remove("sensorTransmission")
            if config.doUseAtmosphereTransmission is not True:
                self.prerequisiteInputs.remove("atmosphereTransmission")
        if config.doIlluminationCorrection is not True:
            self.prerequisiteInputs.remove("illumMaskedImage")
        if config.doDeferredCharge is not True:
            self.prerequisiteInputs.remove("deferredChargeCalib")

        if config.doWrite is not True:
            self.outputs.remove("outputExposure")
            self.outputs.remove("preInterpExposure")
            self.outputs.remove("outputFlattenedThumbnail")
            self.outputs.remove("outputOssThumbnail")
            self.outputs.remove("outputStatistics")

        if config.doSaveInterpPixels is not True:
            self.outputs.remove("preInterpExposure")
        if config.qa.doThumbnailOss is not True:
            self.outputs.remove("outputOssThumbnail")
        if config.qa.doThumbnailFlattened is not True:
            self.outputs.remove("outputFlattenedThumbnail")
        if config.doCalculateStatistics is not True:
            self.outputs.remove("outputStatistics")

    def adjustQuantum(self, inputs, outputs, label, dataId):
        """Adjust the quantum graph to remove arm=n biases"""
        expRefs = inputs["ccdExposure"][1][0]
        arm = expRefs.dataId["arm"]
        adjusted = {}
        if arm == "n":
            # We don't want the bias for arm=n
            connection, oldRefs = inputs["bias"]
            newRefs = [ref for ref in oldRefs if ref.dataId["arm"] != "n"]
            adjusted["bias"] = (connection, newRefs)
            inputs.update(adjusted)
        super().adjustQuantum(inputs, outputs, label, dataId)
        return adjusted, {}


class PfsIsrTaskConfig(ipIsr.IsrTaskConfig, pipelineConnections=PfsIsrConnections):
    """Configuration parameters for PFS's IsrTask.

    Items are grouped in the order in which they are executed by the task.
    """
    doBrokenRedShutter = pexConfig.Field(dtype=bool, default=False,
                                         doc="Attempt to correct for a broken red shutter?")
    brokenRedShutter = pexConfig.ConfigField(
        dtype=BrokenShutterConfig, doc="Broken shutter related configuration options."
    )
    windowed = pexConfig.Field(dtype=bool, default=False,
                               doc="Reduction of windowed data, for real-time acquisition? Implies "
                               "overscanFitType=MEDIAN")
    darkBBoxes = pexConfig.DictField(
        keytype=str, itemtype=str,
        default=dict(r3="[geom.BoxI(geom.PointI(4045, 2660), geom.PointI(4095, 2910))]"),
        doc="List of BBoxes specifying pixels used to determine amplitude of dark frame")
    fitDarkClipPercentiles = pexConfig.ListField(dtype=float, default=[0, 5],
                                                 doc="""\
Percentages used to clip data when using darkBBoxes to estimate dark amplitude; clipping
is between value and 100 - value.

The first value should probably always be zero, as we haven't removed any signal at that point,
but if you have a sufficiently large cosmic ray flux you might want to reconsider.""")
    crosstalk = pexConfig.ConfigurableField(target=PfsCrosstalkTask, doc="Inter-CCD crosstalk correction")
    doIPC = pexConfig.Field(dtype=bool, default=True, doc="Correct for IPC?")
    ipc = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        default=dict(
            n1="pfsIpc-2023-04-17T13:21:09.707-090750-n1.fits",
            n2="pfsIpc-2023-11-27T13:21:09.707-102106-n2.fits",
            n3="pfsIpc-2023-07-15T00:10:01.950-096714-n3.fits",
            n4="pfsIpc-2024-07-09T00:10:01.950-112241-n4.fits",
        ),
        doc="Mapping of detector name to IPC kernel filename",
    )

    def setDefaults(self):
        super().setDefaults()
        self.doLinearize = False
        self.doCrosstalk = False
        self.doBrighterFatter = False
        self.doFringe = False
        self.doStrayLight = False
        self.assembleCcd.retarget(PfsAssembleCcdTask)
        # self.doApplyGains must be False: gains are applied automatically in PFS's overload of assembleCcd
        # for CCDs and in std_raw for H4RGs (as the ASIC gain can be changed)
        self.doApplyGains = False

    def validate(self):
        if self.windowed:
            self.overscanFitType = "MEDIAN"
        if not self.doSaturationInterpolation and "SAT" in self.maskListToInterpolate:  # fixed on LSST master
            self.maskListToInterpolate.remove("SAT")
        super().validate()
        if self.doApplyGains:
            raise RuntimeError("doApplyGains must be False: gains are applied automatically.")

        for detName, bboxStr in self.darkBBoxes.items():
            try:
                eval(bboxStr)           # can't set self.darkBBoxes[detName] here as it fails type validation
                estr = ""
            except Exception as e:
                estr = str(e)           # it's clearer to the user if we don't raise within the try block

            if estr:
                raise ValueError("Malformed isr.darkBBoxes for %s \"%s\": %s" % (detName, bboxStr, estr))
#
# Code to handle the broken shutter matrix matrix inverse cache
#


class BrokenShutterKernelCache:
    __cache = {}
    __maxsize = 0
    log = lsst.log.getLogger("BrokenShutterKernelCache")

    def __init__(self, maxsize=0):
        self.setCacheSize(maxsize)

    @classmethod
    def compute(cls, mat, t_exp):
        """Compute and cache a new kernel

        Kernels are cached and looked up by ``int(t_exp + 0.5)``

        Parameters:
           mat : `np.array`
             The matrix to be inverted to give a new kernel
           t_exp : `float`
             The exposure time associated with `mat`

        Return:
           The new kernel
        """
        key = int(t_exp + 0.5)
        if key not in cls.__cache:
            if cls.__maxsize > 0 and len(cls.__cache) >= cls.__maxsize:  # need to delete an old kernel
                k_min = None                                            # the least popular kernel
                n_min = 0
                for k, v in cls.__cache.items():
                    if k_min is None or v[1] < n_min:
                        k_min = k
                        n_min = v[1]

                if k_min is not None:
                    cls.log.info(f"Clearing cached kernel for key {k_min}")
                    del cls.__cache[k_min]

            cls.log.info(f"Computing ikern for key {key} ({t_exp}s)")
            ikern = np.linalg.inv(mat)
            cls.__cache[key] = [ikern, 0]

        cls.__cache[key][1] += 1
        return cls.__cache[key][0]

    @classmethod
    def setCacheSize(cls, maxsize):
        """Set the maximum number of cached kernels

        N.b. Kernels are not freed until you next
        compute a new one that doesn't fit in the cache;
        see also clear()

        Parameters:
        maxsize : `int`
           The maximum number of saved kernels

        Returns:
           The old value of maxsize
        """
        old = cls.__maxsize
        cls.__maxsize = maxsize

        return old

    @classmethod
    def clear(cls):
        """Empty the cache"""
        cls.__cache = {}

    @classmethod
    def show(cls):
        """Print some statistics about the cache"""
        print("BrokenShutterKernelCache:")
        print(f"Cache size: {'unlimited' if cls.__maxsize <= 0 else cls.__maxsize}")
        print("key   nUsed")
        for k, c in cls.__cache.items():
            print(f"{k:4}  {c[1]}")


try:
    brokenShutterKernelCache
except NameError:
    brokenShutterKernelCache = BrokenShutterKernelCache()


class PfsIsrTask(ipIsr.IsrTask):
    ConfigClass = PfsIsrTaskConfig
    _DefaultName = "isr"

    """Apply common instrument signature correction algorithms to a raw frame.

    The run method calls runCCD() to use the vanilla ISR for the b/r CCDs, with the
    exception of optionally correcting for the effects of a broken red-arm shutter
    on a PFS spectrograph.

    H4RGs are treated very differently as they are very different beasts; see runH4RG.

    Methods have been split into subtasks that can be redirected appropriately.
    """

    def ensureExposure(self, inputExp, camera=None, detectorNum=None):
        if isinstance(inputExp, afwImage.DecoratedImageU):
            inputExp = afwImage.makeExposure(afwImage.makeMaskedImage(inputExp))
        elif isinstance(inputExp, afwImage.ImageF):
            inputExp = afwImage.makeExposure(afwImage.makeMaskedImage(inputExp))
        elif isinstance(inputExp, afwImage.MaskedImageF):
            inputExp = afwImage.makeExposure(inputExp)
        elif isinstance(inputExp, afwImage.Exposure):
            pass
        elif inputExp is None:
            # Assume this will be caught by the setup if it is a problem.
            return inputExp
        else:
            raise TypeError("Input Exposure is not known type in isrTask.ensureExposure: %s." %
                            (type(inputExp), ))
        return inputExp

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        assert self.config.doCrosstalk is False
        assert self.config.doLinearize is False

        if self.config.doDefect is True:
            if "defects" in inputs and inputs['defects'] is not None:
                # defects is loaded as a BaseCatalog with columns
                # x0, y0, width, height. Masking expects a list of defects
                # defined by their bounding box
                if not isinstance(inputs["defects"], Defects):
                    inputs["defects"] = Defects.fromTable(inputs["defects"])

        assert self.config.doBrighterFatter is False

        assert self.config.doFringe is False
        inputs['fringes'] = Struct(fringes=None)

        assert self.config.doStrayLight is False

        if self.config.doHeaderProvenance:
            # Add calibration provenanace info to header.
            exposureMetadata = inputs['ccdExposure'].getMetadata()
            for inputName in sorted(inputs.keys()):
                reference = getattr(inputRefs, inputName, None)
                if reference is not None and hasattr(reference, "run"):
                    runKey = f"LSST CALIB RUN {inputName.upper()}"
                    runValue = reference.run
                    idKey = f"LSST CALIB UUID {inputName.upper()}"
                    idValue = str(reference.id)

                    exposureMetadata[runKey] = runValue
                    exposureMetadata[idKey] = idValue

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def runDataRef(self, sensorRef):
        """Perform instrument signature removal on a ButlerDataRef of a Sensor.

        This method contains the `CmdLineTask` interface to the ISR
        processing.  All IO is handled here, freeing the `run()` method
        to manage only pixel-level calculations.  The steps performed
        are:
        - Read in necessary detrending/isr/calibration data.
        - Process raw exposure in `run()`.
        - Persist the ISR-corrected exposure as "postISRCCD" if
          config.doWrite=True.

        Parameters
        ----------
        sensorRef : `daf.persistence.butlerSubset.ButlerDataRef`
            DataRef of the detector data to be processed

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with component:
            - ``exposure`` : `afw.image.Exposure`
                The fully ISR corrected exposure.

        Raises
        ------
        RuntimeError
            Raised if a configuration option is set to True, but the
            required calibration data does not exist.

        """
        self.log.info("Performing ISR on sensor %s" % (sensorRef.dataId))

        ccdExposure = sensorRef.get(self.config.datasetType)
        camera = sensorRef.get("camera")
        if camera is None and self.config.doAddDistortionModel:
            raise RuntimeError("config.doAddDistortionModel is True "
                               "but could not get a camera from the butler")
        isrData = self.readIsrData(sensorRef, ccdExposure)

        result = self.run(ccdExposure, camera=camera, **isrData.getDict())
        result.exposure.getMetadata().set("PFSCALIB", getCalibPath(sensorRef))

        if self.config.doBrokenRedShutter and \
           ccdExposure.getDetector().getName() in self.config.brokenRedShutter.brokenShutterList:
            if not self.config.brokenRedShutter.useAnalytic:
                #
                # If we're using the `r` arm subtract a neighbouring image if it's a bias
                #
                # If the source of photons is stable, this corrects for the extra
                # light accumulated during CCD wipe and readout
                #
                # We search for the bias in a window around the visit v0 in order
                #   v0 + 1, v0 - 1, v0 + 2, v0 - 2, ...
                # (if brokenRedShutter.causalWindow is True, only use v0 + 1, v0 + 2, ...)
                #
                # Don't fix a bias; this is especially important as if we do, it'll lead
                # to an attempt to correct the bias we're using for the shutter correction
                #
                # N.b. we can't rely on the sensorRef.dataId because we adjusted the visit number
                # before possibly calling this routine recursively to fetch the desired bias
                md = sensorRef.get("raw_md")
                if md.get('DATA-TYP', "").upper() != "BIAS":
                    dataId = dict(arm=sensorRef.dataId["arm"],
                                  filter=sensorRef.dataId["arm"],  # yes, this is necessary
                                  spectrograph=sensorRef.dataId["spectrograph"])
                    dataId0 = sensorRef.dataId

                    for dv in sum([[i, -i] for i in range(1, self.config.brokenRedShutter.window + 1)], []):
                        if self.config.brokenRedShutter.causalWindow and dv < 0:
                            continue

                        exp2 = None
                        try:
                            sensorRef.dataId = dataId
                            dataId["visit"] = dataId0["visit"] + dv

                            md = sensorRef.get("raw_md")
                            if md.get('DATA-TYP', "").upper() == "BIAS":
                                self.log.info("Reading visit %d for broken red shutter correction",
                                              sensorRef.dataId["visit"])

                                exp2 = self.runDataRef(sensorRef).exposure
                        except Exception as e:
                            self.log.warn("Unable to read %d %s%d for broken red shutter correction: %s",
                                          sensorRef.dataId["visit"], sensorRef.dataId["arm"],
                                          sensorRef.dataId["spectrograph"], e)
                        finally:
                            sensorRef.dataId = dataId0

                        if exp2:
                            self.log.info("Correcting for broken red shutter using visit %d", dataId["visit"])
                            result.exposure.maskedImage -= exp2.maskedImage

                            break

                    if exp2 is None:
                        self.log.warn("Unable to find bias frame for broken red shutter correction %s",
                                      sensorRef.dataId)

        if self.config.doWrite:
            sensorRef.put(result.exposure, "postISRCCD")
        if result.ossThumb is not None:
            isrQa.writeThumbnail(sensorRef, result.ossThumb, "ossThumb")
        if result.flattenedThumb is not None:
            isrQa.writeThumbnail(sensorRef, result.flattenedThumb, "flattenedThumb")

        return result

    def readIsrData(self, dataRef, rawExposure):
        """Retrieve necessary frames for instrument signature removal.

        Pre-fetching all required ISR data products limits the IO
        required by the ISR. Any conflict between the calibration data
        available and that needed for ISR is also detected prior to
        doing processing, allowing it to fail quickly.

        Parameters
        ----------
        dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
            Butler reference of the detector data to be processed
        rawExposure : `afw.image.Exposure`
            The raw exposure that will later be corrected with the
            retrieved calibration data; should not be modified in this
            method.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Struct containing the necessary calibs.
        """
        result = super().readIsrData(dataRef, rawExposure)

        if dataRef.dataId["arm"] == "n" and self.config.doIPC:
            try:
                dateObs = rawExposure.getInfo().getVisitInfo().getDate().toPython().isoformat()
            except RuntimeError:
                dateObs = None

            if self.config.doIPC:
                result.ipcCoeffs = self.getIsrExposure(dataRef,
                                                       self.config.ipcDataProductName, dateObs=dateObs)

        return result

    def getIsrExposure(self, dataRef, datasetType, dateObs=None, immediate=True):
        """Retrieve a calibration dataset for removing instrument signature.

        This override refuses to load biases for ``arm=n``: we don't
        use them in ``runH4RG``.

        Parameters
        ----------
        dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
            DataRef of the detector data to find calibration datasets
            for.
        datasetType : `str`
            Type of dataset to retrieve (e.g. 'bias', 'flat', etc).
        dateObs : `str`, optional
            Date of the observation.  Used to correct butler failures
            when using fallback filters.
        immediate : `Bool`
            If True, disable butler proxies to enable error handling
            within this routine.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure` or `None`
            Requested calibration frame, or `None` if not required.

        Raises
        ------
        RuntimeError
            Raised if no matching calibration frame can be found.
        """
        arm = dataRef.dataId["arm"]
        if arm == "n" and datasetType in ("bias",):
            return None
        return super().getIsrExposure(dataRef, datasetType, dateObs=dateObs, immediate=immediate)

    def run(self, ccdExposure, **kwargs):
        """Perform instrument signature removal on an exposure.

        Parameters:
           ccdExposure : `afwImage.Exposure`
             The raw data to be processed
           kwargs : `dict`
             Dict of extra parameters, e.g. combined bias

        Return:
           result : `lsst.pipe.base.Struct`
              Result struct;  see `lsst.ip.isr.isrTask.run`
        """
        exposure = ccdExposure                         # argument must be called "ccdExposure"; PIPE2D-1093

        if exposure.getFilter().bandLabel == 'n':  # treat H4RGs specially
            return self.runH4RG(exposure, **kwargs)
        else:
            return self.runCCD(exposure, **kwargs)

    def runCCD(self, ccdExposure, **kwargs):
        """Perform instrument signature removal on a CCD exposure.

        Parameters:
           exposure : `afwImage.Exposure`
             The raw data to be processed
           kwargs : `dict`
             Dict of extra parameters, e.g. combined bias

        Return:
           result : `lsst.pipe.base.Struct`
              Result struct;  see `lsst.ip.isr.isrTask.run`
        """
        if ccdExposure.getFilter().bandLabel == 'n':  # treat H4RGs specially
            return self.runH4RG(ccdExposure, **kwargs)

        doBrokenRedShutter = self.config.doBrokenRedShutter and \
            ccdExposure.getDetector().getName() in self.config.brokenRedShutter.brokenShutterList

        if doBrokenRedShutter and self.config.brokenRedShutter.checkParallelOverscan:
            excessPOscan = []
            meanFlux = []
            for amp in ccdExposure.getDetector():
                # N.b. MEANCLIP is broken for U16 data
                data = afwMath.makeStatistics(ccdExposure[amp.getRawDataBBox()].convertF().image,
                                              afwMath.MEANCLIP).getValue()
                pOscan = afwMath.makeStatistics(ccdExposure[amp.getRawVerticalOverscanBBox()].image,
                                                afwMath.MEDIAN).getValue()
                sOscan = afwMath.makeStatistics(ccdExposure[amp.getRawHorizontalOverscanBBox()].image,
                                                afwMath.MEDIAN).getValue()

                meanFlux.append(data)
                excessPOscan.append(pOscan - sOscan)

            flux = np.mean(meanFlux)
            excess = np.mean(excessPOscan)
            if excess <= self.config.brokenRedShutter.maximumAllowedParallelOverscanFraction*flux:
                doBrokenRedShutter = False

            self.log.info("Checking parallel overscan: %g %s %g*%g %s",
                          excess, (">" if doBrokenRedShutter else "<="),
                          self.config.brokenRedShutter.maximumAllowedParallelOverscanFraction, flux,
                          ("" if doBrokenRedShutter else "; assuming not broken"))

        detName = ccdExposure.getDetector().getName()
        darkBBoxes = self.config.darkBBoxes
        if self.config.doDark and detName in darkBBoxes:
            kwargs0 = kwargs            # initial kwargs
            kwargs = kwargs.copy()
            kwargs["dark"] = kwargs["dark"].clone()
            kwargs["dark"].maskedImage *= 0         # disable dark correction

        results = super().run(ccdExposure, **kwargs)
        exposure = results.exposure

        isNan = np.isnan(exposure.image.array) | np.isnan(exposure.variance.array)
        if np.any(isNan):
            self.log.warn("Unmasked NaNs in ISR-processed exposure: %d pixels", np.sum(isNan))
            exposure.mask.array[isNan] |= exposure.mask.getPlaneBitMask(["BAD", "UNMASKEDNAN"])

        if self.config.doDark and detName in darkBBoxes:
            bboxes = eval(darkBBoxes[detName])  # we checked that this is OK in PfsIsrTaskConfig.validate()

            scale = self.darkCorrectionFromBBoxes(bboxes, exposure.maskedImage, kwargs0["dark"].maskedImage)

            darktime = exposure.info.getVisitInfo().getDarkTime()
            self.log.info("Scaled dark exposure by %.1f (%.3f/second)", scale, scale/darktime)

        if doBrokenRedShutter and self.config.brokenRedShutter.useAnalytic:
            #
            # Build a model which generates the as-readout data R given
            # the (unknown) true signal S, with integration time t_exp:
            #    R_i = S_i + (t_wipe*S_red + t_stare*S_i + t_read*S_blue,i)/t_exp
            # where the second term gives the photons resulting from reading
            # the chip with the shutter open, and S_blue/red are the photons bluer (redder)
            # than S_i (i.e. sum_0^i S_i and sum_i^N-1 S_i)
            #
            # We write this as R = S + kern@S, or
            #  S = (1 + kern)^{-1} R
            #
            # N.b. we could cache ikern as it only depends on the exposure time
            # and the config parameters specifying times for the readout phases
            #
            t_exp = exposure.getInfo().getVisitInfo().getExposureTime()
            t_wipe = self.config.brokenRedShutter.t_wipe
            t_read = self.config.brokenRedShutter.t_read
            t_stare = self.config.brokenRedShutter.t_stare

            if t_exp == 0:          # we can't correct a bias
                self.log.debug("Not correcting bias for broken red shutter analytically")
            else:
                brokenShutterKernelCache.setCacheSize(self.config.brokenRedShutter.cacheSize)

                N = exposure.getHeight()
                kern = np.zeros((N, N))
                kern += t_wipe/N*(1 - np.tri(N, k=1))  # i.e. upper tri
                kern += t_read/N*np.tri(N, k=-1)
                np.fill_diagonal(kern, t_stare)
                kern /= t_exp

                np.fill_diagonal(kern, 1 + kern.diagonal())  # i.e. add 1 to the diagonal

                ikern = brokenShutterKernelCache.compute(kern, t_exp)

                if not self.config.brokenRedShutter.doCacheCorrections:
                    brokenShutterKernelCache.clear()

                self.log.info("Correcting for broken red shutter analytically")
                exposure.image.array[:] = ikern@exposure.image.array
                exposure.variance.array[:] = (ikern**2)@exposure.variance.array

        return results

    def runH4RG(self, exposure, dark=None, flat=None, defects=None, detectorNum=None, **kwargs):
        """Specialist instrument signature removal for H4RG detectors

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The raw exposure that is to be run through ISR.  The
            exposure is modified by this method.
        dark : `lsst.afw.image.Exposure`, optional
            Dark exposure to subtract.
        flat : `lsst.afw.image.Exposure`, optional
            Flat-field exposure to divide by.
        defects : `lsst.ip.isr.Defects`, optional
            List of defects.
        detectorNum: `int`, optional
            The integer number for the detector to process.
        kwargs : `dict` other keyword arguments specifying e.g. combined biases
            N.b. no values are currently valid
        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with component:
            - ``exposure`` : `afw.image.Exposure`
                The fully ISR corrected exposure.
            - ``outputExposure`` : `afw.image.Exposure`
                An alias for `exposure`
        """
        for k, v in kwargs.items():
            if v is not None:
                if k == "fringes" and v.fringes is None:
                    continue

                if k not in ("camera", "crosstalk", "crosstalkSources", "linearizer"):
                    self.log.warn("Unexpected argument for runH4RG: %s", k)

        if self.config.doDark and dark is None:
            raise RuntimeError("Must supply a dark exposure if config.doDark=True.")
        if self.config.doFlat and flat is None:
            raise RuntimeError("Must supply a flat exposure if config.doFlat=True.")
        if self.config.doDefect and defects is None:
            raise RuntimeError("Must supply defects if config.doDefect=True.")
        if self.config.doIPC:
            if defects is None:
                raise RuntimeError("Must supply defects if config.doIPC=True.")
            ipcCoeffs = self.readIPC(exposure.getDetector().getName())

        # Think about the ordering here.
        # E.g. dark subtraction before defects (but then we'd have to rotate the dark,
        # either here or when writing it)
        assert len(exposure.getDetector()) == 1, "Fix me now we have multiple channels"

        channel = exposure.getDetector()[0]
        gain = channel.getGain()

        exposure.image *= gain        # convert to electrons
        var = exposure.image.clone()  # assumes photon noise -- not true for the persistence
        var += 2*channel.getReadNoise()**2  # 2* comes from CDS
        exposure.variance = var

        nQuarter = exposure.getDetector().getOrientation().getNQuarter()

        if self.config.doIPC:
            self.log.info("Applying IPC correction.")
            self.correctIPC(exposure, defects, ipcCoeffs, -nQuarter)

        if self.config.doDefect:
            super().maskAndInterpolateDefects(exposure, defects)

        if self.config.maskNegativeVariance:
            super().maskNegativeVariance(exposure)

        if nQuarter != 0:
            exposure.maskedImage = afwMath.rotateImageBy90(exposure.maskedImage, nQuarter)

        if self.config.doDark:
            self.log.info("Applying dark correction.")
            super().darkCorrection(exposure, dark)
            super().debugView(exposure, "doDark")

        if self.config.doFlat:
            self.log.info("Applying flat correction.")
            self.flatCorrection(exposure, flat)
            self.debugView(exposure, "doFlat")

        return pipeBase.Struct(exposure=exposure,
                               outputExposure=exposure,  # is this needed? Cargo culted from ip_isr isrTask.py
                               flattenedThumb=None,
                               ossThumb=None,
                               )

    def readIPC(self, detectorName: str) -> dict[int, dict[int, float]]:
        """Read IPC coefficients from file

        We believe these should be static, so haven't gone to the trouble of
        making them calibs. Instead, they are read from a file specified in the
        config, indexed by detector name.

        Parameters
        ----------
        detectorName : `str`
            Name of detector.

        Returns
        ipcCoeffs : `dict` mapping `int` to a `dict` mapping `int` to float
            IPC coefficients. ``ipcCoeffs[dx][dy]`` is the fraction of flux at
            (dx, dy).
        """
        filename = self.config.ipc[detectorName]
        if not os.path.isabs(filename):
            absFilename = os.path.join(getPackageDir("drp_pfs_data"), "ipc", filename)
            if os.path.exists(absFilename):
                filename = absFilename

        data = DecoratedImageF(filename)
        protocol = data.getMetadata()["PROTOCOL"]
        if protocol != 1:
            raise RuntimeError(f"Unexpected IPC protocol version {protocol}")

        arr = data.image.array
        nx, ny = arr.shape

        ipcCoeffs = {}
        for ix in range(-(nx//2), nx//2 + 1):
            ipcCoeffs[ix] = {}
            for iy in range(-(ny//2), ny//2 + 1):
                ipcCoeffs[ix][iy] = arr[iy + ny//2, ix + nx//2]

        return ipcCoeffs

    @staticmethod
    def correctIPC(exposure, defects, ipcCoeffs, nQuarter=0):
        """Correct the exposure for the IPC associated with the defects

        Note that the coefficients may be scalars or images with the same dimensions as the exposure

        exposure:  The image in question
        defects:   List of defects
        ipcCoeffs: dict of dict of coefficients;  ipcCoeffs[dx][dy] is the fraction of flux at (dx, dy)
        nQuarter:  number of pi/2 turns that will be applied to the data to get it into the
                   same orientation as the IPC coefficients
        """
        ipc = exposure.maskedImage.clone()
        ipc.mask[:] = 0x0
        defects.maskPixels(ipc)

        ipcarr = ipc.image.array
        ipcarr[ipc.mask.array == 0] = 0
        ipcarr[np.isnan(ipcarr)] = 0

        ipcmodel = np.zeros_like(ipcarr, dtype=np.float32)

        nx, ny = exposure.getDimensions()
        nc = len(ipcCoeffs)
        nc2 = nc//2
        for y in range(-nc2, nc2 + 1):
            y0, y1 = (y, ny)     if y > 0 else ( 0, ny + y)  # noqa E201, E272
            y2, y3 = (0, ny - y) if y > 0 else (-y, ny)
            for x in range(-nc2, nc2 + 1):
                if x == 0 and y == 0:
                    continue

                x0, x1 = (x, nx)     if x > 0 else ( 0, nx + x)  # noqa E201, E272
                x2, x3 = (0, nx - x) if x > 0 else (-x, nx)
                #
                # The image is rotated by nQuarter relative to the IPC coefficients.
                # Rotate the coefficients to match
                if nQuarter%4 == 0:
                    rx, ry = x, y
                elif nQuarter%4 == 1:
                    rx, ry = y, -x
                elif nQuarter%4 == 2:
                    rx, ry = -x, -y
                else:
                    rx, ry == -y, x

                rx, ry = y, -x
                ipcmodel[y0:y1, x0:x1] += (ipcCoeffs[rx][ry]*ipcarr)[y2:y3, x2:x3]

        exposure.image.array[:, :] -= ipcmodel

        if "IPC" not in exposure.mask.getMaskPlaneDict():
            exposure.mask.addMaskPlane("IPC")
            afwDisplay.setDefaultMaskPlaneColor("IPC", "GREEN")

        exposure.mask.array[ipcmodel != 0] |= exposure.mask.getPlaneBitMask("IPC")

    def maskAmplifier(self, exposure, amp, defects):
        """Mask bad pixels in amplifier

        This PFS override masks pixels with a value of zero. These pixels have
        not really been read out, and only exist in order to pad the image.
        This needs to be done before overscan subtraction, so that the empty
        rows aren't used to measure the overscan.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Raw exposure.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier information.
        defects : `lsst.meas.algorithms.Defects`
            List of defects.

        Returns
        -------
        badAmp : `bool`
            Whether the entire amplifier is covered by defects.
        """
        badAmp = super().maskAmplifier(exposure, amp, defects)
        ampExp = exposure[amp.getRawBBox()]
        zeroData = ampExp.image.array == 0
        if np.any(zeroData):
            # Overscan subtraction only respects SAT
            # NO_DATA signals to our post-assembleCcd code to remove the values
            ampExp.mask.array[zeroData] |= ampExp.mask.getPlaneBitMask(["SAT", "NO_DATA"])

        return badAmp

    def darkCorrectionFromBBoxes(self, bboxes, data, dark):
        """Apply dark correction in place.

        The amplitude is set by a robust estimate from the pixels specified by bboxes,
        assumed to be un-contaminated by non-dark signal

        Parameters
        ----------
        bboxes : `list` of `lsst.geom.Box2I`
           Specify which pixels should be used to estimate amplitude
        data : `lsst.afw.image.MaskedImage`
           Image to process.  The image is modified by this method.
        dark : `lsst.afw.image.MaskedImage`
            Dark image of the same size as ``data``.

        Returns
        -------
        scale : `float`
           The scaling applied to the input dark image
        """

        # Start with a rough linear estimate (well, the MLE ignoring e.g. cosmic rays),
        # then repeat, after clipping out the first and last n-percentiles

        finalScale = 0                  # the scaling we applied to the dark when all iterations have finished
        for i, clip in enumerate(self.config.fitDarkClipPercentiles):
            sumDataDark = 0
            sumDarkDark = 0
            for bbox in bboxes:
                dataArr = data.image[bbox].array
                darkArr = dark.image[bbox].array
                mask = data.mask[bbox].array | dark.mask[bbox].array
                var = data.variance[bbox].array

                keep = (mask & data.mask.getPlaneBitMask(["BAD", "SAT", "NO_DATA"])) == 0x0

                if clip > 0:
                    qa, qb = np.percentile(dataArr, [clip, 100 - clip])
                    keep = np.logical_and(keep, np.logical_and(dataArr > qa, dataArr < qb))

                    self.log.debug("Iteration %d: clipping at %g: %.1f -- %.1f", i, clip, qa, qb)

                ngood = np.sum(keep)
                if ngood == 0:
                    self.log.warn("Iteration %d: There are no good pixels", i)
                    return finalScale

                if ngood < keep.size:
                    dataArr = dataArr[keep]
                    darkArr = darkArr[keep]
                    var = var[keep]

                sumDataDark += np.sum(dataArr*darkArr/var)
                sumDarkDark += np.sum(darkArr**2/var)

            scale = sumDataDark/sumDarkDark
            finalScale += scale

            self.log.debug("Iteration %d: dark scaling %.1f", i, scale)

            data.scaledMinus(scale, dark)

        return finalScale

    def roughZeroPoint(self, exposure):
        """Set an approximate magnitude zero point for the exposure.

        We disable this for PFS, since we don't use zero-points.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        """
        pass

    def _getMetadataName(self):
        return None                     # don't write metadata; requires fix to ip_isr
