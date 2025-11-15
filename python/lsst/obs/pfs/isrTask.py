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
from typing import TYPE_CHECKING, Optional

import os
import time
import warnings

from functools import partial

import numpy as np
import scipy
import ruamel.yaml as yaml

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
from lsst.ip.isr.defects import Defects
from lsst.ip.isr.isrFunctions import updateVariance
import lsst.pipe.base as pipeBase
from lsst.pipe.base.connectionTypes import Input as InputConnection
from lsst.pipe.base.connectionTypes import PrerequisiteInput as PrerequisiteConnection
from lsst.pipe.base.connectionTypes import Output as OutputConnection
from lsst.pipe.base import Struct, PipelineTaskConnections
from lsst.daf.butler import DimensionGroup

from . import imageCube
from . import nirLinearity
from .overscan import PfsOverscanCorrectionTask
from pfs.drp.stella.crosstalk import PfsCrosstalkTask

if TYPE_CHECKING:
    from lsst.afw.image import ExposureF
    from .imageCube import ImageCube
    from .raw import PfsRaw

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


class H4Config(pexConfig.Config):
    """Configuration parameters for H4 reductions"""

    quickCDS = pexConfig.Field(dtype=bool, default=False,
                               doc="Only consider last and first reads instead of the full ramp cube")
    doCR = pexConfig.Field(dtype=bool, default=False,
                           doc="Run ramp-based CR rejection. Requires quickCDS=False")
    crMinReads = pexConfig.Field(dtype=int, default=4,
                                 doc="Minimum number of dark or shutter-open reads for CR rejection")
    useIRP = pexConfig.Field(dtype=bool, default=True,
                             doc="Use Interleaved Reference Pixel planes if available")
    repairAsicSpikes = pexConfig.Field(dtype=bool, default=True,
                                       doc="Detect and mask single-pixel, single-read ASIC glitches")
    repairAsicSpikesSigma = pexConfig.Field(dtype=float, default=3.0,
                                            doc="Sigma threshold from reads for detecting ASIC glitches")
    IRPfilter = pexConfig.Field(dtype=int, default=15,  # 15..31, probably.
                                doc="width of smoothing window for IRP corrections. 0=no smoothing. Odd")
    doIRPbadPixels = pexConfig.Field(dtype=bool, default=True,
                                     doc="Interpolate over known bad IRP row pixels")
    doIRPcrosstalk = pexConfig.Field(dtype=bool, default=False,
                                     doc="Correct for data->IRP crosstalk")

    doIPC = pexConfig.Field(dtype=bool, default=False, doc="Correct for IPC?")
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

    useDarkCube = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Use dark cube for dark subtraction? Disable this to use traditional darks.",
    )
    applyUTRWeights = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Apply UTR weights to the ramp cube? Disable this to use CDS.",
    )
    doWriteRawCube = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="write out raw-ish ISRCube (UTR, e-, with IRP and CR corrections)",
    )
    doLinearize = pexConfig.Field(dtype=bool, default=True, doc="Apply linearity correction?")
    linearizeBeforeUTR = pexConfig.Field(dtype=bool, default=False,
                                         doc="Apply linearity correction before UTR")


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
        metadata = exposure.getMetadata()
        for amp in exposure.getDetector():
            ampImage = exposure.maskedImage[amp.getBBox()]
            ampImage *= amp.getGain()
            metadata[f"PFS AMP{amp.getName()} GAIN ORIG"] = amp.getGain()

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


def lookupBiasDark(datasetType, registry, dataId, collections, allowEmpty=False):
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

    reducedGraph = DimensionGroup(dataId.universe, names=("instrument", "spectrograph"))
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
        if allowEmpty:
            return []
        raise RuntimeError(f"Unable to find {datasetType} for {dataId}")
    if ref1 is None:
        return [ref2]
    return [ref1]  # Same arm as original, so use this.


class PfsIsrConnections(PipelineTaskConnections, dimensions=("instrument", "visit", "arm", "spectrograph")):
    """Connections for IsrTask"""

    ccdExposure = InputConnection(
        name="raw",
        doc="Input exposure to process.",
        storageClass="PfsRaw",
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
        lookupFunction=partial(lookupBiasDark, allowEmpty=True),
        minimum=0,  # allowed to not exist, since we may use nirDark instead
    )
    nirDark = PrerequisiteConnection(
        name="nirDark",
        doc="Input dark calibration for NIR",
        storageClass="ImageCube",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        lookupFunction=partial(lookupBiasDark, allowEmpty=True),
        minimum=0,  # allowed to not exist, since we may use dark instead
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
    outputRawCube = OutputConnection(
        name='rawISRCube',
        doc="Output semi-ISR processed ramp cube.",
        storageClass="ImageCube",
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
            self.prerequisiteInputs.remove("nirDark")
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
            self.outputs.remove("outputRawCube")
            self.outputs.remove("preInterpExposure")
            self.outputs.remove("outputFlattenedThumbnail")
            self.outputs.remove("outputOssThumbnail")
            self.outputs.remove("outputStatistics")
        if config.h4.doWriteRawCube is not True:
            self.outputs.remove("outputRawCube")

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

    h4 = pexConfig.ConfigField(
        dtype=H4Config, doc="H4 related configuration")

    crosstalk = pexConfig.ConfigurableField(target=PfsCrosstalkTask, doc="Inter-CCD crosstalk correction")
    overscan = pexConfig.ConfigurableField(
        target=PfsOverscanCorrectionTask,
        doc="Overscan subtraction task for image segments.",
    )

    def setDefaults(self):
        super().setDefaults()
        self.overscan.fitType = "AKIMA_SPLINE"
        self.overscan.order = 30
        self.fwhm = 2.5
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
        assert self.config.doLinearize is False  # We're not set up to do CCD linearization yet

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

        raw = inputs["ccdExposure"]  # Inheritance requires "ccdExposure", but it's actually a PfsRaw
        isNir = raw.isNir()
        if isNir and self.config.h4.doIPC:
            inputs["ipcCoeffs"] = self.readIPC(raw.detector.getName())

        if self.config.doDark:
            if isNir and self.config.h4.useDarkCube:
                if inputs.get("nirDark") is None:
                    raise RuntimeError(
                        f"No NIR dark cube found for {raw.detector.getName()}; try h4.useDarkCube=False"
                    )
            elif inputs["dark"] is None:
                raise RuntimeError(f"No dark frame found for {raw.detector.getName()}")

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

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
        pfsRaw = ccdExposure  # argument must still be called "ccdExposure"; PIPE2D-1093
        if pfsRaw.isNir():  # treat H4RGs specially
            return self.runH4RG(pfsRaw, **kwargs)
        else:
            return self.runCCD(pfsRaw.getExposure(), **kwargs)

    def updateVariance(self, ampExposure, amp, overscanImage=None, ptcDataset=None):
        """Set the variance plane using the gain and read noise

        This override enforces a gain of unity, since we've already converted
        the data to electrons in ``assembleCcd``. Unfortunately, the detector
        passed to this method is not the updated detector with unit gains, so
        we need to override the gains used here.

        Parameters
        ----------
        ampExposure : `lsst.afw.image.Exposure`
            Exposure to process.
        amp : `lsst.afw.cameraGeom.Amplifier` or `FakeAmp`
            Amplifier detector data.
        overscanImage : `lsst.afw.image.MaskedImage`, optional.
            Image of overscan, required only for empirical read noise.
        ptcDataset : `lsst.ip.isr.PhotonTransferCurveDataset`, optional
            PTC dataset containing the gains and read noise.
        """
        gain = 1.0
        readNoise = amp.getReadNoise()

        metadata = ampExposure.getMetadata()
        metadata[f'LSST GAIN {amp.getName()}'] = gain
        metadata[f'LSST READNOISE {amp.getName()}'] = readNoise
        updateVariance(maskedImage=ampExposure.maskedImage, gain=gain, readNoise=readNoise)

    def runCCD(self, ccdExposure, **kwargs):
        """Perform instrument signature removal on a CCD exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
             The exposure to be processed
        kwargs : `dict`
            Dict of extra parameters, e.g., bias, dark, etc.

        Return
        ------
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

    def runH4RG(
        self,
        pfsRaw,
        dark=None,
        nirDark=None,
        flat=None,
        defects=None,
        detectorNum=None,
        ipcCoeffs=None,
        **kwargs,
    ):
        """Specialist instrument signature removal for H4RG detectors

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
            The raw exposure that is to be run through ISR.  The
            exposure is modified by this method. With the PfsRaw we
            can get access the ramp cubes.
        dark : `lsst.afw.image.Exposure`, optional
            Dark exposure to subtract (if ``config.h4.useDarkCube`` is
            ``False``).
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`, optional
            Dark cube to subtract (if ``config.h4.useDarkCube`` is ``True``).
        flat : `lsst.afw.image.Exposure`, optional
            Flat-field exposure to divide by.
        defects : `lsst.ip.isr.Defects`, optional
            List of defects.
        detectorNum: `int`, optional
            The integer number for the detector to process.
        ipcCoeffs : `dict` of `dict` of `float`, optional
            IPC coefficients.  ipcCoeffs[dx][dy] is the fraction of flux at (dx, dy).
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

                if k in ("crosstalk", "crosstalkSources", "linearizer"):
                    self.log.warn("Unexpected argument for runH4RG: %s", k)

        if self.config.doDark:
            if self.config.h4.useDarkCube:
                if nirDark is None:
                    raise RuntimeError(
                        "Must supply a dark cube if config.doDark=True and config.h4.useDarkCube=True."
                    )
            elif dark is None:
                raise RuntimeError(
                    "Must supply a dark exposure if config.doDark=True and config.h4.useDarkCube=False."
                )
        else:
            nirDark = None
        if self.config.doFlat and flat is None:
            raise RuntimeError("Must supply a flat exposure if config.doFlat=True.")
        if self.config.doDefect and defects is None:
            raise RuntimeError("Must supply defects if config.doDefect=True.")
        if self.config.h4.doIPC:
            if defects is None:
                raise RuntimeError("Must supply defects if config.h4.doIPC=True.")
            if ipcCoeffs is None:
                raise RuntimeError("Must supply IPC coefficients if config.h4.doIPC=True.")
        if self.config.h4.doLinearize:
            linearity = self.resolveNirLinearity(pfsRaw.detector.getName())
            if linearity is None:
                raise RuntimeError("Linearity corrections must be available if config.h4.doLinearize=True.")
            if defects is None:
                self.log.warn('You usually want to supply defects when linearizing. Will avoid worst pixels.')

        # All 3-d ramp-based operations are hidden inside this call. After .makeNirExposure()
        # we operate on a 2-d image.
        exposure, rawRamp = self.makeNirExposure(pfsRaw, nirDark, self.config.h4.doWriteRawCube)

        # We BAD mask the defects now, but do not interpolate.
        if defects is not None:
            super().maskDefect(exposure, defects)
        if self.config.h4.doLinearize:
            self.log.info("Correcting non-linearity.")
            exposure = self.applyNirLinearity(pfsRaw, exposure, linearity)

        assert len(exposure.getDetector()) == 1, "Fix me now we have multiple channels"

        channel = exposure.getDetector()[0]
        gain = channel.getGain()

        exposure.image *= gain  # convert to electrons
        var = exposure.image.array.copy()  # assumes photon noise -- not true for the persistence
        var += 2*channel.getReadNoise()**2  # 2* comes from CDS
        exposure.variance.array[:] = var

        if rawRamp is not None:
            rawRamp *= gain
            rawRamp = imageCube.ImageCube.fromCube(rawRamp, exposure.getMetadata())

        if self.config.doDark and not self.config.h4.useDarkCube:
            self.log.info("Applying simple dark correction.")
            super().darkCorrection(exposure, dark)
            super().debugView(exposure, "doDark")

        # Any nquarter stuff should be removed after PIPE2D-1200
        nQuarter = exposure.getDetector().getOrientation().getNQuarter()

        if self.config.h4.doIPC:
            self.log.info("Applying IPC correction.")
            self.correctIPC(exposure, defects, ipcCoeffs, -nQuarter)

        if self.config.doDefect:
            self.log.info("Masking defects.")
            if not self.config.h4.useDarkCube:
                self.log.warning("Masking defects, but dark cube not used")
            super().maskAndInterpolateDefects(exposure, defects)

        if self.config.maskNegativeVariance:
            super().maskNegativeVariance(exposure)

        if self.config.doFlat:
            self.log.info("Applying flat correction.")
            self.flatCorrection(exposure, flat)
            self.debugView(exposure, "doFlat")

        return pipeBase.Struct(exposure=exposure,
                               outputExposure=exposure,  # is this needed? Cargo culted from ip_isr isrTask.py
                               outputRawCube=rawRamp,
                               flattenedThumb=None,
                               ossThumb=None,
                               )

    def _makeExposure(self, pfsRaw, image):
        """Re-construct an Exposure from an Image.

        Only intended to be used at the end of H4 ramp processing. Uses the PHDU from
        the raw ramp to construct a new one, given the processed 2D image.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
            The raw exposure that we ran through ISR.
        image : `np.ndarray`
            The output of H4 ISR processing.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            The corrected 2-d image
        """
        maskedIm = afwImage.makeMaskedImageFromArrays(image, None, None)
        exposure = afwImage.makeExposure(maskedIm)
        exposure.setDetector(pfsRaw.detector)
        info = exposure.getInfo()
        info.setVisitInfo(pfsRaw.visitInfo)
        info.setId(pfsRaw.visitInfo.id)
        info.setMetadata(pfsRaw.metadata)
        info.setDetector(pfsRaw.detector)
        arm = pfsRaw.obsInfo.ext_arm
        info.setFilter(afwImage.FilterLabel(arm, arm))
        return exposure

    def resolveNirLinearity(self, cam):
        """Get the full path for our linearity corrections.

        Parameters
        ----------
        cam : `str`
           The camera name.

        Returns
        -------
        absFilename : `str`
           The full path of the linearity corrections file, or None if not found.
        """

        filename = f'nirLinearity-{cam}.fits'
        absFilename = os.path.join(getPackageDir("drp_pfs_data"), "nirLinearity", filename)
        if not os.path.exists(absFilename):
            self.log.warn(f'no linearity available for {cam}: {absFilename} not found')
            return None

        return nirLinearity.NirLinearity.readFits(absFilename)

    def nirLinearityChebyshev(self, exposure, linearity):
        """Apply linearity corrections using numpy polynomials.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            exposure to correct.
        linearity : `nirLinearity.NirLinearity`
            all parts of the corrections.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
           the input exposure, corrected in place.
        """

        limits = linearity.limits
        coeffs = linearity.coeffs

        # We only record the max valid range of the linearization.
        # Construct the "domains" and "scales" which numpy needs to
        # normalize the exposure data. [0..maxValid] -> [-1..1]
        domains = np.zeros(shape=(2, 4096, 4096), dtype='f4')
        scales = np.ones(shape=(2, 4096, 4096), dtype='f4')
        domains[1, :, :] = limits
        scales[0, :, :] = -1
        off, scl = np.polynomial.polyutils.mapparms(domains, scales)
        rawIm = exposure.image.array
        normIm = off + scl*rawIm
        corr = np.zeros_like(rawIm)
        BAD = exposure.mask.getPlaneBitMask('BAD')

        # If defects were supplied we have a BAD pixel mask. Use it to
        # avoid linearizing known defects. But in any case add in and
        # avoid correcting the worst unmasked pixels.
        badMask = 0 != (exposure.mask.array & BAD)
        startingBADpixels = badMask.sum()

        # For the moment, do not linearize significantly negative
        # values. We currently mask most of these pixels externally,
        # but leaks are bad.
        assert len(exposure.getDetector()) == 1, "Fix me now we have multiple channels"
        amp = exposure.detector.getAmplifiers()[0]
        tooLow = rawIm < -amp.getReadNoise() * 10  # If permanent logic, add to config...
        badMask[tooLow] |= True
        corr[tooLow] = rawIm[tooLow]

        # Declare that all pixels above the available linearization
        # limits are SATurated.  When we get real data at MKO this
        # will actually be correct. Some lab flats from JHU mask too
        # many pixels, but we don't really care about the brightest
        # ones.  We set the values to the limits found on the
        # unlinearized flats, which is wrong but we cannot trust any
        # linearization. We need to replace the coefficients for these
        # pixels with more sane estimates before we can do anyting
        # better.
        saturated = rawIm > limits
        badMask[saturated] |= True
        corr[saturated] = limits[saturated]
        exposure.mask.array[saturated] |= (exposure.mask.getPlaneBitMask('SAT'))

        # Corrections are only valid for flux >= 0. The chebyshevs are
        # forced to 0 at flux=0, and are pretty low order, and the
        # negative values should be pretty small (for good pixels). We
        # should pretend that the correction is symmetric across 0, but
        # let it extrapolate instead. Hmm, CPL.
        goodMask = ~badMask
        corr[goodMask] = np.polynomial.chebyshev.chebval(normIm[goodMask],
                                                         coeffs[:, goodMask],
                                                         tensor=False)
        exposure.image.array[:] = corr
        exposure.mask.array[badMask] |= BAD

        endingBADpixels = badMask.sum()
        self.log.info(f'linearization added {endingBADpixels - startingBADpixels} BAD pixels')

        # variance is handled later...

        return exposure

    def applyNirLinearity(self, pfsRaw, exposure, linearity):
        """Apply NIR linearity corrections of some defined kind.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
            The raw exposure that we ran through ISR.
        exposure : `lsst.afw.image.Exposure`
            exposure to correct..
        linearity : `nirLinearity.NirLinearity`
            all parts of the linearity corrections

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            input exposure, corrected in place.
        """
        method = linearity.method

        if method == 'identity':
            return exposure
        elif method == 'np.polynomial.chebyshev':
            self.nirLinearityChebyshev(exposure, linearity)
            exposure.metadata.set('PFS ISR LINEARITY METHOD', method)
            # Need some SHA-like ID. And the filename.
        else:
            # blow up here?
            self.log.warn(f'ignoring unknown linearity method: {method}')
            return exposure

        return exposure

    def calcUTRWeights(self, nreads: int) -> np.ndarray:
        """Compute the weights for linear UTR, per Eq 4.24 in RHL

        Parameters
        ----------
        nreads : `int`
            The number of reads to use in the UTR calculation.

        Returns
        -------
        `np.ndarray`
            The weights to apply to the reads.
        """
        n = float(nreads)
        w = []
        for i in range(nreads):
            k = (12*i - 6*(n-1))/(n*(n*n - 1))
            w.append(k)
        return np.array(w)

    def calcUTRrates(self, cube: np.ndarray, nreads: Optional[int] = None) -> np.ndarray:
        """Apply UTR weights to the ramp reads to estimate a per-pixel arrival rate.

        Parameters
        ----------
        cube : `np.ndarray`
            The ramp cube to apply UTR to.
        nreads : `int`, optional
            The number of reads to use in the UTR calculation. default=all cube reads.

        Returns
        -------
        `np.ndarray`
            The estimated per-pixel arrival rate.
        """
        if nreads is None:
            n_i = len(cube)
        else:
            n_i = nreads
        rate_sum = np.zeros_like(cube[0])
        weights = self.calcUTRWeights(n_i)
        for i in range(n_i):
            k = weights[i]
            s1 = k*cube[i]
            rate_sum += s1

        return rate_sum

    def getDarkCube(self, nirDark, nreads: Optional[int] = None) -> np.ndarray:
        """Get the dark cube for the NIR ramp.

        Parameters
        ----------
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`
            The unprocessed to subtract.

        Returns
        -------
        `np.ndarray`
            The dark cube to subtract.
        """

        #
        cube = nirDark.getImageCube(nreads=nreads)
        # If gain was applied, back it out. We used to apply darks to the
        # 2-d image in e-, but have switch to applying it to the raw rampin ADU.
        gain = nirDark.metadata.get("GAIN", 1.0)
        if gain != 1.0:
            cube /= gain
        return cube

    def getDarkRead(self, nirDark, readNum) -> np.ndarray:
        """Get the dark read for the NIR ramp.

        This wraps the fact that we switched from post-gain darks to
        pre-gain darks.

        Parameters
        ----------
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`
            The unprocessed dark cube.
        readNum : `int`
            The read number to get.

        Returns
        -------
        `np.ndarray`
            The dark read.
        """

        dark = nirDark[readNum].array

        # If gain was applied, back it out.
        gain = nirDark.metadata.get("GAIN", 1.0)
        if gain != 1.0:
            dark /= gain
        return dark

    def makeNirExposure(self, pfsRaw, nirDark=None, doReturnRawCube=False):
        """Construct a 2D image from the NIR ramp data.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
            Provides access to the ramp that is to be
            run through ISR.
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`, optional
            Dark cube to subtract (if ``config.h4.useDarkCube`` is ``True``).
        doReturnRawCube : `bool`, optional
            If True, return the raw ramp cube as well as the 2D image.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            The 2-d image from the ramp.
        """

        if self.config.h4.quickCDS:
            self.log.info("creating quick CDS.")
            nirImage = self.makeCDS(pfsRaw)
            flux = None
        else:
            # Too many interacting switches for CDS/UTR/darks. Especially note 2-d doDark is still possible.
            self.log.info("reading full ramp...")
            deltas = self.makeUTRdeltas(pfsRaw)

            # We do not currently use the indices of the masked pixels,
            # but do get them for when we put in the effort.
            deltas, posIdx, negIdx = self.correctCRs(pfsRaw, deltas)

            # Switch to accumulated flux for dark subtraction and
            # UTR weighting. This is actually pretty expensive; should
            # rethink.
            # But do at least save memory by doing this in place.
            flux = np.cumsum(deltas, axis=0, out=deltas)
            del deltas  # Get rid of the name to avoid stupidities

            if self.config.h4.applyUTRWeights:
                if nirDark is not None:
                    self.log.info("subtracting dark cube.")
                    darkCube = self.getDarkCube(nirDark, len(flux))
                    flux -= darkCube

                self.log.info("applying UTR weights.")
                rates = self.calcUTRrates(flux)
                nirImage = rates * len(flux)
            else:
                # CDS, basically.
                if nirDark is not None:
                    nirImage = ((flux[-1] - self.getDarkRead(nirDark, len(flux)-1)) -
                                (flux[0] - self.getDarkRead(nirDark, 0)))
                else:
                    nirImage = flux[-1] - flux[0]

        exposure = self._makeExposure(pfsRaw, nirImage)
        return exposure, flux if doReturnRawCube else None

    def makeRawDataArray(self, pfsRaw, readNum, fromArray=None) -> np.ndarray:
        """Fetch a single data read."""

        if fromArray is not None:
            return fromArray[readNum]
        else:
            return pfsRaw.getRawDataImage(readNum).array

    def makeRawIrpArray(self, pfsRaw, readNum, forceIrp1=True, fromArray=None) -> np.ndarray:
        """Construct a full IRP1 reference image"""

        if fromArray is not None:
            rawIrpArray = fromArray[readNum]
        else:
            rawIrpArray = pfsRaw.getRawIrpImage(readNum).array
        if forceIrp1:
            im = self.constructFullIrp(pfsRaw, rawIrpArray)
        else:
            im = rawIrpArray

        return im

    def makeSimpleStack(
        self,
        pfsRaw,
        reader,
        r0: int = 0,
        r1: int = -1,
        nreads: Optional[int] = None,
        bbox: Optional["geom.Box2I"] = None,
        doDeltas: bool = True,
    ) -> np.ndarray:
        """Return all the raw frames in a single 3d stack.

        Args
        ----
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        reader : `callable`
          A function to read one image.
        r0 : `int`
          The 0-indexed read to start from.
        r1 : `int`
          The 0-indexed read to end with. Inclusive.
        nreads : `int`
          The number of reads to return, evenly spaced between r0 and r1.
          Note that this uses `np.linspace`, so you want to be alert.
        bbox : `Box2I`, optional
          A bounding box to apply to each image.
        doDeltas : `bool`
          Whether to return the ramp as per-read incremental changes, or per-read total flux.

        Returns
        -------
        stack : 3-d float32 numpy array
           the ramp, with axis 0 being the reads.
        """

        r0 = pfsRaw.positiveIndex(r0)
        r1 = pfsRaw.positiveIndex(r1)
        if r1 <= r0:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}):')
        if nreads is None:
            nreads = r1 - r0 + 1
        reads = np.linspace(r0, r1, nreads, dtype='i2')

        for r_i, read1 in enumerate(reads):
            readImg = reader(pfsRaw, read1)
            if bbox is not None:
                readImg = readImg[bbox]
            if r_i == 0:
                h, w = readImg.shape
                stack = np.empty(shape=(nreads, h, w), dtype='f4')
            stack[r_i, :, :] = readImg

        if doDeltas:
            stack = np.diff(stack, axis=0,
                            prepend=np.zeros_like(stack[0:1]))
        return stack

    def makeRawDataDeltas(self, pfsRaw,
                          r0=0, r1=-1, nreads=None, bbox=None) -> np.ndarray:
        """Return the raw image reads in a single 3d stack.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        r0 : `int`
          The 0-indexed read to start from.
        r1 : `int`
          The 0-indexed read to end with. Inclusive.
        nreads : `int`
          The number of reads to return, evenly spaced between r0 and r1.
          Note that this uses `np.linspace`, so you want to be alert.
        bbox : `Box2I`, optional
          A bounding box to apply to each image.
        """
        return self.makeSimpleStack(pfsRaw, reader=self.makeRawDataArray,
                                    r0=r0, r1=r1, nreads=nreads, bbox=bbox)

    def makeRawIrpDeltas(self, pfsRaw,
                         r0=0, r1=-1, nreads=None, bbox=None) -> np.ndarray:
        """Return the raw IRP reads in a single 3d stack.

        See .makeRawDataDeltas.

        Note that IRPn reads are expanded to full IRP1.
        """
        return self.makeSimpleStack(pfsRaw, reader=self.makeRawIrpArray,
                                    r0=r0, r1=r1, nreads=nreads, bbox=bbox)

    def makeRawDataCube(self, pfsRaw,
                        r0=0, r1=-1, nreads=None, bbox=None) -> np.ndarray:
        """Return the raw data reads in a single 3d stack, but not incremental.

        This is used just for the ASIC glitch correction.
        """
        return self.makeSimpleStack(pfsRaw, reader=self.makeRawDataArray,
                                    r0=r0, r1=r1, nreads=nreads, bbox=bbox,
                                    doDeltas=False)

    def makeRawIrpNcube(self, pfsRaw,
                        r0=0, r1=-1, nreads=None, bbox=None) -> np.ndarray:
        """Return the raw IRP reads in a single 3d stack, but not incremental.

        This is used just for the ASIC glitch correction.
        """
        from functools import partial
        reader = partial(self.makeRawIrpArray, forceIrp1=True)
        return self.makeSimpleStack(pfsRaw, reader=reader,
                                    r0=r0, r1=r1, nreads=nreads, bbox=bbox,
                                    doDeltas=False)

    def repairWithWindow(self, deltas, badPixels,
                         corrRad, corrMin, corrMax,
                         repairMask=None,
                         otherIgnore=None) -> np.ndarray:
        """Replace single pixels with a local correction

        Parameters
        ----------
        deltas : `np.ndarray`
           The 3-d cube of the ramp
        badPixels : `np.ndarray`
           A 2-d array of the indices of the reads to replace
        corrRad : `int`
           The width of the window (along reads) to use for replacement values
        corrMin : `int`
           In the deltas cube, the first read we can use for replacements
        corrMax : `int`
           In the deltas cube, the last read we can use for replacements
        repairMask : `np.ndarray`, optional
           If set, only apply corrections to these pixels.
        otherIgnore : `np.ndarray`, optional
           If set, array of read indices we cannot use for replacements.

        Returns
        -------
        correctionWindow : `np.ndarray`
           The windows we use to make replacements. NaNs ignored.
           For engineering. Drop once we believe or replace this.
        """
        # Use neighboring reads, clipped by the ends of the ramp
        lowIdx0 = badPixels - corrRad
        lowIdx = np.maximum(lowIdx0, np.zeros_like(badPixels)+corrMin)

        # If our window hits the end of the valid window, move the
        # start of the window down.
        highIdx0 = badPixels + corrRad
        highIdx = np.minimum(highIdx0, np.zeros_like(badPixels)+corrMax-1)
        highClip = np.where(highIdx != highIdx0)
        lowIdx[highClip] = highIdx[highClip] - 2*corrRad+1

        # Average the nearby reads, ignoring any we are correcting.
        correctionWindow = np.zeros(shape=(2*corrRad+1, deltas.shape[1], deltas.shape[2]))
        for i in range(2*corrRad + 1):
            idx_i = lowIdx + i
            corr1 = np.take_along_axis(deltas, idx_i, axis=0)
            corr1[idx_i == badPixels] = np.nan
            if otherIgnore is not None:
                corr1[idx_i == otherIgnore] = np.nan
            correctionWindow[i, :, :] = corr1
        corrections = np.nanmean(correctionWindow, axis=0)

        # This is not the cheapest way to handle repairMask, but re-writing to
        # make that efficient would be error-prone work.
        if repairMask is not None:
            badPixels[~repairMask] = 0
            corrections[~repairMask] = deltas[0][~repairMask]
        np.put_along_axis(deltas, badPixels, corrections, axis=0)
        return correctionWindow

    def correctCRs(self, pfsRaw, deltas: np.ndarray):
        """Correct cosmic ray hits using the array of flux deltas.

        "Correct" means to detect both the maximum and the minimum increments for each
        pixel, and replace both of those by a local estimate of the rate.

        Since we are already paying the cost of detecting hot pixels, fold correction
        of ASIC-injected bad pixels in here. Currently only n3, channel 24. In the
        delta ramps these appear as paired high/low pixels. We do need to calculate a
        robust sigma, which is expensive enough that we want to skip it if possible.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        deltas : `np.ndarray`
           Corrected in place.

        Returns
        -------
        deltas : `np.array`
           The corrected ramp
        posIdx, negIdx : `np.array`
           The per-pixel indexes into deltas of the corrected reads.
        """

        if not self.config.h4.doCR:
            return deltas, None, None

        # Oh, ultra gross. Just realized that for normal exposures the shutter is
        # always closed for the first read and at least the last two reads, so
        # we do not want to use those for interpolation.
        # Worse, Ar and Ne exposures are illuminated for just one read, so *all*
        # illuminated pixels would be rejected! If we start taking shorter quartz
        # exposures those will lose significant flux.
        # Finally, we do not want to use the closed-shutter reads to identify the
        # minimum pixels: arc lines, etc. always get dinged. This is pretty dangerous
        # for short exposures.
        #
        # Hack that following logic in for now, but need to somehow fix correctly -- CPL.
        #  - if nread is "short", apply no correction: must use CR splits for Ar/Ne
        #  - if exptype is not DARK, do not use first or two last reads
        #
        corrRad, corrMin, corrMax = self.calcCorrectionWindow(pfsRaw, len(deltas))
        if corrMax - corrMin < self.config.h4.crMinReads:
            self.log.warn(f'ramp is too short to correct CRs (need {self.config.h4.crMinReads}, '
                          f'have {corrMax}-{corrMin}): you must use splits!')
            return deltas, None, None

        # Identify both min and max reads for each pixel. This
        # is fairly expensive...
        posIdx = np.argmax(deltas, axis=0, keepdims=True)
        negIdx = np.argmin(deltas, axis=0, keepdims=True)

        if False:  # Instead of repairing pixels, mask them for later processing.
            np.put_along_axis(deltas, posIdx, axis=0, values=np.nan)
            np.put_along_axis(deltas, negIdx, axis=0, values=np.nan)
        else:
            self.repairWithWindow(deltas, posIdx,
                                  corrRad, corrMin, corrMax,
                                  otherIgnore=negIdx)

            # Need to correct for masking every max read. One way is to mask every
            # min read. I hate this.
            if True:
                self.repairWithWindow(deltas, negIdx,
                                      corrRad, corrMin, corrMax,
                                      otherIgnore=posIdx)
            else:
                negIdx *= 0
                self.log.warn('NOT accounting for min pixels in CR step.')

        # we mask/replace two reads for every single pixel. Let the caller decide
        # what to do about that, but return the indices of both.
        return deltas, posIdx, negIdx

    def loadBadIRPpixels(self, detectorName: str) -> np.ndarray:
        """Return the bad IRP row pixel list."""

        filename = 'badRefPixels.yaml'
        absFilename = os.path.join(getPackageDir("drp_pfs_data"), "h4", filename)
        if not os.path.exists(absFilename):
            self.log.warn(f'{absFilename} not found')
            return np.zeros(dtype=np.int16, shape=())

        with open(absFilename) as f:
            cfg = yaml.YAML(typ="safe", pure=True).load(f)

        if detectorName not in cfg:
            self.log.warn(f'{detectorName} not defined in {absFilename}')
            return np.zeros(dtype=np.int16, shape=())

        return np.array(cfg[detectorName])

    def correctBadIrpPixels(self, pfsRaw, refImg):
        """Given an IRP1 image, interpolate over known bad IRP row pixels.

        Given the large offsets between the reference pixels and the drifts across the full reads,
        it only really makes sense to run this on two subtracted IRP reads.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        refImg : `numpy.ndarray`
           An IRP1 image. Corrected in place.
        """
        if not self.config.h4.doIRPbadPixels:
            return
        if pfsRaw.irpN != 1:
            warnings.warn('Do not know how to correct for bad IRP pixels in non-IRP1 images')
            return

        badPixels = self.loadBadIRPpixels(pfsRaw.detector.getName())

        # We really want to be doing this on IRP diffs, given the huge common per-pixel offsets.
        # If not, need to somehow flatten out those offsets.
        #
        for i in range(32):
            pixLow = i*128
            pixHigh = (i+1)*128
            chan0 = refImg[pixLow:pixHigh, :]

            chanBadPix = badPixels[(badPixels >= pixLow) & (badPixels < pixHigh)] - pixLow
            if len(chanBadPix) > 0:
                # Use median of entire column in channel: trying to correct for
                # fairly local (as seen in column) 1/f
                chan0[chanBadPix, :] = np.nanmedian(chan0, axis=0)

    def loadBadAsicChannels(self, detectorName: str) -> tuple:
        # Yes, yes, get this from a config file...
        badAsicChannels = dict(n3=(24,))
        return badAsicChannels.get(detectorName, ())

    def getSimpleDiffIrp(self, pfsRaw, rawDiffIrp) -> np.ndarray:
        """Return a diff IRP image which is just the median across the channel columns"""
        nchan = pfsRaw.nchan
        h, w = rawDiffIrp.shape
        chan_w = h//nchan

        out = np.zeros_like(rawDiffIrp)
        for chan_i in range(nchan):
            rowLow = chan_i*chan_w
            rowHigh = (chan_i + 1)*chan_w
            chan0 = rawDiffIrp[rowLow:rowHigh, :]

            chanVec = np.nanmedian(chan0, axis=0, keepdims=True)
            out[rowLow:rowHigh, :] = np.tile(chanVec, (h, 1))
        return out

    def getFinalDiffIrp(self, pfsRaw, rawDiffIrp, useFft=True) -> np.ndarray:
        """Given a simple IRP difference image, return the final IRP image.

        We do a few things here:
         - replace the detector columns from known bad IRP row pixels with
           the median of their channel
         - per-channel, smooth the IRP image rows over {filterWidth} pixels

        For the moment, we assume that the following has been done before
        we are called:
          - IRPn images have been filled out into IRP1 images.
          - Any crosstalk/leakage from the data images has been removed.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        rawDiffIrp : `numpy.ndarray`
            The raw difference IRP image. Full-size (i.e. IRP1)
        filterWidth : `int`
            The width of the smoothing filter to use. Should be odd.
            If 0, do not apply any filter: simply return the input image
            If -1, return an image where each channel is replaced by its per-col median.

        Returns
        -------
        finalIrp : `np.ndarray`
            The processed IRP image.
           """

        self.correctBadIrpPixels(pfsRaw, rawDiffIrp)

        skipChannels = self.loadBadAsicChannels(pfsRaw.detector.getName())

        filterWidth = self.config.h4.IRPfilter
        if filterWidth == 0:
            return rawDiffIrp
        if filterWidth == -1:
            return self.getSimpleDiffIrp(pfsRaw, rawDiffIrp)
        if filterWidth % 2 == 0:
            raise ValueError(f'filterWidth must be odd: {filterWidth}')
        if skipChannels is None:
            skipChannels = []

        nchan = pfsRaw.nchan
        h, w = rawDiffIrp.shape
        chan_w = h//nchan

        # Construct a cos bell filter to run over each channel.
        # Pad filter with zeros to same width as padded channel.
        ww_half = filterWidth//2
        padWidth = ww_half*2
        filtCore = scipy.signal.windows.hann(filterWidth)
        filtCore /= filtCore.sum()
        filt1 = np.zeros(shape=(chan_w+padWidth, 1), dtype='f4')
        fc = (chan_w+filterWidth)//2
        filt1[fc-ww_half:fc+ww_half+1, 0] = filtCore
        filt = np.tile(filt1, (1, w))

        # scratch space for the padded channel
        chan1 = np.zeros(shape=(chan_w+padWidth, w), dtype=rawDiffIrp.dtype)
        out = np.zeros_like(rawDiffIrp)
        for chan_i in range(nchan):
            rowLow = chan_i*chan_w
            rowHigh = (chan_i + 1)*chan_w
            chan0 = rawDiffIrp[rowLow:rowHigh, :].copy()

            if chan_i in skipChannels:
                out[rowLow:rowHigh, :] = chan0
                self.log.debug(f'skipping repairs on channel {chan_i}')
            else:
                # Pad channel with median of channel out to half the filter width
                chan1[ww_half:chan_w+ww_half, :] = chan0
                chan1[:ww_half, :] = np.nanmedian(chan0[:ww_half, :], axis=0, keepdims=True)
                chan1[-ww_half:, :] = np.nanmedian(chan0[-ww_half:, :], axis=0, keepdims=True)

                if useFft:
                    chan_f = scipy.signal.fftconvolve(chan1, filt, mode='same', axes=0)
                else:
                    chan_f = scipy.signal.oaconvolve(chan1, filt, mode='same', axes=0)
                out[rowLow:rowHigh, :] = chan_f[ww_half:chan_w+ww_half, :]

        return out

    def interpolateChannelIrp(self, pfsRaw,
                              rawChan: np.ndarray,
                              doFlip: bool) -> np.ndarray:
        """Given an IRPn channel from a PFSB file, return IRP1-sized channel image.

        Parameters
        ----------
        rawChan : array
           The raw IRP channel, with the columns possibly subsampled by a factor of self.irpN.
        doFlip : `bool`
           Whether we need to treat this channel as read out B-to-T (R-to-L in H4). Only
           meaningful if we are filtering in time.

        Returns
        -------
        im : `np.ndarray`
           the full-sized interpolated reference pixel channel.

        We do not yet know how best to interpolate, so simply repeat the pixels refRatio times.
        We do not handle the IRP row bad pixel mask, which is a *major* flaw.
        """

        irpN = pfsRaw.irpN
        if irpN <= 1:
            return rawChan

        # Allow the future possibility of using readout order.
        temporalFilter = False

        irpHeight, irpWidth = rawChan.shape
        refChan = np.empty(shape=(irpHeight * irpN, irpWidth), dtype=rawChan.dtype)

        if doFlip and temporalFilter:
            rawChan = rawChan[::-1, :]

        # For now, simply repeat reference pixels
        for i in range(0, irpN):
            refChan[i::irpN, :] = rawChan

        if doFlip and temporalFilter:
            refChan = refChan[::-1, :]

        return refChan

    def constructFullIrp(self, pfsRaw, refImg) -> np.ndarray:
        """Given an IRPn image, return IRP1 image.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        refImg : ndarray
          A raw reference read from the ASIC.

        Returns
        -------
        img : `np.ndarray`
          full 4096x4096 image.

        Notes
        -----
        - The detector was read out in nChannel channels, usually 32, but possibly 16, 4, or 1.

        - By default, the read order from the detector of pairs of channels is --> <--, but
          it can any other pairing.
          The ASIC "corrects" that order so that the image always "looks" right: the
          column-order of the image from the ASIC is spatially, not temporally, correct.

        - the N:1 ratio of science to reference pixels is deduced from the size of the image.

        - self.irpOffset tells us the position of the reference pixel within the N
          science pixels. It must be >= 1 (there must be at least one
          science pixel before the reference pixel). The ASIC default is
          for it to be the last pixel in the group, but we usually try to
          put the reference pixel in the middle of the block of science pixels.
        """

        height, width = refImg.shape
        nchan = pfsRaw.nchan

        # If we are a full frame, no interpolation is necessary.
        if height == pfsRaw.h4Size:
            return refImg

        refChanHeight = height // nchan

        refChans = []
        readOrders = pfsRaw.getH4channelReadOrder()
        self.log.debug(f'filling out IRP{pfsRaw.h4Size // height} image, with {readOrders=}')
        for c_i in range(nchan):
            rawChan = refImg[c_i*refChanHeight:(c_i+1)*refChanHeight, :]
            doFlip = readOrders[c_i%2]

            # This is where we would intelligently interpolate.
            refChan = self.interpolateChannelIrp(pfsRaw, rawChan, doFlip)
            refChans.append(refChan)

        fullRefImg = np.vstack(refChans)

        return fullRefImg

    def applyIRPcrosstalk(self, pfsRaw, ref, image) -> np.ndarray:
        """Correct crosstalk from the data pixels to the IRP

        It appears that crosstalk varies "per channel", but we have
        not looked into what happens with non-32channel ramps.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        ref : `np.ndarray`
           The IRP image to modify.
           Corrected in place.
        image : `np.ndarray`
           The crosstalk source image.

        Returns
        -------
        ref : `np.ndarray`
           The corrected image. Same as `ref`.

        """
        if not self.config.h4.doIRPcrosstalk:
            return ref
        raise NotImplementedError('not loading crosstalk yet')

        crosstalkFactors = self.loadIRPcrosstalk(pfsRaw)

        nchan = pfsRaw.nchan
        chanStep = 32//nchan
        chanWidth = pfsRaw.h4Size / nchan

        for c in range(nchan):
            xslice = slice(c*chanWidth, (c+1)*chanWidth)
            ref[xslice, :] -= image[xslice, :] * crosstalkFactors[c*chanStep]

        return ref

    def borderCorrect(self, pfsRaw, image, colWindow=4,
                      doRows=True, doCols=True):
        """This is the "standard" Teledyne 'refPixel4' border reference pixel scheme.

        Step 1:
           For each channel, average all 8 top&bottom rows to one number.
           Subtract that from the channel.

        Step 2:
            Take a 9-row running average of the left&right columns.
            Subtract that from each row.

        The output is not lovely but is functional. We are mainly here just to duplicate their
        standard logic. So please do *not* modify/improve this particular routine: we want to have
        it for any support discussions. Ok, maybe this should have been named 'borderCorrectRefPixel4'.
        """
        refPix = 4

        imHeight, imWidth = image.shape
        if imHeight != pfsRaw.h4Size or imWidth != pfsRaw.h4Size:
            raise RuntimeError('not ready to deal with non-full images')

        if doCols:
            top = image[imHeight-refPix:, :]
            bottom = image[:refPix, :]

            chanWidth = imWidth // pfsRaw.nchan
            for c_i in range(pfsRaw.nchan):
                xlow = c_i * chanWidth
                xhigh = xlow + chanWidth
                # Leave ref columns unmolested.
                if c_i == 0:
                    xlow = refPix
                elif c_i == pfsRaw.nchan-1:
                    xhigh = imWidth-refPix
                ichan = image[:, xlow:xhigh]

                chanOffset = (top[:, xlow:xhigh].mean()
                              + bottom[:, xlow:xhigh].mean()) / 2
                ichan -= chanOffset

        if doRows:
            sideRefImage = np.ndarray((imHeight, refPix*2), dtype=image.dtype)
            sideRefImage[:, :refPix] = image[:, :refPix]
            sideRefImage[:, refPix:] = image[:, imWidth-refPix:]
            sideCorr = np.zeros((imHeight, 1))

            # Not quite right at the ends -- should not be including the reference rows.
            for row_i in range(colWindow, imHeight-colWindow+1):
                sideCorr[row_i] = sideRefImage[row_i-colWindow:row_i+colWindow, :].mean()
                image[row_i, :] -= sideCorr[row_i]

        return image

    def makeCDS(self, pfsRaw, r0=0, r1=-1):
        """Get a single reference-corrected CDS from a ramp. Does not read full ramp.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        r0 : `int`
            The read to start from. 0-indexed.
        r1 : `int`
            The read to end at. 0-indexed.

        Returns
        -------
        image : `afw.image.ImageF`
            reference-corrected image.
        """

        r0 = pfsRaw.positiveIndex(r0)
        r1 = pfsRaw.positiveIndex(r1)
        if r1-r0 < 1:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}) + 1:')

        image0 = self.makeRawDataArray(pfsRaw, r0)
        image1 = self.makeRawDataArray(pfsRaw, r1)

        if self.config.h4.useIRP and pfsRaw.irpN > 0:
            ref0 = self.makeRawIrpArray(pfsRaw, r0)
            self.applyIRPcrosstalk(pfsRaw, ref0, image0)
            ref1 = self.makeRawIrpArray(pfsRaw, r1)
            self.applyIRPcrosstalk(pfsRaw, ref1, image1)

            dref = ref1 - ref0
            dref = self.getFinalDiffIrp(pfsRaw, dref)

            image1 -= image0
            image1 -= dref
        else:
            image0 = self.borderCorrect(pfsRaw, image0)
            image1 = self.borderCorrect(pfsRaw, image1)
            image1 -= image0

        return image1

    def calcCorrectionWindow(self, pfsRaw, nread, corrRad=2):
        """Figure out the bounds of any interpolation/correction window for this ramp

        If illuminated, we cannot use the reads when the shutter is closed;
        currently (2025-01) the first and the last two.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        nread : `int`
          The number of actual reads we have kept.
        corrRad : `int`
          How many reads on either side we want to use

        Returns
        -------
        corrRad : `int`
           The radius of the window we keep.
        corrMin, corrMax : `int`
           The first and last reads we can use
        """

        exptype = pfsRaw.visitInfo.observationType

        if exptype == 'dark' or nread < pfsRaw.getNumReads():
            corrMin = 0
            corrMax = nread-1
        else:
            # TODO: figure the shutter-closed reads out from the header. Need to add cards for that.
            # This is correct for 2024-ish.
            corrMin = 1
            corrMax = nread-3

        return corrRad, corrMin, corrMax

    def repairAsicSpikes(self, pfsRaw, cube: np.ndarray, channel: int, sigClip=None, doTest=False):
        """Replace bad pixels from bad ASIC channels

        The bad pixels are positive going spikes from the ASIC, which can hit
        both data and IRP pixels. We correct the pixels early, on the separate
        data and IRP cubes, before even being converted to IRP1 or incremental UTR.

        We expect to only run this for a few (1?) channel.

        Need to iterate until no more spikes are found.....
        Need to handle spikes at the ends of the ramp.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        cube : `np.ndarray`
           The 3-d cube of the ramp. Corrected in place.
        channel : `int`
           The index of the (32-channel read) channel we are correcting
        """

        if not self.config.h4.repairAsicSpikes:
            return

        if pfsRaw.nchan != 32:
            self.log.warn(f'nchan ({pfsRaw.nchan}) != 32, so not (yet) correcting ASIC glitches')
            return

        if sigClip is None:
            sigClip = self.config.h4.repairAsicSpikesSigma
        badChan = cube[:, channel*128:(channel+1)*128, :]

        diffChan = np.diff(badChan, axis=0, prepend=np.zeros_like(badChan[0:1]))
        p25, meds, p75 = np.percentile(diffChan, [25, 50, 75], axis=0)
        iqrSig = 0.741*(p75 - p25)

        # the reads we are correcting have a single hot pixel, just like CRs. But the flux
        # levels return on the following read, unlike CRs. So, in delta stacks look for
        # positive spikes just before a matching negative spike. Or in cumulative stacks,
        # look for a positive spike between two "normal" reads
        brightMask = np.abs(diffChan) > (meds + sigClip*iqrSig)
        bright_w = np.where(brightMask)

        # Todo: Handle spikes at the ends of the ramp separately -- CPL
        atEnd = (bright_w[0] >= len(diffChan)-1 | (bright_w[0] == 0))
        bright_w1 = bright_w[0][~atEnd], bright_w[1][~atEnd], bright_w[2][~atEnd]
        next_w1 = (bright_w1[0]+1), bright_w1[1], bright_w1[2]

        peaks1 = diffChan[bright_w1]
        next1 = diffChan[next_w1]
        diffs = peaks1 + next1

        meds1 = meds[bright_w1[1:]]
        iqrSig1 = iqrSig[bright_w1[1:]]

        fix_i = np.abs(diffs) < (meds1 + sigClip*iqrSig1)
        fix_w = tuple([w1[fix_i] for w1 in bright_w1])
        if True:  # Use the average of the two neighboring reads
            fix_wm1 = (fix_w[0]-1, fix_w[1], fix_w[2])
            fix_wp1 = (fix_w[0]+1, fix_w[1], fix_w[2])
            repairVals = (badChan[fix_wm1] + badChan[fix_wp1]) / 2
        else:     # Use the median of the entire ramp
            repairVals = meds1[fix_w[1:]]
        badChan[fix_w] = repairVals

        # Return for testing
        if doTest:
            return fix_w, badChan, diffChan, bright_w1, iqrSig, meds, repairVals
        else:
            return fix_w

    def makeUTRdeltas(self, pfsRaw, r0=0, r1=-1, nreads=None, bbox=None,
                      showTimes=False) -> np.ndarray:
        """Return all the fully IRP-corrected frames in a single 3d stack.

        Given two raw data images d0 and d1, and two raw IRP images i0 and i1, the net CDS image
        can be either (d1 - i1) - (d0 - i0), or (d1 - d0) - (i1 - i0). The IRP row has various
        artifacts which make using the latter "nicer", or at least easier to make sense of. So that
        is the way we do it.
        In particular there are:
          - bad IRP row pixels
          - fixed pixel-to-pixel offsets
          - bad ASIC channels

        Note that there is also up to ~1% printthrough from data->IRP, and presumably from IRP->data

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
          Access to the raw ramp
        r0 : `int`
          The 0-indexed read to start from.
        r1 : `int`
          The 0-indexed read to end with. Inclusive.
        nreads : `int`, optional.
          The number of reads to return, spaced evenly between r0 and r1.
        filterIrp : `bool`
            Whether to apply known corrections to the IRP images.
        bbox : `lsst.geom.Box2I`
            The region of the image to return. If None, return the whole image.

        Returns
        -------
        deltas : 3-d float32 numpy array
           the UTR stack, with axis 0 being the reads up the ramp. We return
           the flux increments between reads, not the actual fluxes.
        """

        r0 = pfsRaw.positiveIndex(r0)
        r1 = pfsRaw.positiveIndex(r1)
        if r1-r0 < 1:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}) + 1:')
        if nreads is None:
            nreads = r1 - r0 + 1
        reads = np.linspace(r0, r1, nreads, dtype='i2')

        # n3/18321 channel 24 is covered with pixels which have spikes from the ASIC.
        # We need to repair both the data and IRP pixels separately, and furthermore
        # want to repair the non-interpolated IRP images. So read both in first and
        # correct them before building the proper deltas. Since I/O is significant,
        # we keep the raw data and IRP cubes, then re-read from those instead of
        # from disk.
        #
        badChans = self.loadBadAsicChannels(pfsRaw.detector.getName())
        if badChans is not None:
            rawData = self.makeRawDataCube(pfsRaw=pfsRaw, r0=r0, r1=r1, nreads=nreads)
            rawIrps = self.makeRawIrpNcube(pfsRaw=pfsRaw, r0=r0, r1=r1, nreads=nreads)

            for badChan in badChans:
                fixedData_w = self.repairAsicSpikes(pfsRaw, rawData, badChan)
                fixedIrp_w = self.repairAsicSpikes(pfsRaw, rawIrps, badChan)
                self.log.info(f'repaired {len(fixedData_w[0])} data and {len(fixedIrp_w[0])} IRP '
                              f'ASIC pixels in channel {badChan}')
        else:
            rawData = rawIrps = None

        # Grab the components of read 0, which we will subtract from all the others.
        if rawData is not None:
            data0 = self.makeRawDataArray(pfsRaw, 0, fromArray=rawData)
            irp0 = self.makeRawIrpArray(pfsRaw, 0, fromArray=rawIrps)
        else:
            data0 = self.makeRawDataArray(pfsRaw, r0, fromArray=rawData)
            irp0 = self.makeRawIrpArray(pfsRaw, r0, fromArray=rawIrps)
        self.applyIRPcrosstalk(pfsRaw, irp0, data0)

        # We are not squirreling away the bbox, but really should for the final Exposure
        if bbox is None:
            stackShape = (nreads-1, *data0.shape)
        else:
            stackShape = (nreads-1, bbox.getHeight(), bbox.getWidth())
        stack = np.empty(shape=stackShape, dtype='f4')
        for r_idx, r_i in enumerate(reads):
            if r_idx == 0:
                continue
            t0 = time.time()
            if rawData is not None:
                data1 = self.makeRawDataArray(pfsRaw, r_idx, fromArray=rawData)
                irp1 = self.makeRawIrpArray(pfsRaw, r_idx, fromArray=rawIrps)
            else:
                data1 = self.makeRawDataArray(pfsRaw, r_i, fromArray=rawData)
                irp1 = self.makeRawIrpArray(pfsRaw, r_i, fromArray=rawIrps)
            t1 = time.time()
            self.applyIRPcrosstalk(pfsRaw, irp1, data1)
            dirp = irp1 - irp0
            ddata = data1 - data0
            ddata -= self.getFinalDiffIrp(pfsRaw, dirp)
            if bbox is None:
                stack[r_idx-1, :, :] = ddata
            else:
                stack[r_idx-1, :, :] = ddata[bbox.getBeginY():bbox.getEndY(),
                                             bbox.getBeginX():bbox.getEndX()]
            t2 = time.time()
            if showTimes:
                print(f'cds {r_i} io1={t1-t0:0.3f} proc={t2-t1:0.3f}')

        # Warning: prepend=0 causes promotion to float64
        # This should arguably be constructed on-the-fly, read by read.
        return np.diff(stack, axis=0, prepend=np.zeros_like(stack[0:1]))

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
        filename = self.config.h4.ipc[detectorName]
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

        ipcCoeffs: dict[int, dict[int, float]] = {}
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

    def nirDarkCorrection(self, pfsRaw: "PfsRaw", exposure: "ExposureF", dark: "ImageCube") -> None:
        """Apply NIR dark correction to the exposure

        The NIR dark correction is applied to the exposure in place.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
            Access to the raw ramp
        exposure : `lsst.afw.image.ExposureF`
            Exposure to correct.
        dark : `lsst.obs.pfs.dark.ImageCube`
            Dark correction to apply.
        """

        nObjRead = pfsRaw.getNumReads()
        nDarkRead = dark.getNumReads()
        if nObjRead > nDarkRead:
            darkScale = nObjRead/nDarkRead
            self.log.warn(f"More reads in object ({nObjRead}) than dark ({nDarkRead}); "
                          f"scaling dark by {darkScale:0.3f}")
        else:
            darkScale = 1.0
            nDarkRead = nObjRead
        try:
            darkImage = dark[nDarkRead-1]
        except Exception as e:
            self.log.warn(f"No dark image available for NIR dark correction of read {nDarkRead}: {e}")
            return

        darkArray = darkImage.array.copy()
        darkArray *= darkScale

        # This should be pulled out into some NirDark class
        # The dark imageCube itself has no variance, so we need to add it
        # Readnoise needs to be worked out: the dark cube is from a stack of darks, and we use UTR to
        #    get to our selected read. But as of 2025-03 there is observably more read noise than expected.
        readNoise = dark.metadata["READNOISE"]
        varArray = darkArray.copy()  # assumes photon noise
        varArray += (readNoise/np.sqrt(nDarkRead))**2
        darkImage = afwImage.makeMaskedImageFromArrays(darkArray, None, varArray)

        exposure.maskedImage -= darkImage

    def _getMetadataName(self):
        return None                     # don't write metadata; requires fix to ip_isr
