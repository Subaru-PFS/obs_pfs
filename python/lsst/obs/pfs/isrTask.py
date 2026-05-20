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
from types import SimpleNamespace

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
from . import h4Linearity
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

    quickCDS = pexConfig.Field(
        dtype=bool, default=None, optional=True,
        doc="Only consider last and first reads instead of the full ramp cube. "
            "None (default) = dispatch by observation type (CDS for arcs, UTR otherwise).")
    useIRP = pexConfig.Field(dtype=bool, default=True,
                             doc="Use Interleaved Reference Pixel planes if available")
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

    doCR = pexConfig.Field(
        dtype=bool, default=True,
        doc="Run the iterative UTR-rate-based CR + ASIC-glitch detector on "
            "the linearized cube after dark subtraction. Replaces all prior "
            "H4 CR-correction code paths.",
    )
    rateCRnSigma = pexConfig.Field(
        dtype=float, default=5.0,
        doc="Per-pixel sigma threshold for the iterative CR/glitch detector.",
    )
    rateCRsigmaFloorADU = pexConfig.Field(
        dtype=float, default=8.0,
        doc="Minimum sigma in ADU/read; protects faint pixels from collapsing the threshold "
            "to MAD-noise.",
    )
    doDeglitch = pexConfig.Field(
        dtype=bool, default=True,
        doc="Detect ASIC glitches alongside CRs in the iterative detector. "
            "Default True: run glitch detection on all pixels, all channels. "
            "Detection is needed even when glitches are not corrected — a "
            "glitch up-spike must be recognized as part of a pair so it is "
            "not misclassified as a CR. Set False to disable glitch detection "
            "entirely; CRs are still detected. Whether detected glitches are "
            "*corrected* is controlled separately by ``correctGlitches``. "
            "Only meaningful when ``doCR`` is True.",
    )
    deglitchAmplitudeMinADU = pexConfig.Field(
        dtype=float, default=0.0,
        doc="Minimum residual amplitude (in ADU) required for an ASIC-glitch "
            "pair to be classified. 0 (default) means the only floor is the "
            "CR threshold (``nSigma * sigma``). Set higher (e.g. 100) to "
            "suppress faint-end deglitching where the classifier is less "
            "reliable; bright glitches above this floor are still picked up. "
            "Implemented as: at least one of the two deltas in the pair must "
            "have |residual| above this value.",
    )
    correctGlitches = pexConfig.Field(
        dtype=bool, default=False,
        doc="Correct interior ASIC-glitch pairs: subtract them from the "
            "linearized cube and exclude their deltas from the UTR rate. "
            "Default False — interior pairs are still detected (``ASIC_GLITCH`` "
            "stamped, up-spike not misclassified as a CR) but left in place; "
            "a real symmetric +A/-A pair cancels on its own in the mean rate, "
            "so not correcting avoids acting on an unreliable glitch "
            "classification. End glitches (a lone glitch at the first or last "
            "delta, with no pair partner) are always corrected regardless of "
            "this flag. Only meaningful when ``doDeglitch`` is True.",
    )
    rateCRiterMax = pexConfig.Field(
        dtype=int, default=5,
        doc="Maximum iteration count for the iterative UTR CR detector.",
    )

    firstRead = pexConfig.Field(
        dtype=int, default=None, optional=True,
        doc="0-indexed first read of the ramp to include (inclusive). "
            "None (default) = dispatch by observation type; an explicit value is used as-is.",
    )
    lastRead = pexConfig.Field(
        dtype=int, default=None, optional=True,
        doc="0-indexed last read of the ramp to include (inclusive; -1 = last read). "
            "None (default) = dispatch by observation type; an explicit value is used as-is.",
    )


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


def _stampBorderMask(exposure, *, borderWidth: int = 4) -> None:
    """Stamp the BORDER + BAD mask planes on the outer ``borderWidth``-pixel ring.

    Registers a new ``BORDER`` plane on first use and ORs it into the
    border rows/cols on top, bottom, left, right of the exposure mask.
    The BAD bit is set as well so downstream code that skips BAD pixels
    keeps working; BORDER is the additive informational distinction
    between "edge" and "broken-mid-detector" pixels.

    Width defaults to 4 to match ``h4Linearity.BORDER_PIX`` placement in
    ``h4Linearity/fit.py``.
    """
    if "BORDER" not in exposure.mask.getMaskPlaneDict():
        exposure.mask.addMaskPlane("BORDER")
        afwDisplay.setDefaultMaskPlaneColor("BORDER", "CYAN")
    borderBit = exposure.mask.getPlaneBitMask("BORDER")
    badBit = exposure.mask.getPlaneBitMask("BAD")
    bits = exposure.mask.array.dtype.type(borderBit | badBit)
    exposure.mask.array[:borderWidth, :] |= bits
    exposure.mask.array[-borderWidth:, :] |= bits
    exposure.mask.array[:, :borderWidth] |= bits
    exposure.mask.array[:, -borderWidth:] |= bits


def _stampReadRangeMetadata(exposure, *, r0, r1, nTotal, applyUTRWeights):
    """Add H4 read-range header keys to ``exposure.metadata``.

    Round-trips into the FITS header on persist, so downstream code (and
    humans inspecting an exposure) can tell a partial-ramp postISRCCD
    apart from a full-ramp one.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Modified in place.
    r0, r1 : int
        First/last reads used, absolute 0-indexed, inclusive.
    nTotal : int
        Total read count of the original ramp on disk (``pfsRaw.getNumReads()``).
    applyUTRWeights : bool
        Whether the image is UTR-weighted (True) or CDS-style (False).
    """
    md = exposure.getMetadata()
    md.set('H4READ0', int(r0), 'First H4 read used (absolute, 0-indexed)')
    md.set('H4READ1', int(r1), 'Last H4 read used (absolute, 0-indexed, inclusive)')
    md.set('H4NREAD', int(r1 - r0 + 1), 'Number of H4 reads spanned by this exposure')
    md.set('H4NTOT', int(nTotal), 'Total H4 reads in the original ramp')
    md.set('H4UTRWT', bool(applyUTRWeights),
           'applyUTRWeights state (True=UTR-weighted image, False=CDS)')


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
        else:
            linearity = None
        # All 3-d ramp-based operations are hidden inside this call. After .makeNirExposure()
        # we operate on a 2-d image.
        exposure, rawRamp = self.makeNirExposure(pfsRaw, nirDark,
                                                 linearity=linearity,
                                                 defects=defects,
                                                 doReturnRawCube=self.config.h4.doWriteRawCube)

        # We BAD mask the defects now, but do not interpolate.
        if defects is not None:
            super().maskDefect(exposure, defects)

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

        # Discriminate h4Linearity files (carrying a MODEL keyword whose
        # value is in h4Linearity.MODEL_REGISTRY) from the legacy
        # nirLinearity.NirLinearity format by inspecting the primary header.
        if h4Linearity.isH4LinearityFile(absFilename):
            return h4Linearity.loadFits(absFilename)
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

    def getDarkCube(self, nirDark, nreads: Optional[int] = None, r0: int = 0) -> np.ndarray:
        """Get the dark cube for the NIR ramp.

        Parameters
        ----------
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`
            The unprocessed to subtract.
        nreads : `int`, optional
            Number of dark frames to return. If None, return all from ``r0``.
        r0 : `int`
            0-indexed offset into the dark cube. Returns ``dark[r0:r0+nreads]``.
            Used when processing a sub-range of the data ramp; the dark slice
            then aligns with the data reads being processed.

        Returns
        -------
        `np.ndarray`
            The dark cube to subtract, shape ``(nreads, H, W)``.
        """

        if r0 > 0:
            end = nirDark.nreads if nreads is None else r0 + nreads
            full = nirDark.getImageCube(nreads=end)
            cube = full[r0:end].copy()
        else:
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

    def rampParams(self, pfsRaw):
        """Ramp-processing parameters, dispatched by observation type.

        Per-obstype defaults: arcs (``comparison``) → full-ramp CDS;
        darks → full-ramp UTR; science → UTR over reads ``1:-3`` (drop
        the shutter-closed read 0 and the trailing 3 shutter-closed /
        transitional reads). Anything else falls back to full-ramp UTR.

        Each `H4Config` field (``quickCDS`` / ``firstRead`` / ``lastRead``)
        that is set to an explicit (non-``None``) value overrides the
        dispatched default.

        Returns
        -------
        quickCDS : `bool`
        firstRead, lastRead : `int`
            0-indexed inclusive read bounds (negative counts from the end).
        """
        obsType = (pfsRaw.obsInfo.observation_type or "").lower()
        if obsType == "comparison":
            quickCDS, firstRead, lastRead = True, 0, -1
        elif obsType == "science":
            # Shutter open/close FITS cards are not yet written, so the
            # exact illuminated read range is unknown. Once those cards
            # are available the range can be computed correctly; until
            # then use reads 1:-3, which empirically drops the shutter-
            # closed read 0 and the trailing 3 shutter-closed reads.
            quickCDS, firstRead, lastRead = False, 1, -4
        else:  # dark, bias, flat, unknown — full-ramp UTR
            quickCDS, firstRead, lastRead = False, 0, -1

        cfg = self.config.h4
        if cfg.quickCDS is not None:
            quickCDS = cfg.quickCDS
        if cfg.firstRead is not None:
            firstRead = cfg.firstRead
        if cfg.lastRead is not None:
            lastRead = cfg.lastRead
        self.log.info(
            f"ramp params for obsType={obsType!r}: quickCDS={quickCDS} "
            f"firstRead={firstRead} lastRead={lastRead}"
        )
        return quickCDS, firstRead, lastRead

    def makeNirExposure(self,
                        pfsRaw,
                        nirDark=None,
                        linearity=None,
                        defects=None,
                        doReturnRawCube=False,
                        intermediates=None):
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

        # Dispatch CDS/UTR and the read range by observation type, then
        # resolve to absolute, 0-indexed inclusive bounds.
        quickCDS, firstRead, lastRead = self.rampParams(pfsRaw)
        r0 = pfsRaw.positiveIndex(firstRead)
        r1 = pfsRaw.positiveIndex(lastRead)
        if r1 - r0 < 1:
            raise ValueError(
                f"firstRead={firstRead} (->{r0}) and "
                f"lastRead={lastRead} (->{r1}) leave no readable range "
                f"(need r1 > r0)."
            )
        if quickCDS:
            # CDS from two reads, with two-read linearization: linearize
            # the absolute cumulative frame at each endpoint and
            # difference them (the model is nonlinear, so linearize then
            # subtract — never the reverse). This still yields the
            # linearity-derived mask planes (defects, bad fit,
            # saturation). The endpoint frames are dark-subtracted first
            # because the linearity model is calibrated on dark-
            # subtracted cumulative ADU. CR/glitch detection needs the
            # full ramp, so crResult stays None.
            self.log.info(f"creating CDS over reads [{r0}, {r1}].")
            frameR1 = self.makeCDS(pfsRaw, r0=0, r1=r1)
            frameR0 = self.makeCDS(pfsRaw, r0=0, r1=r0) if r0 > 0 else None
            if nirDark is not None:
                frameR1 -= self.getDarkRead(nirDark, r1 - 1)
                if frameR0 is not None:
                    frameR0 -= self.getDarkRead(nirDark, r0 - 1)

            flux = None
            crResult = None
            defectMask = None
            if self.config.h4.doLinearize and linearity is not None:
                self.log.info("Correcting non-linearity (two-read CDS).")
                linR1, _ = h4Linearity.applyFrame(linearity, frameR1)
                if frameR0 is not None:
                    linR0, _ = h4Linearity.applyFrame(linearity, frameR0)
                    nirImage = linR1 - linR0
                else:
                    nirImage = linR1
                # Merge the fit-time flags (defects, bad fit) with the
                # freshly computed out-of-range flags, as h4Linearity.apply
                # does for the full ramp.
                newMask = linearity.badPixelMask.copy()
                frames = [frameR1] if frameR0 is None else [frameR0, frameR1]
                for frame in frames:
                    newMask[frame > linearity.fitMax] |= h4Linearity.ABOVE_VALID_RANGE
                    newMask[frame < linearity.fitMin] |= h4Linearity.BELOW_VALID_RANGE
            else:
                nirImage = frameR1 if frameR0 is None else frameR1 - frameR0
                newMask = None
        else:
            # Too many interacting switches for CDS/UTR/darks. Especially note 2-d doDark is still possible.
            self.log.info(f"reading ramp over reads [{r0}, {r1}]...")
            # ``flux`` is the cumulative IRP-corrected ramp zero-anchored
            # at r0 — no diff→cumsum roundtrip; makeUTRcumulative returns
            # it directly.
            flux = self.makeUTRcumulative(pfsRaw, r0=r0, r1=r1)

            # makeUTRcumulative zero-anchors at r0, so flux[i] = read[r0+i+1] - read[r0].
            # Add the cumulative flux from read 0 to read r0 so flux is bias-relative
            # absolute (read[r0+i+1] - read[0]). Linearity is calibrated on absolute
            # cumulative ADU, so this anchoring is required for correctness when r0>0.
            offsetRaw = None
            if r0 > 0:
                self.log.info(f"adding absolute baseline from reads [0, {r0}].")
                offsetRaw = self.makeCDS(pfsRaw, r0=0, r1=r0).astype(flux.dtype, copy=False)
                flux += offsetRaw[None]

            if intermediates is not None:
                intermediates['raw'] = flux.copy()  # absolute cumulative, pre-dark

            if self.config.h4.applyUTRWeights:
                if nirDark is not None:
                    self.log.info("subtracting dark cube.")
                    darkCube = self.getDarkCube(nirDark, nreads=len(flux), r0=r0)
                    flux -= darkCube
                    del darkCube  # ~2.4 GB on a 4096²×40 ramp; no longer needed

            if intermediates is not None:
                intermediates['darkSubbed'] = flux.copy()  # input to linearity

            if self.config.h4.doLinearize and linearity is not None:
                self.log.info("Correcting non-linearity.")
                if defects is not None:
                    mi = afwImage.MaskedImageF(geom.Extent2I(4096, 4096))
                    defects.maskPixels(mi)
                    defectMask = (mi.mask.array > 0)
                    self.log.info(f"   starting with {defectMask.sum()} defect pixels")
                else:
                    defectMask = None

                ramp = h4Linearity.Ramp(reads=flux, validMask=defectMask)
                linearizedRamp = h4Linearity.apply(linearity, ramp)
                flux = linearizedRamp.cumulativeLinear
                newMask = linearizedRamp.badPixelMask
                # The pre-linearization cube is now superseded by the
                # linearized output. Drop the Ramp / LinearizedRamp
                # holders so the input (~2.4 GB) can be GC'd before the
                # CR/glitch detector copies a fresh ``cubeOriginal``.
                del ramp, linearizedRamp

                # Re-anchor the linearized cube at r0: subtract the linearized value
                # at read r0 from every read so the output represents only the flux
                # accumulated during (r0, r1]. Without this rebase, postISRCCD over
                # disjoint sub-ranges would not sum to the full-ramp postISRCCD.
                # No-op for r0=0 because the linearity model maps 0 -> 0.
                if r0 > 0 and offsetRaw is not None:
                    self.log.info(f"re-anchoring linearized cube at read {r0}.")
                    if nirDark is not None and self.config.h4.applyUTRWeights:
                        # Linearization saw flux already dark-subtracted; the rebase
                        # offset must match. offsetRaw is the absolute cumulative
                        # at read r0, so its dark is nirDark[r0-1] (nirDark[j] is
                        # the dark for read j+1 — the same k->k-1 indexing
                        # getDarkCube uses).
                        offsetForLin = offsetRaw - self.getDarkRead(nirDark, r0 - 1)
                    else:
                        offsetForLin = offsetRaw
                    linearizedOffset, _ = h4Linearity.applyFrame(linearity, offsetForLin)
                    flux -= linearizedOffset[None]

                if intermediates is not None:
                    intermediates['linearized'] = flux.copy()  # post-lin, pre-CR

                # Switch from cumulative flux to delta-space for everything
                # downstream of linearization. The CR/glitch detector
                # operates on deltas directly (no diff/cumsum dance inside),
                # and the science 2-D image is just ``deltas.sum(axis=0)``
                # — for linearized data the per-pixel mean of deltas IS
                # the optimal UTR rate. We only reconstruct the cumulative
                # cube on demand (intermediates / doReturnRawCube).
                read0 = flux[0:1].copy()
                deltas = np.diff(flux, axis=0)
                del flux
                iterResult = None  # filled in if doCR runs

                if self.config.h4.doCR:
                    crGood = (newMask == 0)
                    if defectMask is not None:
                        crGood &= ~defectMask

                    if self.config.h4.doDeglitch:
                        self.log.info("Correcting CRs and ASIC glitches.")
                        # Glitch detection runs on all pixels. The matched-pair
                        # cancellation criterion (opposite signs, sum within
                        # threshold) is what discriminates a glitch from a CR;
                        # restricting by ASIC channel misses glitches on
                        # channels other than the historically-noted ones.
                        glitchChanMask = np.ones(deltas.shape[1:], dtype=bool)
                    else:
                        self.log.info("Correcting CRs (ASIC deglitching disabled).")
                        glitchChanMask = None
                    iterResult = h4Linearity.cr.iterativeUtrDetectAndRepair(
                        deltas,
                        goodPixelMask=crGood,
                        glitchPixelMask=glitchChanMask,
                        sigmaFloorADU=self.config.h4.rateCRsigmaFloorADU,
                        nSigma=self.config.h4.rateCRnSigma,
                        maxIterations=self.config.h4.rateCRiterMax,
                        correctGlitches=self.config.h4.correctGlitches,
                        glitchAmplitudeMinADU=self.config.h4.deglitchAmplitudeMinADU,
                    )
                    crFlagMask2D = iterResult.crFlagMask.any(axis=0)
                    glitchFlagMask2D = iterResult.glitchFlagMask.any(axis=0)
                    crResult = SimpleNamespace(
                        flagMask=crFlagMask2D,
                        nFlagged=int(crFlagMask2D.sum()),
                        glitchFlagMask=glitchFlagMask2D,
                        nGlitchFlagged=int(glitchFlagMask2D.sum()),
                    )
                    self.log.info(
                        f"iterative CR/glitch step: "
                        f"{iterResult.nCRs} CR flag entries, "
                        f"{iterResult.nGlitchPairs} glitch pairs, "
                        f"{crResult.nFlagged} unique CR pixels, "
                        f"{crResult.nGlitchFlagged} unique glitch pixels "
                        f"in {iterResult.nIterations} iterations."
                    )
                else:
                    crResult = None

                # Reconstruct the cumulative cube ONLY when the debug
                # ``intermediates`` dict or ``doReturnRawCube`` asks for
                # it. The default production path is delta-only from
                # here on, which is the whole point of this restructure.
                flux = None
                if (intermediates is not None) or doReturnRawCube:
                    flux = np.empty(
                        (deltas.shape[0] + 1,) + deltas.shape[1:],
                        dtype=deltas.dtype,
                    )
                    flux[0:1] = read0
                    np.cumsum(deltas, axis=0, out=flux[1:])
                    flux[1:] += read0
                    if intermediates is not None:
                        intermediates['crCorrected'] = flux.copy()  # post-CR
            else:
                newMask = None
                crResult = None
                # No linearization → no delta-space switch. Leave flux
                # as-is (post-dark cumulative); fall through to the
                # legacy flux-based nirImage path below.
                deltas = None
                read0 = None

            if deltas is not None:
                # Delta-space science image: total accumulated flux is
                # the per-pixel UTR rate times the number of deltas.
                # On linearized deltas the UTR weights collapse to
                # uniform, so the optimal rate is just the median/mean
                # of the repaired deltas — which the CR iteration
                # already produced as ``iterResult.rate`` (or
                # equivalently the median saved in ``rateFull``).
                # Skipping a separate ``deltas.sum(axis=0)`` saves one
                # full-cube reduction on the (N-1, H, W) buffer.
                if self.config.h4.doCR and iterResult is not None:
                    pixelRate = iterResult.rate
                else:
                    pixelRate = deltas.mean(axis=0)
                nirImage = pixelRate * deltas.shape[0]
                if not self.config.h4.applyUTRWeights and nirDark is not None:
                    # CDS mode: dark cube was NOT subtracted upstream,
                    # so apply the 2-D dark CDS correction here.
                    nirImage = nirImage - (
                        self.getDarkRead(nirDark, r1 - 1)
                        - self.getDarkRead(nirDark, r0)
                    )
                del deltas, read0
            elif self.config.h4.applyUTRWeights:
                self.log.info("applying UTR weights.")
                rates = self.calcUTRrates(flux)
                nirImage = rates * len(flux)
            else:
                # CDS, basically. flux[i] aligns with nirDark[r0+i] (the dark
                # cube has one entry per cumsum frame, indexed as
                # nirDark[i]==dark-for-read[i+1]). So flux[0] uses nirDark[r0]
                # and flux[-1] uses nirDark[r1-1]. For the full-ramp default
                # (r0=0, r1=N-1) this reproduces the prior 0 and len(flux)-1.
                if nirDark is not None:
                    nirImage = ((flux[-1] - self.getDarkRead(nirDark, r1 - 1)) -
                                (flux[0] - self.getDarkRead(nirDark, r0)))
                else:
                    nirImage = flux[-1] - flux[0]

        exposure = self._makeExposure(pfsRaw, nirImage)
        _stampReadRangeMetadata(
            exposure, r0=r0, r1=r1,
            nTotal=int(pfsRaw.getNumReads()),
            applyUTRWeights=bool(self.config.h4.applyUTRWeights),
        )
        # BORDER + BAD on the outer 4-pixel ring, unconditionally. Done
        # before the linearity-mask sweep below so the BAD bit is set even
        # when doLinearize=False (in which case ``newMask`` is None and
        # the sweep block is skipped). When linearization HAS run, the
        # sweep also catches h4Linearity.BORDER_PIX and folds it into BAD;
        # that's harmless overlap with the same bits we set here.
        _stampBorderMask(exposure)
        if newMask is not None:
            # Declare that all pixels above the available linearization
            # limits are SATurated.
            saturated = (newMask & h4Linearity.ABOVE_VALID_RANGE) > 0
            exposure.mask.array[saturated] |= exposure.mask.getPlaneBitMask('SAT')

            # Declare that everything else is BAD. Note that this gets the border pixels. We
            # no not want to do that, but they are causing trouble downstream so we mark them
            # here for now.
            lowVal = (newMask & h4Linearity.BELOW_VALID_RANGE) > 0
            badMask = (newMask & ~saturated & ~lowVal) > 0
            exposure.mask.array[badMask] |= exposure.mask.getPlaneBitMask('BAD')

            defectMask2 = (newMask & h4Linearity.MASKED_BY_INPUT) > 0
            exposure.mask.array[defectMask2] |= exposure.mask.getPlaneBitMask('BAD')
            nCR = crResult.nFlagged if crResult is not None else 0
            self.log.info(f'nSat={saturated.sum()} '
                          f'nLow={lowVal.sum()} '
                          f'nBad={badMask.sum()} '
                          f'nDefects={"none" if defectMask is None else defectMask.sum()} '
                          f'nDefects2={defectMask2.sum()} '
                          f'nCR={nCR}')

        if crResult is not None and crResult.nFlagged > 0:
            exposure.mask.array[crResult.flagMask] |= exposure.mask.getPlaneBitMask("CR")

        glitchFlagMask = getattr(crResult, 'glitchFlagMask', None)
        if glitchFlagMask is not None and glitchFlagMask.any():
            if "ASIC_GLITCH" not in exposure.mask.getMaskPlaneDict():
                exposure.mask.addMaskPlane("ASIC_GLITCH")
                afwDisplay.setDefaultMaskPlaneColor("ASIC_GLITCH", "ORANGE")
            exposure.mask.array[glitchFlagMask] |= exposure.mask.getPlaneBitMask("ASIC_GLITCH")

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

    def asicBadChannelMask(self, detectorName: str, shape: tuple,
                           nChannels: int = 32) -> np.ndarray:
        """Return a ``(H, W)`` bool mask True at rows in known-bad ASIC channels.

        H4 ASIC channels run along the y (rows) axis; each channel
        spans ``H // nChannels`` rows. Useful as the ``glitchPixelMask``
        argument to ``cr.iterativeUtrDetectAndRepair``: glitch detection
        is restricted to the rows where ASIC glitches are known to occur.

        Parameters
        ----------
        detectorName : str
            E.g. ``"n3"``.
        shape : tuple
            ``(H, W)`` detector dimensions.
        nChannels : int
            Number of horizontal ASIC channels stacked along Y; 32 for H4.

        Returns
        -------
        mask : np.ndarray
            ``(H, W)`` bool. All-False if the detector has no known-bad
            channels.
        """
        channels = self.loadBadAsicChannels(detectorName)
        H, W = shape
        channelHeight = H // nChannels
        mask = np.zeros((H, W), dtype=bool)
        for ch in channels:
            mask[ch * channelHeight:(ch + 1) * channelHeight, :] = True
        return mask

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

    def makeUTRcumulative(self, pfsRaw, r0=0, r1=-1, nreads=None, bbox=None,
                          showTimes=False) -> np.ndarray:
        """Return the IRP-corrected cumulative ramp as a single 3-D stack.

        Given two raw data images d0 and d1, and two raw IRP images i0 and i1, the net CDS image
        can be either (d1 - i1) - (d0 - i0), or (d1 - d0) - (i1 - i0). The IRP row has various
        artifacts which make using the latter "nicer", or at least easier to make sense of. So that
        is the way we do it.
        In particular there are:
          - bad IRP row pixels
          - fixed pixel-to-pixel offsets

        ASIC glitches are handled later by the iterative CR/glitch detector
        operating on the linearized cube, not here.

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
        flux : 3-d float32 numpy array
           Cumulative IRP-corrected ADU, zero-anchored at read ``r0``.
           ``flux[k]`` is read ``r0+k+1`` minus read ``r0`` minus the
           per-channel IRP-filtered diff. Caller uses this directly as
           the bias-relative cumulative ramp — no diff/cumsum roundtrip
           needed.
        """

        r0 = pfsRaw.positiveIndex(r0)
        r1 = pfsRaw.positiveIndex(r1)
        if r1-r0 < 1:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}) + 1:')
        if nreads is None:
            nreads = r1 - r0 + 1
        reads = np.linspace(r0, r1, nreads, dtype='i2')

        # Grab the components of read r0, which we will subtract from all the others.
        # ASIC glitches are now handled later by the iterative CR/glitch detector
        # operating on the linearized cube, so no longer need to preload raw data
        # and IRP cubes here for pre-linearization spike repair.
        data0 = self.makeRawDataArray(pfsRaw, r0)
        irp0 = self.makeRawIrpArray(pfsRaw, r0)
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
            data1 = self.makeRawDataArray(pfsRaw, r_i)
            irp1 = self.makeRawIrpArray(pfsRaw, r_i)
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

        # ``stack`` is already the cumulative IRP-corrected ramp the
        # caller wants. The previous implementation did an extra
        # ``np.diff(..., prepend=zeros)`` here only to have the caller
        # ``np.cumsum`` it back — an exact roundtrip costing 5.85 GB of
        # transient alloc and ~10 s of CPU on a 4096²×88 ramp.
        return stack

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
