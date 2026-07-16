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
from typing import Optional

import os
import time
import warnings

from functools import partial
from types import SimpleNamespace

import numpy as np
import scipy

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
from . import h4Linearity
from .overscan import PfsOverscanCorrectionTask
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


class H4Config(pexConfig.Config):
    """Configuration parameters for H4 reductions"""

    quickCDS = pexConfig.Field(
        dtype=bool, default=None, optional=True,
        doc="Only consider last and first reads instead of the full ramp cube. "
            "None (default) = dispatch by observation type "
            "(CDS for arcs and flats; UTR for darks and science).")
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
            "the linearized cube after dark subtraction.",
    )
    repairCR = pexConfig.Field(
        dtype=bool, default=True,
        doc="When ``doCR`` is True, replace flagged deltas in the "
            "cube with the per-pixel rate so the integrated image "
            "reflects the CR-corrected ramp. Set False to leave the "
            "raw flagged values in place (diagnostic use; the rate "
            "image is still computed from the detector's UTR-weighted "
            "estimator, but the cube exposed via ``intermediates`` / "
            "``doReturnRawCube`` shows uncorrected reads).",
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
    rateCRmaxDropFraction = pexConfig.Field(
        dtype=float, default=0.5,
        doc="Cumulative-drop check for CR classification — amplitude "
            "criterion. A candidate CR fails this criterion when the "
            "per-pixel running cumulative residual from the CR delta "
            "onward drops by more than this fraction of the CR "
            "amplitude. Default 0.5; set very high to disable. The "
            "candidate is REJECTED only when BOTH this and "
            "``rateCRnDropSigma`` (the noise criterion) fail.",
    )
    rateCRnDropSigma = pexConfig.Field(
        dtype=float, default=3.0,
        doc="Cumulative-drop check for CR classification — noise "
            "criterion. A candidate CR fails this criterion when the "
            "cumulative residual drop exceeds ``nDropSigma * σ * "
            "sqrt(N-1-k)`` ADU. Real CRs leave only random-walk noise "
            "behind, with an expected cumulative drift ~σ√(N-1-k) — "
            "that drift is not a transient signature and shouldn't "
            "demote the candidate. Transients leak the deposited "
            "charge back, producing a drop that vastly exceeds the "
            "noise floor. The candidate is REJECTED only when both "
            "this AND ``rateCRmaxDropFraction`` fail. Default 3.0; "
            "set very high to fall back to the amplitude criterion "
            "alone.",
    )
    badPixelMinOutliers = pexConfig.Field(
        dtype=int, default=4,
        doc="BAD-pixel gate inside the CR detector (count criterion). "
            "A pixel is OR'd into the UNSTABLE + BAD mask planes when "
            "its ramp has at least this many delta residuals exceeding "
            "``badPixelOutlierSigma × σ_IQR`` from the per-pixel "
            "median. Set 0 to disable.",
    )
    badPixelOutlierSigma = pexConfig.Field(
        dtype=float, default=4.0,
        doc="BAD-pixel gate inside the CR detector (sigma criterion). "
            "Per-delta outlier threshold in units of ``σ_IQR``. Default "
            "4.0 — 3.0 catches many merely-noisy pixels because the "
            "IQR-σ underestimates the true scatter by ~30 % on real "
            "H4 deltas (mild non-Gaussian tails from shot/read-noise "
            "mixture + linearity residuals).",
    )
    asicGlitchHeightMaskADU = pexConfig.Field(
        dtype=float, default=0.0,
        doc="ASIC-glitch repair/mask threshold. Pixels with at least "
            "one detected glitch pair whose height ``|A| = (δ[k] − "
            "δ[k+1]) / 2`` exceeds this value (ADU) are OR'd into BAD "
            "on the published mask. Default 0 = mask every detected "
            "glitch pixel. Set to a few hundred ADU to leave the "
            "noise-pair detection tail un-masked while still catching "
            "genuine bit-flip events.",
    )

    doRateStability = pexConfig.Field(
        dtype=bool, default=True,
        doc="Run the per-pixel rate-stability rejection on the "
            "linearized, CR/glitch-repaired UTR delta cube. Splits the "
            "ramp in half, compares the per-half mean delta (the "
            "production rate estimator), and OR's RATE_UNSTABLE into "
            "the internal mask for pixels whose halves disagree "
            "fractionally above the threshold. Inactive on the CDS "
            "path or when linearization is skipped.",
    )
    rateStabilityThreshold = pexConfig.Field(
        dtype=float, default=0.20,
        doc="Rejection threshold for the rate-stability test, applied "
            "to ``sqrt(max(0, (r1-r2)² − sem1² − sem2²)) / max(|r1|, "
            "|r2|, rateStabilityRateFloorADU)``. Pixels above this "
            "fraction are flagged RATE_UNSTABLE.",
    )
    rateStabilityRateFloorADU = pexConfig.Field(
        dtype=float, default=5.0,
        doc="Lower bound on the denominator of the rate-stability "
            "fraction (ADU/read). Protects dark / near-zero-rate "
            "pixels from blowing up the fraction from noise.",
    )
    rateStabilityMinDeltasPerSegment = pexConfig.Field(
        dtype=int, default=3,
        doc="Minimum number of un-flagged deltas for a half to count "
            "as testable in the rate-stability test.",
    )
    maskRateUnstable = pexConfig.Field(
        dtype=bool, default=False,
        doc="If True, pixels with the internal RATE_UNSTABLE bit set "
            "are also OR'd into the published UNSTABLE and BAD planes. "
            "If False (default), RATE_UNSTABLE-only pixels are flagged "
            "in the dedicated RATE_UNSTABLE plane but kept out of "
            "UNSTABLE and BAD — downstream consumers see the marking "
            "without losing the pixel. A pixel that is RATE_UNSTABLE and "
            "also carries another internal bit still lands in BAD via "
            "that other bit regardless of this knob.",
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


def _makeInternalMask(
    shape, *, linearity=None, defects=None, borderWidth: int = 4,
):
    """Build the H4 ISR internal mask (uint16, ``(H, W)``).

    The internal-mask alphabet is the union of:

      - ``h4Linearity.BORDER_PIX`` on the outer ``borderWidth`` ring,
        seeded here so even doLinearize=False paths get it.
      - ``h4Linearity.MASKED_BY_INPUT`` on pixels carried by the
        ``defects`` calib.
      - Any bits already set in ``linearity.badPixelMask``
        (typically also BORDER_PIX + MASKED_BY_INPUT + DEAD-group bits
        the linearity fit() emitted).

    Downstream stages (``h4Linearity.apply``, the CR detector,
    ``cr.iterativeUtrDetectAndRepair``) read this mask, skip already-
    flagged pixels, and OR their own findings back in.
    ``_projectInternalMask`` lifts the final mask into ``Exposure.mask``
    planes at the end.
    """
    H, W = shape
    internal = np.zeros((H, W), dtype=np.uint16)
    internal[:borderWidth, :] |= h4Linearity.BORDER_PIX
    internal[-borderWidth:, :] |= h4Linearity.BORDER_PIX
    internal[:, :borderWidth] |= h4Linearity.BORDER_PIX
    internal[:, -borderWidth:] |= h4Linearity.BORDER_PIX
    if defects is not None:
        defImg = afwImage.MaskedImageF(geom.Extent2I(W, H))
        defects.maskPixels(defImg)
        internal[defImg.mask.array > 0] |= h4Linearity.MASKED_BY_INPUT
    if linearity is not None:
        internal |= linearity.badPixelMask.astype(np.uint16, copy=False)
    return internal


def _projectInternalMask(exposure, internalMask, *, crResult=None,
                         maskRateUnstable: bool = False) -> None:
    """Lift the H4 internal mask into ``Exposure.mask`` planes.

    Single projection point for the canonical published set:

      - ``DARK_DEFECT``      ← ``MASKED_BY_INPUT``                  (also BAD)
      - ``LINEARITY_DEFECT`` ← ``DEAD`` group                       (also BAD)
      - ``SAT``              ← ``ABOVE_VALID_RANGE``                (also BAD)
      - ``UNSTABLE``         ← CR-stage ``UNSTABLE`` always;
        ∪ ``RATE_UNSTABLE`` only when ``maskRateUnstable=True``     (also BAD)
      - ``RATE_UNSTABLE``    ← ``RATE_UNSTABLE`` always (independent
        of ``maskRateUnstable``); included in BAD only when
        ``maskRateUnstable=True``
      - ``ASIC_GLITCH``      ← ``ASIC_GLITCH`` always (NOT BAD on its
        own; a corrected glitch carries only this, an uncorrectable one
        also carries ``GLITCH_MASKED`` and is BAD)
      - ``BAD`` ← anywhere any internal bit is set, except a lone
        ``ASIC_GLITCH`` (always) or a lone ``RATE_UNSTABLE`` (when
        ``maskRateUnstable=False``); so BORDER / BELOW_VALID_RANGE etc.
        land in BAD without an externally distinguished plane.
      - ``CR`` ← from ``crResult.crFlagMask`` (per-delta).

    Per the "first-reason-wins" rule, ABOVE_VALID_RANGE only fires on
    pixels that survived defects + fit, so SAT is now the clean
    "genuinely saturating" signal rather than a side effect of dead
    pixels with stale ``fitMax``. ``UNCLASSIFIED`` is not published as a
    standalone plane -- it lives on ``IterativeRepairResult`` for internal
    diagnostic use and folds into BAD via the catch-all. The CR-stage
    ``UNSTABLE`` bit (outlier-count gate) always projects to the published
    ``UNSTABLE`` plane; the rate-stability ``RATE_UNSTABLE`` bit
    (half-vs-half disagreement) gets its own published plane and is *also*
    added to UNSTABLE + BAD only when the caller passes
    ``maskRateUnstable=True``. ``ASIC_GLITCH`` is a published, non-BAD
    plane marking every glitch pixel (corrected or masked). A pixel that
    carries another internal bit still lands in BAD via that other bit
    regardless of ``maskRateUnstable`` or ASIC_GLITCH.
    """
    mask = exposure.mask
    for plane in ("DARK_DEFECT", "LINEARITY_DEFECT", "UNSTABLE",
                  "RATE_UNSTABLE", "ASIC_GLITCH"):
        if plane not in mask.getMaskPlaneDict():
            mask.addMaskPlane(plane)
    bit = mask.getPlaneBitMask

    darkDefect = (internalMask & h4Linearity.MASKED_BY_INPUT) != 0
    linDefect = (internalMask & h4Linearity.DEAD) != 0
    sat = (internalMask & h4Linearity.ABOVE_VALID_RANGE) != 0
    rateUnstable = (internalMask & h4Linearity.RATE_UNSTABLE) != 0
    unstable = (internalMask & h4Linearity.UNSTABLE) != 0
    asicGlitch = (internalMask & h4Linearity.ASIC_GLITCH) != 0
    # ASIC_GLITCH alone (a corrected glitch) is usable and never implies
    # BAD; an uncorrectable glitch also carries GLITCH_MASKED. RATE_UNSTABLE
    # alone is BAD only when maskRateUnstable is set.
    notBad = np.uint16(h4Linearity.ASIC_GLITCH)
    if maskRateUnstable:
        unstable = unstable | rateUnstable
    else:
        notBad = notBad | np.uint16(h4Linearity.RATE_UNSTABLE)
    anyBad = (internalMask & ~notBad) != 0

    arr = mask.array
    if darkDefect.any():
        arr[darkDefect] |= bit("DARK_DEFECT")
    if linDefect.any():
        arr[linDefect] |= bit("LINEARITY_DEFECT")
    if sat.any():
        arr[sat] |= bit("SAT")
    if unstable.any():
        arr[unstable] |= bit("UNSTABLE")
    if rateUnstable.any():
        arr[rateUnstable] |= bit("RATE_UNSTABLE")
    if asicGlitch.any():
        arr[asicGlitch] |= bit("ASIC_GLITCH")
    if anyBad.any():
        arr[anyBad] |= bit("BAD")

    if crResult is not None:
        crFlag = getattr(crResult, "flagMask", None)
        if crFlag is not None and crFlag.any():
            arr[crFlag] |= bit("CR")


def _stampRampMetadata(exposure, *, r0, r1, nTotal, appliedUTR):
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
    appliedUTR : bool
        Whether UTR weights were actually applied to produce
        ``exposure.image`` (True for the linearized + legacy-UTR arms,
        False for quickCDS and the non-linearized non-UTR fallback).
        Authoritative for downstream consumers — e.g. variance
        estimation chooses between the CDS and UTR noise formulas on
        this flag.
    """
    md = exposure.getMetadata()
    md.set('H4READ0', int(r0), 'First H4 read used (absolute, 0-indexed)')
    md.set('H4READ1', int(r1), 'Last H4 read used (absolute, 0-indexed, inclusive)')
    md.set('H4NREAD', int(r1 - r0 + 1), 'Number of H4 reads spanned by this exposure')
    md.set('H4NTOT', int(nTotal), 'Total H4 reads in the original ramp')
    md.set('H4UTRWT', bool(appliedUTR),
           'UTR weights applied to exposure image (False=CDS)')


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


def lookupNirDark(datasetType, registry, dataId, collections):
    """Resolve a NIR dark cube (``nirDark`` or ``nirDark_irp4``).

    Ramps taken with different dataPixel-to-IRP ratios (``W_H4IRPN``, 1 or 4)
    have different per-read times and so need different dark cubes. The
    ``nirDark`` and ``nirDark_irp4`` dataset types hold the otherwise-identical
    cubes for ratios 1 and 4.

    The registry does not (yet) carry the ratio, so the choice cannot be made
    here at qgraph-build time. This shared lookup therefore resolves *whichever*
    type it is asked for whenever that dataset exists; both connections resolve
    independently. The ratio-matched dark is chosen later, during runtime
    initialization, by `selectNirDark` from the loaded ramp's ``W_H4IRPN``.

    Parameters
    ----------
    datasetType : `lsst.daf.butler.DatasetType`
        The dataset type to look up (``nirDark`` or ``nirDark_irp4``).
    registry : `lsst.daf.butler.Registry`
        The butler registry.
    dataId : `lsst.daf.butler.DataCoordinate`
        The data identifier.
    collections : `list` of `str`
        The collections to search.

    Returns
    -------
    refs : `list` of `lsst.daf.butler.Reference`
        A single reference to the dark cube, or empty if none exists.
    """
    ref = registry.findDataset(
        datasetType, collections=collections, dataId=dataId, timespan=dataId.timespan
    )
    return [ref] if ref is not None else []


def selectNirDark(inputs):
    """Pick the ratio-matched NIR dark and collapse it to a single ``nirDark``.

    Interim runtime hack while the registry lacks the IRP-ratio dimension:
    `lookupNirDark` resolves both ``nirDark`` (ratio 1) and ``nirDark_irp4``
    (ratio 4) into ``inputs``; here we choose the one matching the loaded ramp's
    ``W_H4IRPN`` (``PfsRaw.irpN``) and leave the rest of the task ratio-agnostic
    — it sees only ``nirDark``. For a non-NIR (CCD) exposure neither applies and
    both are dropped.

    When ``irp_ratio`` becomes a populated registry dimension, `lookupNirDark`
    can select at build time and this collapses to a no-op.

    Parameters
    ----------
    inputs : `dict`
        The task inputs from ``butlerQC.get``; mutated in place. Must contain
        ``"ccdExposure"`` (a `~lsst.obs.pfs.PfsRaw`); ``"nirDark"`` and
        ``"nirDarkIrp4"`` are consumed when present.
    """
    raw = inputs["ccdExposure"]
    nirDark1 = inputs.pop("nirDark", None)
    nirDark4 = inputs.pop("nirDarkIrp4", None)
    if raw.isNir():
        inputs["nirDark"] = nirDark4 if raw.irpN == 4 else nirDark1


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
        doc="Input dark calibration for NIR (IRP ratio 1)",
        storageClass="ImageCube",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        lookupFunction=lookupNirDark,
        minimum=0,  # allowed to not exist, since we may use dark instead
    )
    badRefPixels = PrerequisiteConnection(
        name="badRefPixels",
        doc="Bad IRP reference-row pixels (per-detector)",
        storageClass="NirBadRefPixels",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        minimum=0,  # only required for NIR when h4.doIRPbadPixels
    )
    linearity = PrerequisiteConnection(
        name="h4Linearity",
        doc="H4 NIR nonlinearity correction (per-detector)",
        storageClass="H4Linearity",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        minimum=0,  # only required for NIR when h4.doLinearize
    )
    nirDarkIrp4 = PrerequisiteConnection(
        name="nirDark_irp4",
        doc="Input dark calibration for NIR (IRP ratio 4)",
        storageClass="ImageCube",
        dimensions=["instrument", "arm", "spectrograph"],
        isCalibration=True,
        lookupFunction=lookupNirDark,
        minimum=0,  # allowed to not exist; resolved only for IRP ratio 4 ramps
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
            self.prerequisiteInputs.remove("nirDarkIrp4")
        if config.h4.doIRPbadPixels is not True:
            self.prerequisiteInputs.remove("badRefPixels")
        if config.h4.doLinearize is not True:
            self.prerequisiteInputs.remove("linearity")
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

        # Choose the ratio-matched NIR dark from the ramp's W_H4IRPN and
        # collapse to a single "nirDark" (see selectNirDark); the rest of the
        # task is ratio-agnostic.
        selectNirDark(inputs)

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
            if isNir:
                if inputs.get("nirDark") is None:
                    raise RuntimeError(f"No NIR dark cube found for {raw.detector.getName()}")
            elif inputs.get("dark") is None:
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
        nirDark=None,
        flat=None,
        defects=None,
        detectorNum=None,
        ipcCoeffs=None,
        badRefPixels=None,
        linearity=None,
        **kwargs,
    ):
        """Specialist instrument signature removal for H4RG detectors

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
            The raw exposure that is to be run through ISR.  The
            exposure is modified by this method. With the PfsRaw we
            can get access the ramp cubes.
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`, optional
            Dark cube to subtract.
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

        # Stash for correctBadIrpPixels(), which is reached deep in the ramp path.
        self._badRefPixels = badRefPixels

        if self.config.doDark:
            if nirDark is None:
                raise RuntimeError("Must supply a NIR dark cube if config.doDark=True.")
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
            if linearity is None:
                self.log.warn(
                    f'no usable linearity for {pfsRaw.detector.getName()}; '
                    'proceeding without linearity correction for this exposure'
                )
            elif defects is None:
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
        # CDS vs UTR variance. Driven by ``H4UTRWT`` (records whether
        # UTR weights were actually applied during ``makeNirExposure``)
        # and ``H4NREAD`` (the actual reads integrated after the
        # firstRead/lastRead trim — *not* ``pfsRaw.getNumReads()``).
        md = exposure.getMetadata()
        if md.getScalar("H4UTRWT"):
            nread = int(md.getScalar("H4NREAD"))
            # UTR noise per RHL Eq. 4.45.
            var *= 6 * (nread * nread + 1) / (5 * nread * (nread + 1))
            var += (12 * (nread - 1) / (nread * (nread + 1))
                    * channel.getReadNoise()**2)
        else:
            var += 2 * channel.getReadNoise()**2  # CDS
        exposure.variance.array[:] = var

        if rawRamp is not None:
            rawRamp *= gain
            # The cube is now in electrons; label it with the applied gain so
            # consumers can recover the e-/ADU scaling. deepCopy so we don't
            # mutate the metadata shared with the (also-electron) postISRCCD.
            rawMd = exposure.getMetadata().deepCopy()
            rawMd.set("GAIN", gain)
            rawRamp = imageCube.ImageCube.fromCube(rawRamp, rawMd)

        # Any nquarter stuff should be removed after PIPE2D-1200
        nQuarter = exposure.getDetector().getOrientation().getNQuarter()

        if self.config.h4.doIPC:
            self.log.info("Applying IPC correction.")
            self.correctIPC(exposure, defects, ipcCoeffs, -nQuarter)

        if self.config.doDefect:
            self.log.info("Masking defects.")
            self._maskDefects(exposure, defects)

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
        """Load the shipped H4 linearity correction for a camera from disk.

        This file-based loader exists for offline validation tooling
        (`lsst.obs.pfs.h4Linearity.validate`). The ISR pipeline itself
        consumes the ``h4Linearity`` calibration through the butler
        (the ``linearity`` prerequisite input), not this method.

        Parameters
        ----------
        cam : `str`
           The camera name.

        Returns
        -------
        linearity : `lsst.obs.pfs.h4Linearity.LinearityCorrection` or `None`
           The linearity correction, or `None` if it is missing or not in the
           h4Linearity FITS format.
        """
        filename = f'nirLinearity-{cam}.fits'
        absFilename = os.path.join(getPackageDir("drp_pfs_data"), "nirLinearity", filename)
        if not os.path.exists(absFilename):
            self.log.warn(f'no linearity available for {cam}: {absFilename} not found')
            return None
        if not h4Linearity.isH4LinearityFile(absFilename):
            self.log.warn(
                f'no usable linearity for {cam}: {absFilename} is not an h4Linearity-format file'
            )
            return None
        return h4Linearity.loadFits(absFilename)

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

    @staticmethod
    def calcUTRrateFromDeltas(deltas: np.ndarray) -> np.ndarray:
        """UTR-weighted per-pixel rate, computed directly from a delta
        cube. Equivalent to :meth:`calcUTRrates` applied to the
        reconstructed cumulative ramp ``read0 + cumsum(deltas, axis=-1)``,
        but skips the cumulative reconstruction by using the closed-form
        delta weights ``u[j] = 6(j+1)(N-1-j) / (N(N-1)(N+1))`` for
        j = 0..N-2. These sum to 1 and reproduce the read-space UTR
        weights ``w[i] = (12i - 6(N-1)) / (N(N²-1))`` after the
        cumsum/diff change of variables.

        Independent of task state; exposed as a staticmethod so
        diagnostic code can call the production rate formula without
        instantiating ``PfsIsrTask``.

        Parameters
        ----------
        deltas : `np.ndarray`
            ``(H, W, N-1)`` per-pixel delta cube with the time axis last.
            Also accepts lower-rank arrays as long as the time axis is
            last; the rate has shape ``deltas.shape[:-1]``.

        Returns
        -------
        `np.ndarray`
            Per-pixel rate (ADU/read), shape ``deltas.shape[:-1]``.
        """
        nDeltas = deltas.shape[-1]
        nReads = nDeltas + 1
        ks = np.arange(nDeltas, dtype=np.float32)
        weights = (
            6.0 * (ks + 1.0) * (nReads - 1.0 - ks)
            / (nReads * (nReads - 1.0) * (nReads + 1.0))
        ).astype(np.float32, copy=False)
        # Slice-wise accumulation avoids an (H, W, N-1) transient.
        rate = np.zeros(deltas.shape[:-1], dtype=np.float32)
        for k in range(nDeltas):
            rate += weights[k] * deltas[..., k]
        return rate

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

    def _darkGain(self, nirDark) -> float:
        """The nirDark's electrons/ADU gain, guarded against a bad label.

        ``combineNirDark`` labels the electron-valued dark with the applied
        gain, which the dark methods divide out to reach ADU. A non-physical
        value -- notably the ``9999`` raw-ASIC placeholder -- means the dark was
        never gain-labeled, so ``dark / gain`` would silently mis-scale the
        subtraction (``/9999`` leaves essentially no dark). Fail loudly rather
        than guess the units. ``1.0`` (no ``GAIN`` header) is allowed and means
        the dark is already in ADU.
        """
        gain = nirDark.metadata.get("GAIN", 1.0)
        if not (0.0 < gain < 100.0):
            raise RuntimeError(
                f"nirDark GAIN={gain} is non-physical (9999 is the raw ASIC "
                f"placeholder); rebuild the dark with combineNirDark."
            )
        return gain

    def subtractDarkCube(self, nirDark, cube: np.ndarray, r0: int = 0) -> None:
        """Subtract the per-read dark frames from ``cube`` in place.

        Iterates absolute read indices ``[r0, r0 + cube.shape[0])`` and
        subtracts each dark frame from the corresponding read of
        ``cube``. No transient ramp-sized buffer: each ``(H, W)`` dark
        frame is fetched via ``nirDark.getReadArray(...)``, gain-corrected,
        and applied directly to ``cube[k]``. The per-read subtract is a
        contiguous 2-D write because ``cube`` is ``(N, H, W)`` C-order,
        which is cache-friendly and avoids the ~14 s / ~6.7 GB cost of
        materializing-then-transposing a full dark cube via
        :meth:`getDarkCube`.

        Parameters
        ----------
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`
            The unloaded dark cube — frames live on disk and are
            fetched on demand via ``nirDark.getReadArray``.
            ``nirDark[i]`` is the dark for absolute read ``i + 1``
            (the convention :meth:`getDarkRead` and :meth:`getDarkCube`
            use).
        cube : `np.ndarray`
            ``(N, H, W)`` cumulative ramp to dark-subtract, modified in
            place. ``cube[k]`` corresponds to absolute read ``r0 + k + 1``
            of the original ramp.
        r0 : `int`
            Absolute index of the first read processed; ``nirDark[r0+k]``
            is paired with ``cube[k]``.
        """
        gain = self._darkGain(nirDark)
        N = cube.shape[0]
        for k in range(N):
            darkFrame = nirDark.getReadArray(r0 + k)
            if gain != 1.0:
                # Back out the gain that ``ImageCube`` applied when
                # writing the dark in electrons; the cube here is in ADU.
                cube[k] -= darkFrame / gain
            else:
                cube[k] -= darkFrame

    def getDarkCube(self, nirDark, nreads: Optional[int] = None, r0: int = 0) -> np.ndarray:
        """Get the dark cube for the NIR ramp.

        .. note::

           Production code uses :meth:`subtractDarkCube`, which fuses the
           per-read fetch with an in-place subtract — no transient
           dark-cube allocation, no transpose. ``getDarkCube`` is kept
           for diagnostic paths (``validate.collectPixelRampData``'s
           ``cubeDark`` and tests) that want to inspect the dark cube
           in the project's ``(H, W, N)`` layout.

        Parameters
        ----------
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`
            The unloaded dark cube — frames live on disk and are
            fetched on demand via ``nirDark.getReadArray``.
        nreads : `int`, optional
            Number of dark frames to return. If None, return all from ``r0``.
        r0 : `int`
            0-indexed offset into the dark cube. Returns ``dark[r0:r0+nreads]``.
            Used when processing a sub-range of the data ramp; the dark slice
            then aligns with the data reads being processed.

        Returns
        -------
        `np.ndarray`
            The dark cube to subtract, shape ``(H, W, nreads)`` — the
            time axis is last to match the H4 ISR cube convention. The
            underlying ``ImageCube`` stores frames per-read in
            ``(N, H, W)`` form, so the transpose happens here at the
            boundary.
        """

        if r0 > 0:
            end = nirDark.nreads if nreads is None else r0 + nreads
            full = nirDark.getImageCube(nreads=end)
            cube = full[r0:end].copy()
        else:
            cube = nirDark.getImageCube(nreads=nreads)
        # If gain was applied, back it out. We used to apply darks to the
        # 2-d image in e-, but have switch to applying it to the raw rampin ADU.
        gain = self._darkGain(nirDark)
        if gain != 1.0:
            cube /= gain
        return np.ascontiguousarray(cube.transpose(1, 2, 0))

    def getDarkRead(self, nirDark, readNum) -> np.ndarray:
        """Get one read of the dark cube, in ADU.

        The dark is stored in electrons, so the gain it was taken with is backed
        out to match the raw ramp it is subtracted from.

        Parameters
        ----------
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`
            The unloaded dark cube — frames live on disk and are
            fetched on demand via ``nirDark.getReadArray``.
        readNum : `int`
            The read number to get.

        Returns
        -------
        `np.ndarray`
            The dark read. Always a new array: `ImageCube.__getitem__` caches the
            image it returns, so scaling it in place would corrupt the cube for
            every later read of the same index.
        """

        # Back out the applied gain (electrons -> ADU). Divide into a fresh
        # array: ImageCube caches the read, so an in-place /= would corrupt it.
        gain = self._darkGain(nirDark)
        return nirDark[readNum].array / gain

    def rampParams(self, pfsRaw):
        """Ramp-processing parameters, dispatched by observation type.

        Per-obstype defaults: arcs (``comparison``) and flats → full-ramp
        CDS; darks → full-ramp UTR; science → UTR over reads ``1:-3``
        (drop the shutter-closed read 0 and the trailing 3 shutter-closed
        / transitional reads). Anything else falls back to full-ramp UTR.

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
        if obsType in ("comparison", "flat"):
            # Arcs and flats use CDS. We may want to process longer
            # flats UTR — deferred until we have a principled
            # read-range choice.
            quickCDS, firstRead, lastRead = True, 0, -1
        elif obsType == "science":
            # Shutter open/close FITS cards are not yet written, so the
            # exact illuminated read range is unknown. Once those cards
            # are available the range can be computed correctly; until
            # then use reads 1:-3, which empirically drops the shutter-
            # closed read 0 and the trailing 3 shutter-closed reads.
            quickCDS, firstRead, lastRead = False, 1, -4
        else:  # dark, unknown — full-ramp UTR
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

    _INTERMEDIATE_KEYS = (
        'raw', 'darkSubbed', 'linearized', 'crCorrected', 'crResult',
        'rateStabilityResult',
    )

    def checkNirDark(self, pfsRaw, nirDark, nReadsNeeded):
        """Verify the NIR dark cube matches the ramp before subtraction.

        Two independent hard requirements, each raised as a fatal error
        rather than surfacing as a cryptic ``KeyError: 'Extension
        IMAGE_<n> not found'`` deep in :meth:`subtractDarkCube`:

        1. The dark's IRP ratio (``W_H4IRPN``) must equal the ramp's.
           Ramps taken at different data-to-IRP ratios have different
           per-read cadence, so a ratio-mismatched dark's reads do not
           correspond to the ramp's — the usual symptom of the
           ratio-matched ``nirDark_irp4`` not having been resolved. A dark
           with no ``W_H4IRPN`` is treated as irp1 (old darks predate the
           card), so it is still accepted for an irp1 ramp.
        2. The dark must have at least ``nReadsNeeded`` reads, enough to
           cover every read index the subtraction accesses
           (``0 .. nReadsNeeded - 1``).

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
            The raw H4 ramp being reduced.
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`
            The dark cube about to be subtracted.
        nReadsNeeded : `int`
            One past the highest absolute read index the subtraction will
            access (i.e. ``r1``).

        Raises
        ------
        RuntimeError
            If the IRP ratios differ, or the dark has too few reads.
        """
        rampRatio = pfsRaw.irpN
        # Old darks predate the W_H4IRPN card and are implicitly irp1.
        darkRatio = nirDark.metadata.get("W_H4IRPN")
        if darkRatio is None:
            darkRatio = 1
        if darkRatio != rampRatio:
            raise RuntimeError(
                f"nirDark IRP ratio ({darkRatio}) does not match ramp IRP "
                f"ratio ({rampRatio}); the ratio-matched dark was not resolved."
            )
        numReads = nirDark.getNumReads()
        if numReads < nReadsNeeded:
            raise RuntimeError(
                f"nirDark has only {numReads} reads, too few to cover the "
                f"ramp's {nReadsNeeded} reads."
            )

    def _maskDefects(self, exposure, defects):
        """Mask defects as BAD; interpolate over them only if configured.

        The base `~lsst.ip.isr.IsrTask.maskAndInterpolateDefects` always
        interpolates. This wrapper honors ``config.doInterpolate``: when it is
        False the defects are still BAD-masked but the pixel values are left
        untouched (no INTRP), so a "mask but don't interpolate" mode is possible
        on the H4 NIR path.
        """
        if self.config.doInterpolate:
            super().maskAndInterpolateDefects(exposure, defects)
        else:
            super().maskDefect(exposure, defects)

    def makeNirExposure(self,
                        pfsRaw,
                        nirDark=None,
                        linearity=None,
                        defects=None,
                        doReturnRawCube=False,
                        intermediates=None):
        """Build a 2-D image from the H4 ramp.

        The read range and CDS-vs-UTR choice come from
        :meth:`rampParams` (driven by ``config.h4`` and the observation
        type). The CDS arm linearizes and differences two endpoint
        frames; the UTR arm reads the cumulative ramp, subtracts the
        dark cube, optionally linearizes the whole cube, and
        — when ``config.h4.doCR`` is enabled — runs the iterative
        UTR-rate CR / ASIC-glitch detector on the per-read deltas.
        The returned exposure carries ``CR``, ``BAD``, ``SAT``,
        ``DARK_DEFECT``, ``LINEARITY_DEFECT``, ``UNSTABLE``, and
        ``RATE_UNSTABLE`` mask planes as appropriate. ``RATE_UNSTABLE``
        is always populated when the gate runs; whether those pixels
        also feed ``UNSTABLE`` + ``BAD`` is controlled by
        ``config.h4.maskRateUnstable``.

        Parameters
        ----------
        pfsRaw : `lsst.obs.pfs.PfsRaw`
            Raw H4 exposure giving access to the ramp.
        nirDark : `lsst.obs.pfs.imageCube.ImageCube`, optional
            Per-read dark cube to subtract from the ramp. Required when
            ``config.doDark`` is True.
        linearity : `lsst.obs.pfs.h4Linearity.H4Linearity`, optional
            Linearity solution. Required when ``config.h4.doLinearize``
            is True.
        defects : `lsst.ip.isr.Defects`, optional
            Known bad pixels; excluded from the CR/glitch detector's
            statistics and folded into the linearity-time mask.
        doReturnRawCube : `bool`, optional
            Reconstruct and return the (post-linearization, post-CR)
            cumulative ramp cube alongside the 2-D image. Defaults to
            False because the cube is large (~6.7 GB on a 4096²×100
            ramp) and the production path is delta-only.
        intermediates : `dict`, optional
            If provided, the keys *already present* select which stages
            to snapshot; each captured cube costs one full ``.copy()``.
            Recognized keys:

              ``'raw'``         pre-dark cumulative ramp
              ``'darkSubbed'``  ramp after dark subtraction (input to
                                linearity)
              ``'linearized'``  post-linearization, pre-CR cumulative
                                ramp
              ``'crCorrected'`` post-CR cumulative ramp (equals
                                ``'linearized'`` when
                                ``config.h4.doCR`` is False)
              ``'crResult'``    the iterative CR detector's
                                ``IterativeRepairResult`` — per-pixel
                                ``rate`` / ``sigma``, per-delta
                                ``crFlagMask`` / ``glitchFlagMask`` /
                                ``unclassifiedFlagMask``, per-pixel
                                ``badPixelMask`` (UNSTABLE / RTS),
                                iteration counts. None when
                                ``config.h4.doCR`` is False.
              ``'rateStabilityResult'``
                                the rate-stability detector's
                                ``RateStabilityResult`` — per-pixel
                                ``rejectMask`` / ``fraction`` /
                                ``nTestable`` / ``segmentRates``, scalar
                                ``nRejected`` / ``nUntestable``. None
                                when ``config.h4.doRateStability`` is
                                False or the path is not the linearized
                                UTR arm.

            Pass ``{'linearized': None}`` to capture just one stage;
            pass ``dict.fromkeys(['raw', 'darkSubbed', 'linearized',
            'crCorrected', 'crResult', 'rateStabilityResult'])`` for
            everything. The values are overwritten in place.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            The 2-D image with mask planes stamped.
        rawCube : `numpy.ndarray` or `None`
            The post-linearization, post-CR cumulative ramp when
            ``doReturnRawCube`` is True; otherwise ``None``. This is
            the same data the rate computation operates on — the
            production path keeps it in delta form, but the cube is
            rebuilt via ``cumsum`` for callers that need the
            ``(N, H, W)`` cumulative layout (e.g. for persistence to
            ``rawISRCube``).
        """

        # Capture set is whichever recognized keys the caller pre-seeded
        # in ``intermediates``. Skipping unrequested stages avoids their
        # full-cube ``.copy()`` overhead.
        if intermediates is None:
            captureKeys = frozenset()
        else:
            captureKeys = frozenset(intermediates)
            unknown = captureKeys - frozenset(self._INTERMEDIATE_KEYS)
            if unknown:
                raise ValueError(
                    f"unknown intermediates keys: {sorted(unknown)}; "
                    f"expected subset of {self._INTERMEDIATE_KEYS}."
                )

        # Dispatch CDS/UTR and the read range by observation type, then
        # resolve to absolute, 0-indexed inclusive bounds.
        quickCDS, firstRead, lastRead = self.rampParams(pfsRaw)
        if quickCDS and self.config.h4.applyUTRWeights:
            self.log.warn(
                "applyUTRWeights=True but obstype logic switched that to CDS; "
                "UTR weighting will not be applied for this exposure."
            )
        r0 = pfsRaw.positiveIndex(firstRead)
        r1 = pfsRaw.positiveIndex(lastRead)
        if r1 - r0 < 1:
            raise ValueError(
                f"firstRead={firstRead} (->{r0}) and "
                f"lastRead={lastRead} (->{r1}) leave no readable range "
                f"(need r1 > r0)."
            )
        if nirDark is not None:
            self.checkNirDark(pfsRaw, nirDark, r1)
        # Seed the H4 internal mask once, before the CDS-vs-UTR
        # dispatch — BORDER + DARK_DEFECT (input calib) + the linearity
        # calib's fit-time bits are the universal first step, regardless
        # of reduction path. Downstream stages (apply, the CR detector,
        # the projection at the end) read this mask and OR in their own
        # findings under the "first-reason-wins" rule.
        detBBox = pfsRaw.detector.getBBox()
        internalMask = _makeInternalMask(
            (detBBox.getHeight(), detBBox.getWidth()),
            linearity=linearity, defects=defects,
        )

        # Tracks whether ``exposure.image`` ends up UTR-weighted. The
        # quickCDS branch and the non-linearized non-UTR fallback are
        # the only paths that produce a CDS-style image; everything
        # else applies UTR weights (the linearized arm via the closed-
        # form delta weights, the legacy carve-out via
        # ``calcUTRrates``). Set per-branch below and stamped into
        # ``H4UTRWT`` by ``_stampRampMetadata``.
        appliedUTR = False
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
            if self.config.h4.doLinearize and linearity is not None:
                self.log.info("Correcting non-linearity (two-read CDS).")
                linR1, _ = h4Linearity.applyFrame(linearity, frameR1)
                if frameR0 is not None:
                    linR0, _ = h4Linearity.applyFrame(linearity, frameR0)
                    nirImage = linR1 - linR0
                else:
                    nirImage = linR1
                # Compute runtime range bits on the endpoint frames,
                # gating on the "first-reason-wins" rule: only set
                # ABOVE/BELOW range bits where the internal mask is
                # otherwise empty.
                goodPixels = internalMask == 0
                frames = [frameR1] if frameR0 is None else [frameR0, frameR1]
                for frame in frames:
                    internalMask[(frame > linearity.fitMax)
                                 & goodPixels] |= h4Linearity.ABOVE_VALID_RANGE
                    internalMask[(frame < linearity.fitMin)
                                 & goodPixels] |= h4Linearity.BELOW_VALID_RANGE
            else:
                nirImage = frameR1 if frameR0 is None else frameR1 - frameR0
        else:
            self.log.info(f"reading ramp over reads [{r0}, {r1}]...")
            # ``flux`` is the cumulative IRP-corrected ramp zero-anchored
            # at r0, shape ``(N, H, W)`` — frames first so the per-read
            # ``subtractDarkCube`` below writes contiguously. The single
            # transpose to ``(H, W, N)`` happens at the ``apply()``
            # boundary further down (in the linearization arm) or stays
            # ``(N, H, W)`` all the way through the no-linearization
            # fallback.
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

            if 'raw' in captureKeys:
                # PixelRampData / the rest of the diagnostic chain wants
                # the captured cubes in (H, W, N) form to match
                # cubeLin / cubeCR. Transpose at capture time — paid
                # only on the diagnostic path.
                intermediates['raw'] = np.ascontiguousarray(flux.transpose(1, 2, 0))

            if nirDark is not None:
                self.log.info("subtracting dark cube (per-read, in place).")
                self.subtractDarkCube(nirDark, flux, r0=r0)

            if 'darkSubbed' in captureKeys:
                intermediates['darkSubbed'] = np.ascontiguousarray(
                    flux.transpose(1, 2, 0)
                )

            if self.config.h4.doLinearize and linearity is not None:
                self.log.info("Correcting non-linearity.")
                # h4Linearity.apply, the CR detector, the diff and the
                # cumsum reconstruction all want the time axis last and
                # contiguous — per-pixel Horner / partition / IQR /
                # cumsum then stride 1 along reads instead of crossing
                # the slowest axis. Single transpose here; everything
                # downstream of this point stays ``(H, W, N)``.
                flux = np.ascontiguousarray(flux.transpose(1, 2, 0))
                # Hand the pre-seeded internal mask to apply() as
                # ``validMask``; apply() OR's its own findings (range
                # bits) back in under the first-reason-wins rule.
                ramp = h4Linearity.Ramp(reads=flux, validMask=internalMask)
                linearizedRamp = h4Linearity.apply(linearity, ramp)
                flux = linearizedRamp.cumulativeLinear
                internalMask = linearizedRamp.badPixelMask
                self.log.info(
                    f"   starting with {int((internalMask != 0).sum())} "
                    f"already-flagged pixels (border + defects + fit-time bits)"
                )
                # The pre-linearization cube is now superseded by the
                # linearized output. Drop the Ramp / LinearizedRamp
                # holders so the input (~6.7 GB) can be GC'd before the
                # CR/glitch detector copies a fresh ``cubeOriginal``.
                del ramp, linearizedRamp

                # Re-anchor the linearized cube at r0: subtract the linearized value
                # at read r0 from every read so the output represents only the flux
                # accumulated during (r0, r1]. Without this rebase, postISRCCD over
                # disjoint sub-ranges would not sum to the full-ramp postISRCCD.
                # No-op for r0=0 because the linearity model maps 0 -> 0.
                if r0 > 0 and offsetRaw is not None:
                    self.log.info(f"re-anchoring linearized cube at read {r0}.")
                    if nirDark is not None:
                        # Linearization saw flux already dark-subtracted; the rebase
                        # offset must match. offsetRaw is the absolute cumulative
                        # at read r0, so its dark is nirDark[r0-1] (nirDark[j] is
                        # the dark for read j+1 — the same k->k-1 indexing
                        # getDarkCube uses).
                        offsetForLin = offsetRaw - self.getDarkRead(nirDark, r0 - 1)
                    else:
                        offsetForLin = offsetRaw
                    linearizedOffset, _ = h4Linearity.applyFrame(linearity, offsetForLin)
                    flux -= linearizedOffset[..., None]

                if 'linearized' in captureKeys:
                    intermediates['linearized'] = flux.copy()  # (H, W, N) post-lin, pre-CR

                # Switch from cumulative flux to delta-space for everything
                # downstream of linearization. The CR/glitch detector
                # operates on deltas directly (no diff/cumsum dance inside).
                # We only reconstruct the cumulative cube on demand
                # (intermediates / doReturnRawCube).
                read0 = flux[..., 0:1].copy()
                deltas = np.diff(flux, axis=-1)
                del flux
                iterResult = None  # filled in if doCR runs
                rsResult = None    # filled in if doRateStability runs

                if self.config.h4.doCR:
                    # First-reason-wins: CR detector only runs on
                    # pixels with no pre-existing reason in the
                    # internal mask.
                    crGood = (internalMask == 0)

                    if self.config.h4.doDeglitch:
                        self.log.info("Correcting CRs and ASIC glitches.")
                        # Glitch detection runs on all pixels. The matched-pair
                        # cancellation criterion (opposite signs, sum within
                        # threshold) is what discriminates a glitch from a CR;
                        # restricting by ASIC channel misses glitches on
                        # channels other than the historically-noted ones.
                        glitchChanMask = np.ones(deltas.shape[:-1], dtype=bool)
                    else:
                        self.log.info("Correcting CRs (ASIC deglitching disabled).")
                        glitchChanMask = None
                    # deltas is already (H, W, N-1) — pass directly to
                    # the CR detector; no boundary transpose here.
                    iterResult = h4Linearity.cr.iterativeUtrDetectAndRepair(
                        deltas,
                        goodPixelMask=crGood,
                        glitchPixelMask=glitchChanMask,
                        sigmaFloorADU=self.config.h4.rateCRsigmaFloorADU,
                        nSigma=self.config.h4.rateCRnSigma,
                        maxIterations=self.config.h4.rateCRiterMax,
                        repair=self.config.h4.repairCR,
                        correctGlitches=self.config.h4.correctGlitches,
                        glitchAmplitudeMinADU=self.config.h4.deglitchAmplitudeMinADU,
                        maxDropFraction=self.config.h4.rateCRmaxDropFraction,
                        nDropSigma=self.config.h4.rateCRnDropSigma,
                        badPixelMinOutliers=self.config.h4.badPixelMinOutliers,
                        badPixelOutlierSigma=self.config.h4.badPixelOutlierSigma,
                    )
                    # Suppress CR / glitch flags at pixels classified
                    # BAD — RTS pixels often paint up as one or two CRs
                    # by accident; downstream consumers should treat the
                    # whole pixel as bad, not attempt CR repair.
                    badPix2D = iterResult.badPixelMask
                    if badPix2D.any():
                        iterResult.crFlagMask[badPix2D] = False
                        iterResult.glitchFlagMask[badPix2D] = False
                    crFlagMask2D = iterResult.crFlagMask.any(axis=-1)
                    glitchFlagMask2D = iterResult.glitchFlagMask.any(axis=-1)
                    unclassMask2D = iterResult.unclassifiedFlagMask.any(axis=-1)
                    crResult = SimpleNamespace(
                        flagMask=crFlagMask2D,
                        nFlagged=int(crFlagMask2D.sum()),
                        glitchFlagMask=glitchFlagMask2D,
                        nGlitchFlagged=int(glitchFlagMask2D.sum()),
                        badPixelMask=badPix2D,
                        nBadPixels=int(badPix2D.sum()),
                    )
                    # OR CR-stage findings back into the internal mask.
                    # crResult.flagMask (CR-only) stays out of the
                    # internal mask — projected directly to the CR
                    # plane by _projectInternalMask.
                    internalMask[badPix2D] |= h4Linearity.UNSTABLE
                    internalMask[unclassMask2D] |= h4Linearity.UNCLASSIFIED
                    # Flag every ASIC-glitch pixel (pair height above
                    # ``asicGlitchHeightMaskADU``, default 0) with the
                    # published, non-BAD ASIC_GLITCH plane. Heights come
                    # straight from the saved deltas because interior
                    # pairs are not repaired (correctGlitches=False by
                    # default). Clean glitches in correction-scoped
                    # channels are kept -- result.rate is the mean of the
                    # un-flagged deltas, so a clean matched pair already
                    # cancels there and the flux is unbiased. Messy or
                    # out-of-scope glitch pixels additionally get
                    # GLITCH_MASKED, which promotes them to BAD.
                    glitchMask3D = iterResult.glitchFlagMask
                    if glitchMask3D.any():
                        pairStart = (glitchMask3D[..., :-1]
                                     & glitchMask3D[..., 1:])
                        halfDiff = 0.5 * (deltas[..., :-1]
                                          - deltas[..., 1:])
                        maxHeight = np.where(
                            pairStart, np.abs(halfDiff), 0.0,
                        ).max(axis=-1)
                        glitchPix = (
                            (maxHeight
                             > self.config.h4.asicGlitchHeightMaskADU)
                            & glitchMask3D.any(axis=-1)
                        )
                        correctable = self.correctableGlitchMask(
                            pfsRaw.detector.getName(), glitchMask3D,
                            iterResult.crFlagMask, deltas, badPix2D,
                        )
                        internalMask[glitchPix] |= h4Linearity.ASIC_GLITCH
                        masked = glitchPix & ~correctable
                        if masked.any():
                            internalMask[masked] |= h4Linearity.GLITCH_MASKED
                    if 'crResult' in captureKeys:
                        intermediates['crResult'] = iterResult
                    self.log.info(
                        f"iterative CR/glitch step: "
                        f"{iterResult.nCRs} CR flag entries, "
                        f"{iterResult.nGlitchPairs} glitch pairs, "
                        f"{crResult.nFlagged} unique CR pixels, "
                        f"{crResult.nGlitchFlagged} unique glitch pixels, "
                        f"{crResult.nBadPixels} BAD (RTS) pixels "
                        f"in {iterResult.nIterations} iterations."
                    )
                else:
                    crResult = None

                if self.config.h4.doRateStability:
                    # First-reason-wins: rate-stability tests only
                    # pixels with no pre-existing internal-mask reason
                    # (BORDER / DEFECT / DEAD / range / CR-stage
                    # UNSTABLE / HIGH_FIT_RESIDUAL / ASIC_GLITCH).
                    # CR- and glitch-flagged deltas are excluded from
                    # each half's mean rate.
                    rsGood = (internalMask == 0)
                    if iterResult is not None:
                        rsFlagMask = (iterResult.crFlagMask
                                      | iterResult.glitchFlagMask)
                    else:
                        rsFlagMask = np.zeros(deltas.shape, dtype=bool)
                    rsResult = h4Linearity.rateStability.detectRateInstability(
                        deltas, rsFlagMask,
                        goodPixelMask=rsGood,
                        threshold=self.config.h4.rateStabilityThreshold,
                        rateFloorADU=self.config.h4.rateStabilityRateFloorADU,
                        minDeltasPerSegment=(
                            self.config.h4.rateStabilityMinDeltasPerSegment
                        ),
                    )
                    if rsResult.nRejected:
                        internalMask[rsResult.rejectMask] |= h4Linearity.RATE_UNSTABLE
                    if 'rateStabilityResult' in captureKeys:
                        intermediates['rateStabilityResult'] = rsResult
                    self.log.info(
                        f"rate-stability step: {rsResult.nRejected} "
                        f"pixels rejected, {rsResult.nUntestable} untestable."
                    )

                # Reconstruct the cumulative cube ONLY when a debug
                # capture (``'crCorrected'``) or ``doReturnRawCube`` asks
                # for it. The default production path is delta-only from
                # here on. cumsum along the contiguous time axis is
                # ~6 s vs ~40 s along the strided axis on a 4096²×100
                # cube.
                flux = None
                if ('crCorrected' in captureKeys) or doReturnRawCube:
                    fluxHWN = np.empty(
                        deltas.shape[:-1] + (deltas.shape[-1] + 1,),
                        dtype=deltas.dtype,
                    )
                    fluxHWN[..., 0:1] = read0
                    np.cumsum(deltas, axis=-1, out=fluxHWN[..., 1:])
                    fluxHWN[..., 1:] += read0
                    if 'crCorrected' in captureKeys:
                        intermediates['crCorrected'] = fluxHWN.copy()  # (H, W, N)
                    if doReturnRawCube:
                        # ``doReturnRawCube`` is consumed by ``runH4RG``
                        # → ``ImageCube.fromCube``, which iterates the
                        # cube frame-by-frame and so still needs the
                        # ``(N, H, W)`` layout. Transpose only here.
                        flux = np.ascontiguousarray(fluxHWN.transpose(2, 0, 1))
                    del fluxHWN
            else:
                # No linearization path: the pre-seeded internalMask
                # (BORDER + DARK_DEFECT only — linearity bits were skipped
                # since no calib was loaded) carries through unchanged.
                crResult = None
                # No linearization → no apply-boundary transpose ran;
                # flux is still ``(N, H, W)`` from makeUTRcumulative,
                # which is the form the legacy ``calcUTRrates`` and
                # ``flux[-1] - flux[0]`` paths below expect.
                deltas = None
                read0 = None

            if deltas is not None:
                # Delta-space science image: UTR-weighted rate × nReads,
                # equivalent to ``calcUTRrates`` on the reconstructed
                # (CR-corrected) flux ramp ``read0 + cumsum(deltas)``.
                # The CR detector already computes this for us in
                # ``iterResult.rate`` (delta-space closed form of the
                # read-space UTR weights). When CR is disabled, apply
                # the same delta-space weights directly. ``deltas`` is
                # (H, W, N-1) here — the time axis is last.
                if self.config.h4.doCR and iterResult is not None:
                    pixelRate = iterResult.rate
                else:
                    pixelRate = self.calcUTRrateFromDeltas(deltas)
                nirImage = pixelRate * deltas.shape[-1]
                appliedUTR = True
                del deltas, read0
            elif self.config.h4.applyUTRWeights:
                # Carve-out for non-linearizable detectors: they have no
                # linearity solution but still need a UTR-weighted rate
                # path. The linearized arm above is the long-term home
                # for applyUTRWeights once a delta-form calcUTRrates
                # exists, but for cameras without a linearity curve this
                # branch stays load-bearing.
                self.log.info("applying UTR weights.")
                rates = self.calcUTRrates(flux)
                nirImage = rates * len(flux)
                appliedUTR = True
            else:
                # No linearization and no UTR weighting: two-read
                # difference on the already dark-subtracted ramp.
                nirImage = flux[-1] - flux[0]

        exposure = self._makeExposure(pfsRaw, nirImage)
        _stampRampMetadata(
            exposure, r0=r0, r1=r1,
            nTotal=int(pfsRaw.getNumReads()),
            appliedUTR=appliedUTR,
        )
        # Single projection point: lift the internal mask + CR result
        # into Exposure.mask planes (DARK_DEFECT, LINEARITY_DEFECT, SAT,
        # UNSTABLE, RATE_UNSTABLE, BAD, CR). See ``_projectInternalMask``.
        _projectInternalMask(
            exposure, internalMask, crResult=crResult,
            maskRateUnstable=self.config.h4.maskRateUnstable,
        )

        nCR = crResult.nFlagged if crResult is not None else 0
        self.log.info(
            f"nSat={int(((internalMask & h4Linearity.ABOVE_VALID_RANGE) != 0).sum())} "
            f"nLow={int(((internalMask & h4Linearity.BELOW_VALID_RANGE) != 0).sum())} "
            f"nDefects={int(((internalMask & h4Linearity.MASKED_BY_INPUT) != 0).sum())} "
            f"nLinDefects={int(((internalMask & h4Linearity.DEAD) != 0).sum())} "
            f"nUnstable={int(((internalMask & h4Linearity.UNSTABLE) != 0).sum())} "
            f"nRateUnstable={int(((internalMask & h4Linearity.RATE_UNSTABLE) != 0).sum())} "
            f"nUnclassified={int(((internalMask & h4Linearity.UNCLASSIFIED) != 0).sum())} "
            f"nGlitchMasked={int(((internalMask & h4Linearity.ASIC_GLITCH) != 0).sum())} "
            f"nCR={nCR}"
        )

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

    def getBadIRPpixels(self) -> np.ndarray:
        """Return the bad IRP row pixel list from the butler-provided calib."""

        calib = getattr(self, "_badRefPixels", None)
        if calib is None:
            self.log.warn('no badRefPixels calib available; skipping bad-IRP-pixel repair')
            return np.zeros(dtype=np.int16, shape=())

        return np.asarray(calib.pixels, dtype=np.int32)

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

        badPixels = self.getBadIRPpixels()

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

    def loadCorrectGlitchChannels(self, detectorName: str, nChannels: int = 32) -> tuple:
        """ASIC channels whose clean glitches are *corrected* (kept) rather than
        masked, on pixels with a clean glitch/CR pattern.

        Every channel is in scope except the known-bad ASIC channels
        (:meth:`loadBadAsicChannels`) -- e.g. n3 channel 24, whose IRP is
        bypassed under Hann filtering, injecting large digital glitches that
        can't be cleanly corrected, so it stays detect-and-mask. Under
        ``IRPfilter == -1`` the per-column median cleans those channels (their
        glitches become ordinary correctable outlier+return pairs), so nothing
        is excluded.
        """
        # The bad-ASIC exclusion only applies while the IRP is bypassed / raw
        # (IRPfilter != -1); the -1 median cleans those channels.
        bad = set() if self.config.h4.IRPfilter == -1 else set(self.loadBadAsicChannels(detectorName))
        return tuple(c for c in range(nChannels) if c not in bad)

    def correctableGlitchMask(self, detectorName: str, glitchFlagMask: np.ndarray,
                              crFlagMask: np.ndarray, deltas: np.ndarray, badMask: np.ndarray,
                              nChannels: int = 32, returnFraction: float = 0.5) -> np.ndarray:
        """Return a ``(H, W)`` bool mask True at pixels eligible for glitch
        *correction* rather than masking.

        An ASIC glitch is a two-delta event: a one-read outlier followed by a
        roughly-equal, opposite-sign return. A pixel is correctable when it is
        in a correction-scoped channel (:meth:`loadCorrectGlitchChannels`), not
        already flagged bad, and every glitch behaves that way:

        - **outlier** -- a glitch-flagged delta whose residual from the pixel's
          robust rate exceeds the detector's own threshold,
          ``rateCRnSigma * max(sigma, rateCRsigmaFloorADU)`` (``sigma`` the IQR
          estimate). This is the same test that flagged it.
        - **return** -- every outlier must have an adjacent glitch-flagged delta
          of opposite sign and magnitude at least ``returnFraction`` of the
          outlier. A lone outlier with no such return (a persistent step, i.e.
          CR-like) makes the pixel non-correctable.
        - Sub-threshold glitch flags (the small strays the matched-pair logic
          sweeps in around a real pair) are **not** outliers, so they need no
          return and never by themselves disqualify a pixel.
        - a CR (dropped from the rate separately) must have **no flagged
          neighbour**, so glitch-vs-CR is unambiguous.

        This keeps a clean pair, any number of well-separated pairs, and a pair
        carrying a small spurious third flag (the common "triple-run"), while
        masking lone persistent spikes and entangled multi-outlier structures.
        CR detection and repair are untouched -- the CR flags are read only for
        the adjacency guard.

        Parameters
        ----------
        detectorName : `str`
            E.g. ``"n4"``.
        glitchFlagMask, crFlagMask : `np.ndarray`
            ``(H, W, nDeltas)`` bool per-read glitch / CR flags.
        deltas : `np.ndarray`
            ``(H, W, nDeltas)`` per-read deltas of the linearized ramp; used to
            find the outlier and test that its return roughly cancels it.
        badMask : `np.ndarray`
            ``(H, W)`` bool; pixels already classified bad (RTS/UNSTABLE).
        nChannels : `int`
            Number of horizontal ASIC channels stacked along Y; 32 for H4.
        returnFraction : `float`
            An outlier's return must be at least this fraction of its magnitude
            (opposite sign) to count. ~0.5 cleanly separates the real ~equal
            return from a small spurious stray.
        """
        channels = self.loadCorrectGlitchChannels(detectorName)
        H, W = glitchFlagMask.shape[:2]
        glitch = glitchFlagMask.astype(bool)

        inScope = np.zeros((H, W), dtype=bool)
        channelHeight = H // nChannels
        for ch in channels:
            inScope[ch * channelHeight:(ch + 1) * channelHeight, :] = True

        # Restrict the per-delta work to candidate pixels (in scope, not bad,
        # at least one glitch flag) -- on the full detector this is a tiny
        # fraction, so the IQR/percentile pass stays cheap and memory-light.
        candidate = inScope & ~badMask.astype(bool) & glitch.any(-1)
        out = np.zeros((H, W), dtype=bool)
        ys, xs = np.where(candidate)
        if len(ys) == 0:
            return out

        g = glitch[ys, xs]                             # (Ncand, nDeltas)
        c = crFlagMask.astype(bool)[ys, xs]
        d = np.asarray(deltas)[ys, xs].astype(np.float64)
        absd = np.abs(d)
        sgn = np.sign(d)

        # Detector-equivalent outlier test: |delta - rate| above nSigma*sigma.
        p25 = np.percentile(d, 25, axis=-1)
        p75 = np.percentile(d, 75, axis=-1)
        rate = np.median(d, axis=-1)
        sigma = (p75 - p25) / 1.349
        thresh = self.config.h4.rateCRnSigma * np.maximum(
            sigma, self.config.h4.rateCRsigmaFloorADU)
        outlier = g & (np.abs(d - rate[:, None]) > thresh[:, None])

        # Each outlier needs an adjacent glitch delta, opposite sign, magnitude
        # >= returnFraction * |outlier| -- its return.
        leftReturn = np.zeros_like(g)
        leftReturn[:, 1:] = (g[:, :-1] & (sgn[:, :-1] * sgn[:, 1:] < 0)
                             & (absd[:, :-1] >= returnFraction * absd[:, 1:]))
        rightReturn = np.zeros_like(g)
        rightReturn[:, :-1] = (g[:, 1:] & (sgn[:, 1:] * sgn[:, :-1] < 0)
                               & (absd[:, 1:] >= returnFraction * absd[:, :-1]))
        hasReturn = leftReturn | rightReturn
        unpairedOutlier = (outlier & ~hasReturn).any(-1)

        # A CR must not sit next to any flagged delta.
        flagged = g | c
        flaggedNeighbor = np.zeros_like(flagged)
        flaggedNeighbor[:, :-1] |= flagged[:, 1:]
        flaggedNeighbor[:, 1:] |= flagged[:, :-1]
        crAdjacent = (c & flaggedNeighbor).any(-1)

        out[ys, xs] = outlier.any(-1) & ~unpairedOutlier & ~crAdjacent
        return out

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
        """Return a diff IRP image: each channel column replaced by the median
        across its reference rows.

        For IRP1, bad reference rows (the ``badRefPixels`` calib) are excluded
        from the per-column median so they cannot bias it. The calib is indexed
        in the IRP1 frame, so the masking is applied only when ``irpN == 1``;
        for IRP4 (any ``irpN > 1``) the median runs over all rows -- bad-ref
        masking there needs an IRP4-aware calib mapping and is deferred.
        """
        nchan = pfsRaw.nchan
        h, w = rawDiffIrp.shape
        chan_w = h//nchan

        # badRefPixels is indexed in the IRP1 frame, so only mask for irpN == 1.
        badRows = np.zeros(0, dtype=np.int32)
        calib = getattr(self, "_badRefPixels", None)
        if calib is not None and self.config.h4.doIRPbadPixels and pfsRaw.irpN == 1:
            badRows = np.asarray(calib.pixels, dtype=np.int32)

        out = np.zeros_like(rawDiffIrp)
        for chan_i in range(nchan):
            rowLow = chan_i*chan_w
            rowHigh = (chan_i + 1)*chan_w
            chan0 = rawDiffIrp[rowLow:rowHigh, :]

            chanBad = badRows[(badRows >= rowLow) & (badRows < rowHigh)] - rowLow
            if len(chanBad) > 0:
                good = np.ones(chan_w, dtype=bool)
                good[chanBad] = False
                chanVec = np.nanmedian(chan0[good, :], axis=0, keepdims=True)
            else:
                chanVec = np.nanmedian(chan0, axis=0, keepdims=True)
            out[rowLow:rowHigh, :] = chanVec  # broadcasts (1, w) down the channel's rows
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

        # The reference-pixel filter runs on the DIFF IRP (irpN - irp0), never
        # on a single IRP read -- and that is essential. The IRP pixels carry a
        # strong, fixed pixel-to-pixel capacitance pattern; on a single read
        # that pattern dominates the rank ordering of a channel's 128 reference
        # rows, so the per-column median degenerates to repeatedly selecting the
        # same fixed-pattern pixel (measured n4/ch18: one row is the median for
        # 92% of columns; only ~10-12 of 128 rows are ever chosen). Differencing
        # against read0 cancels the fixed pattern, leaving smooth noise, so the
        # median draws ~uniformly across all 128 rows and is a genuine robust
        # estimator.
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
           Cumulative IRP-corrected ADU, zero-anchored at read ``r0``,
           shape ``(nreads-1, H, W)``. ``flux[k]`` is read ``r0+k+1``
           minus read ``r0`` minus the per-channel IRP-filtered diff.
           The frame-first layout is what the production caller wants
           for the in-place ``subtractDarkCube`` step; a single
           ``ascontiguousarray(flux.transpose(1, 2, 0))`` at the
           ``apply()`` boundary is the only transpose between here and
           the linearization step.
        """

        r0 = pfsRaw.positiveIndex(r0)
        r1 = pfsRaw.positiveIndex(r1)
        if r1-r0 < 1:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}) + 1:')
        if nreads is None:
            nreads = r1 - r0 + 1
        reads = np.linspace(r0, r1, nreads, dtype='i2')

        # Reference-correct via interleaved reference pixels when available,
        # otherwise via the Teledyne refPixel4 border pixels -- matching makeCDS.
        useIrp = self.config.h4.useIRP and pfsRaw.irpN > 0

        # Grab the components of read r0, which we will subtract from all the others.
        data0 = self.makeRawDataArray(pfsRaw, r0)
        if useIrp:
            irp0 = self.makeRawIrpArray(pfsRaw, r0)
            self.applyIRPcrosstalk(pfsRaw, irp0, data0)
        else:
            data0 = self.borderCorrect(pfsRaw, data0)

        # We are not squirreling away the bbox, but really should for the final Exposure.
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
            if useIrp:
                irp1 = self.makeRawIrpArray(pfsRaw, r_i)
                t1 = time.time()
                self.applyIRPcrosstalk(pfsRaw, irp1, data1)
                ddata = data1 - data0
                ddata -= self.getFinalDiffIrp(pfsRaw, irp1 - irp0)
            else:
                t1 = time.time()
                ddata = self.borderCorrect(pfsRaw, data1) - data0
            if bbox is None:
                stack[r_idx-1, :, :] = ddata
            else:
                stack[r_idx-1, :, :] = ddata[bbox.getBeginY():bbox.getEndY(),
                                             bbox.getBeginX():bbox.getEndX()]
            t2 = time.time()
            if showTimes:
                print(f'cds {r_i} io1={t1-t0:0.3f} proc={t2-t1:0.3f}')

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

    def _getMetadataName(self):
        return None                     # don't write metadata; requires fix to ip_isr
