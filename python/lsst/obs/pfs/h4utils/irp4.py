"""Determine the IRP bad-pixel alignment (offset) for IRP4 data from real darks.

H4RG detectors can read out interleaved reference pixels (IRP) at a 1:N data:IRP
ratio. In IRP1 every reference pixel is read; in IRP4 only one of every four is
read, at phase ``W_H4IRPO`` within each group. The known-bad reference pixels in
``drp_pfs_data/h4/badRefPixels.yaml`` were measured from IRP1 darks, in IRP1
full-frame coordinates (0..4095 along the channel-stacked axis). To reuse that
list for IRP4 we need the phase (``offset``) that maps a subsampled IRP sample
``s`` back to its IRP1 index ``p = s*irpN + offset``.

The bad pixels are found with a quick CDS test: difference the first two IRP1
reads and flag reference rows whose robust scatter across the width is large::

    cube = task.makeRawIrpNcube(raw, r0=0, r1=1)
    irp1 = cube[-1] - cube[0]
    rows = np.where(iqrStd(irp1, axis=1) > threshold)[0]

In IRP4 a single bad reference pixel is repeat-expanded, so it shows up as
``irpN`` adjacent noisy rows. The alignment test then picks the ``offset`` whose
correction (the known IRP1 bad list, mapped to samples and expanded) best covers
those observed noisy rows. Meant to be imported and run from a notebook;
``main()`` is a thin command-line wrapper.

Worked example (n1, an IRP4 dark): rows 3896-3899 read as strikingly noisy.
Only 3899 is in the IRP1 bad list, and ``3899 = 974*4 + 3``, so only ``offset=3``
maps it to sample 974 and repairs the whole 3896-3899 block.

Typical notebook use::

    from lsst.daf.butler import Butler
    from lsst.obs.pfs.h4utils import irp4

    butler = Butler("/repo", collections="PFS/raw/all")
    raw = butler.get("raw", dataId=dict(visit=144800, arm="n", spectrograph=1))
    task = irp4.makeQuickIsrTask()
    result = irp4.alignByCorrection(task, raw)
    print(result.report())
"""

import os
from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import yaml
from lsst.utils import getPackageDir

__all__ = [
    "iqrStd",
    "excessRowNoise",
    "makeQuickIsrTask",
    "loadBadRefPixels",
    "noisyIrpRows",
    "correctedRowsForOffset",
    "correctedRows",
    "collapseToIrp1Pixels",
    "offsetCoverage",
    "offsetCoverageByParity",
    "alignByCorrection",
    "validateCorrection",
    "CorrectionValidation",
    "rampExcess",
    "scanBadRefPixels",
    "scanBadRefPixelsAcrossVisits",
    "proposedBadRefPixelsYaml",
    "plotScan",
    "runBadRefPixelSurvey",
    "BadRefPixelScan",
    "IrpAlignment",
    "plotIrpRowNoise",
    "main",
]

H4_SIZE = 4096


def iqrStd(a: np.ndarray, axis=None) -> np.ndarray:
    """Robust standard deviation from the interquartile range (0.741 * IQR)."""
    q1, q3 = np.nanpercentile(a, [25, 75], axis=axis)
    return (q3 - q1) * 0.741


def _robustLinearBaseline(x: np.ndarray, y: np.ndarray,
                          nsig: float = 3.0, iters: int = 3) -> np.ndarray:
    """Sigma-clipped linear fit of ``y`` vs ``x``, evaluated at ``x``.

    Bad rows (which sit above the trend) are clipped from the fit so they stand
    out in the residual rather than tilting the baseline.
    """
    good = np.isfinite(y)
    if good.sum() < 2:
        return np.full_like(y, np.nanmedian(y) if good.any() else 0.0)
    coef = np.polyfit(x[good], y[good], 1)
    for _ in range(iters):
        resid = y - np.polyval(coef, x)
        med = np.nanmedian(resid[good])
        sig = 1.4826 * np.nanmedian(np.abs(resid[good] - med))
        if not np.isfinite(sig) or sig == 0:
            break
        newgood = np.isfinite(y) & (np.abs(resid - med) < nsig * sig)
        if newgood.sum() < 2 or newgood.sum() == good.sum():
            break
        good = newgood
        coef = np.polyfit(x[good], y[good], 1)
    return np.polyval(coef, x)


def excessRowNoise(std: np.ndarray, nchan: int = 32) -> np.ndarray:
    """Per-row noise above each channel's own baseline.

    Subtracts, from each row's ``iqrStd``, a robust per-channel baseline fitted
    as a *line* along the channel (a sigma-clipped linear fit). Much of the
    per-row scatter is a fixed property of the ASIC channel and varies between
    FPAs; some channels also ramp along their length, so a linear baseline (not
    just a constant pedestal) is needed to keep the ends of a sloped channel from
    being spuriously flagged. The residual is what a single threshold can flag
    consistently across channels and detectors.

    Parameters
    ----------
    std : `numpy.ndarray`
        Per-row metric (length ``H4_SIZE``), e.g. from ``iqrStd(cds, axis=1)``.
    nchan : `int`
        Number of readout channels.

    Returns
    -------
    excess : `numpy.ndarray`
        ``std`` minus the per-channel linear baseline, same length as ``std``.
    """
    chanHeight = len(std) // nchan
    x = np.arange(chanHeight)
    out = np.empty(len(std), dtype=float)
    for c in range(nchan):
        lo = c * chanHeight
        chan = std[lo:lo + chanHeight]
        out[lo:lo + chanHeight] = chan - _robustLinearBaseline(x, chan)
    return out


def makeQuickIsrTask(**configOverrides):
    """Build a `PfsIsrTask` configured for IRP reference-pixel inspection.

    Turns off the steps that would obscure the raw reference structure when
    eyeballing IRP corrections: dark subtraction, H4 linearization, and IRP
    smoothing, and uses quick CDS. The config is set before construction so the
    task's prerequisite inputs (dark, linearizer) are dropped accordingly.

    Parameters
    ----------
    **configOverrides
        Extra ``config`` attributes set as dotted paths, e.g.
        ``makeQuickIsrTask(**{"h4.IRPfilter": 15})``.

    Returns
    -------
    task : `lsst.obs.pfs.isrTask.PfsIsrTask`
    """
    from lsst.obs.pfs.isrTask import PfsIsrTask

    config = PfsIsrTask.ConfigClass()
    config.h4.quickCDS = True
    config.h4.IRPfilter = 0
    config.doDark = False
    config.h4.doLinearize = False
    for key, value in configOverrides.items():
        target = config
        *parents, leaf = key.split(".")
        for name in parents:
            target = getattr(target, name)
        setattr(target, leaf, value)
    return PfsIsrTask(config=config)


def loadBadRefPixels(detectorName: str) -> np.ndarray:
    """Load the IRP1-coordinate bad reference pixel list for a detector.

    Reads ``drp_pfs_data/h4/badRefPixels.yaml``. Returns a 1-d int array of
    indices along the channel-stacked axis (0..4095), or an empty array if the
    detector is absent.
    """
    absFilename = os.path.join(getPackageDir("drp_pfs_data"), "h4", "badRefPixels.yaml")
    with open(absFilename) as f:
        cfg = yaml.safe_load(f)
    if detectorName not in cfg or cfg[detectorName] is None:
        return np.zeros(0, dtype=int)
    return np.array(cfg[detectorName], dtype=int)


def noisyIrpRows(task, raw, threshold: float = 2.0,
                 r0: int = 0, r1: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """Flag noisy IRP1 rows from a CDS of the first reads, above channel baseline.

    Differences two IRP1 reads (default 0 and 1) and flags rows whose robust
    scatter across the width exceeds their channel's baseline by more than
    ``threshold``. The CDS cancels the large constant per-pixel reference
    offsets, so a bad reference pixel stands out as a high-scatter row; the
    per-channel baseline subtraction (see `excessRowNoise`) removes the fixed
    ASIC-channel pedestal, which varies between FPAs, so one threshold works
    across channels and detectors. In IRP4 a bad pixel is repeat-expanded into
    ``irpN`` adjacent noisy rows.

    Bad-pixel replacement is disabled for the measurement: otherwise the
    ``constructFullIrp`` correction would repair already-known bad pixels before
    the CDS and hide them.

    Parameters
    ----------
    task : `lsst.obs.pfs.isrTask.PfsIsrTask`
        A task whose ``makeRawIrpNcube`` builds the IRP1 cube.
    raw : `lsst.obs.pfs.PfsRaw`
        The raw ramp.
    threshold : `float`
        Rows whose iqrStd exceeds their channel-median iqrStd by more than this
        (ADU) are flagged.
    r0, r1 : `int`
        Reads to difference (default 0 and 1, i.e. CDS).

    Returns
    -------
    rows : `numpy.ndarray`
        Flagged IRP1 row indices (0..4095).
    excess : `numpy.ndarray`
        Per-row iqrStd above the channel baseline, length ``H4_SIZE``.
    """
    cfg = task.config.h4
    prev = cfg.doIRPbadPixels
    cfg.doIRPbadPixels = False
    try:
        cube = task.makeRawIrpNcube(raw, r0=r0, r1=r1)
    finally:
        cfg.doIRPbadPixels = prev
    excess = excessRowNoise(iqrStd(cube[-1] - cube[0], axis=1), raw.nchan)
    rows = np.where(excess > threshold)[0]
    return rows, excess


def correctedRowsForOffset(badPixels: np.ndarray, irpN: int, offset: int) -> np.ndarray:
    """IRP1 rows that the IRP bad-pixel correction repairs at a given offset.

    A known bad IRP1 pixel ``p`` is only sampled in IRPn when
    ``(p - offset) % irpN == 0``; that sample ``s = (p - offset) // irpN`` is
    repeat-expanded back across its whole group, so the repaired rows are
    ``[s*irpN, s*irpN + irpN)``. For ``irpN == 1`` this is just ``badPixels``.
    """
    if irpN <= 1:
        return np.unique(badPixels)
    relevant = badPixels[(badPixels - offset) % irpN == 0]
    samples = (relevant - offset) // irpN
    rows = (samples[:, None] * irpN + np.arange(irpN)).ravel()
    return np.unique(rows)


def correctedRows(badPixels: np.ndarray, irpN: int, irpOffset: int,
                  readOrders: Sequence[int], nchan: int = 32) -> np.ndarray:
    """IRP1 rows actually repaired by the per-channel bad-pixel correction.

    Mirrors ``PfsIsrTask.constructFullIrp``: the reference phase within each
    group is ``irpOffset - 1`` for forward-read channels and ``irpN - irpOffset``
    for reverse-read ones, with the direction taken per channel from
    ``readOrders`` (``getH4channelReadOrder()``). This is the set the production
    correction repairs, unlike `correctedRowsForOffset`, which assumes a single
    global offset.

    Parameters
    ----------
    badPixels : `numpy.ndarray`
        Bad pixels in IRP1 coordinates.
    irpN : `int`
        Data:IRP ratio.
    irpOffset : `int`
        The header ``W_H4IRPO`` value.
    readOrders : sequence of `int`
        ``(even, odd)`` flip flags; channel ``c`` is reversed if
        ``readOrders[c % 2]``.
    nchan : `int`
        Number of readout channels.
    """
    badPixels = np.asarray(badPixels, dtype=int)
    if irpN <= 1:
        return np.unique(badPixels)
    forward = irpOffset - 1
    reversed_ = irpN - irpOffset
    chanRows = H4_SIZE // nchan
    blocks = []
    for c in range(nchan):
        low = c * chanRows
        offset = reversed_ if readOrders[c % 2] else forward
        inchan = badPixels[(badPixels >= low) & (badPixels < low + chanRows)] - low
        relevant = inchan[(inchan - offset) % irpN == 0]
        samples = (relevant - offset) // irpN
        blocks.append(low + (samples[:, None] * irpN + np.arange(irpN)).ravel())
    return np.unique(np.concatenate(blocks)) if blocks else np.zeros(0, dtype=int)


def collapseToIrp1Pixels(rows, irpN: int, irpOffset: int,
                         readOrders: Sequence[int], nchan: int = 32) -> np.ndarray:
    """Collapse IRP4-expanded noisy-row blocks to single IRP1 bad-pixel ordinates.

    A bad reference pixel in IRP4 data appears as a block of ``irpN`` adjacent
    rows (the repeat-expanded sample). This returns, for each distinct block, the
    single IRP1 coordinate of the physically sampled reference pixel -- the value
    to add to ``badRefPixels.yaml`` -- using the per-channel reference phase
    (forward ``irpOffset-1``, reversed ``irpN-irpOffset``).

    Parameters
    ----------
    rows : array-like
        Noisy IRP1 rows (e.g. ``validateCorrection(...).remainingRows``).
    irpN, irpOffset : `int`
        Data:IRP ratio and the header ``W_H4IRPO``.
    readOrders : sequence of `int`
        ``(even, odd)`` flip flags from ``getH4channelReadOrder()``.
    nchan : `int`
        Number of readout channels.

    Returns
    -------
    pixels : `numpy.ndarray`
        Sorted, de-duplicated IRP1 bad-pixel ordinates.
    """
    rows = np.unique(np.asarray(rows, dtype=int))
    if irpN <= 1:
        return rows
    forward = irpOffset - 1
    reversed_ = irpN - irpOffset
    chanRows = H4_SIZE // nchan
    out = set()
    for r in rows.tolist():
        chanLow = (r // chanRows) * chanRows
        offset = reversed_ if readOrders[(r // chanRows) % 2] else forward
        group = (r - chanLow) // irpN
        out.add(chanLow + group * irpN + offset)
    return np.array(sorted(out), dtype=int)


def offsetCoverage(noisyRows, badPixels, irpN: int) -> Dict[int, dict]:
    """Score each candidate offset by how many noisy rows its correction repairs.

    For every ``offset`` in ``0..irpN-1`` (just 0 for IRP1), find the rows the
    known bad-pixel list would repair and intersect with the observed noisy
    rows.

    Returns
    -------
    coverage : `dict`
        ``offset -> {"nCovered", "coveredFrac", "uncovered"}``, where
        ``uncovered`` lists the noisy rows that offset fails to repair.
    """
    noisySet = set(np.asarray(noisyRows, dtype=int).tolist())
    coverage = {}
    for offset in range(max(irpN, 1)):
        fixed = set(correctedRowsForOffset(badPixels, irpN, offset).tolist())
        covered = noisySet & fixed
        coverage[offset] = dict(
            nCovered=len(covered),
            coveredFrac=(len(covered) / len(noisySet)) if noisySet else 0.0,
            uncovered=sorted(noisySet - fixed),
        )
    return coverage


def offsetCoverageByParity(noisyRows, badPixels, irpN: int,
                           nchan: int = 32) -> Dict[str, Dict[int, dict]]:
    """`offsetCoverage` computed separately for even- and odd-numbered channels.

    The H4 temporal readout direction alternates between even and odd channels.
    If the IRP reference pixels inherit that flip, the within-group phase is
    ``irpOffset`` in one parity and ``irpN-1-irpOffset`` in the other, so even
    and odd channels prefer *different* offsets and no single global offset fits.
    Compare the best offset of the two returned coverages to test for that: if
    they differ (and sum to ``irpN-1``), the readout flip is the culprit.

    Parameters
    ----------
    noisyRows, badPixels : array-like
        IRP1 row indices (0..H4_SIZE-1).
    irpN : `int`
        Data:IRP ratio.
    nchan : `int`
        Number of readout channels.

    Returns
    -------
    coverage : `dict`
        ``{"even": coverage, "odd": coverage}`` (see `offsetCoverage`).
    """
    chanHeight = H4_SIZE // nchan
    noisyRows = np.asarray(noisyRows, dtype=int)
    badPixels = np.asarray(badPixels, dtype=int)
    out = {}
    for name, parity in (("even", 0), ("odd", 1)):
        nr = noisyRows[(noisyRows // chanHeight) % 2 == parity]
        bp = badPixels[(badPixels // chanHeight) % 2 == parity]
        out[name] = offsetCoverage(nr, bp, irpN)
    return out


@dataclass
class IrpAlignment:
    """Result of `alignByCorrection`."""

    detectorName: str
    irpN: int
    headerOffset: int
    threshold: float
    noisyRows: np.ndarray
    badPixels: np.ndarray
    coverage: Dict[int, dict]
    byParity: Optional[Dict[str, Dict[int, dict]]] = None
    readOrders: Optional[Sequence[int]] = None
    nchan: int = 32
    dataId: Optional[dict] = None

    @property
    def bestOffset(self) -> int:
        """Global offset whose correction covers the most observed noisy rows."""
        return max(self.coverage, key=lambda o: self.coverage[o]["nCovered"])

    def correctedRows(self) -> np.ndarray:
        """IRP1 rows the actual per-channel correction repairs (see `correctedRows`)."""
        return correctedRows(self.badPixels, self.irpN, self.headerOffset,
                             self.readOrders, self.nchan)

    @staticmethod
    def _best(cov: Dict[int, dict]) -> int:
        return max(cov, key=lambda o: cov[o]["nCovered"])

    def report(self) -> str:
        n = len(self.noisyRows)
        lines = [
            f"detector {self.detectorName}: irpN={self.irpN}, "
            f"header W_H4IRPO={self.headerOffset}, threshold={self.threshold}",
            f"flagged {n} noisy IRP rows; known IRP1 bad pixels: {len(self.badPixels)}",
        ]

        if self.byParity is not None:
            # Forward/reversed channels alternate, so the alignment is per parity.
            evenBest = self._best(self.byParity["even"])
            oddBest = self._best(self.byParity["odd"])
            nEven = self.byParity["even"][evenBest]["nCovered"]
            nOdd = self.byParity["odd"][oddBest]["nCovered"]
            lines.append(f"per-parity best offsets: even={evenBest} ({nEven} rows), "
                         f"odd={oddBest} ({nOdd} rows)")
            fwd, rev = self.headerOffset - 1, self.irpN - self.headerOffset
            if evenBest != oddBest:
                lines.append(
                    f"even/odd offsets differ -> temporal readout flip; expected from "
                    f"header: forward={fwd}, reversed={rev} "
                    f"({'matches' if {evenBest, oddBest} == {fwd, rev} else 'does NOT match'})."
                )
            covered = nEven + nOdd
            uncovered = sorted(set(self.byParity['even'][evenBest]['uncovered']) |
                               set(self.byParity['odd'][oddBest]['uncovered']))
            if uncovered:
                lines.append(
                    f"{n - covered} noisy rows unrepaired by the per-parity correction "
                    f"(e.g. {uncovered[:12]}): candidate new bad pixels for badRefPixels.yaml."
                )
            else:
                lines.append("every noisy row is repaired by the per-parity correction.")
        else:
            best = self.bestOffset
            for offset in sorted(self.coverage):
                c = self.coverage[offset]
                mark = "  <-- best" if offset == best else ""
                lines.append(f"  offset {offset}: corrects {c['nCovered']}/{n} noisy rows "
                             f"({c['coveredFrac']:.0%}){mark}")
        return "\n".join(lines)


def alignByCorrection(task, raw, threshold: float = 2.0,
                      r0: int = 0, r1: int = 1) -> IrpAlignment:
    """Pick the IRP offset whose correction best repairs the observed noisy rows.

    Finds the strikingly noisy IRP rows (`noisyIrpRows`), then for each candidate
    ``offset`` computes the rows the known bad-pixel list would repair
    (`correctedRowsForOffset`) and scores by how many noisy rows that covers. The
    offset with the most coverage is the alignment; rows left uncovered at every
    offset are candidate new bad pixels.

    Parameters
    ----------
    task : `lsst.obs.pfs.isrTask.PfsIsrTask`
        Task used to build the IRP cube.
    raw : `lsst.obs.pfs.PfsRaw`
        The raw IRP4 dark ramp.
    threshold : `float`
        Noise threshold (ADU) passed to `noisyIrpRows`.
    r0, r1 : `int`
        Reads to difference for the CDS.

    Returns
    -------
    result : `IrpAlignment`

    Raises
    ------
    ValueError
        If ``raw`` is IRP1 data: every reference pixel is read, so there is no
        sampling phase to align. Use `noisyIrpRows` directly to find bad pixels.
    """
    irpN = raw.irpN
    if irpN <= 1:
        raise ValueError(
            f"alignByCorrection needs IRP4+ data, but this ramp is IRP1 "
            f"(irpN={irpN}): every reference pixel is read, so there is no "
            f"offset to determine. Use noisyIrpRows() to find bad pixels."
        )
    detectorName = raw.detector.getName()
    noisyRows, _ = noisyIrpRows(task, raw, threshold=threshold, r0=r0, r1=r1)
    badPixels = loadBadRefPixels(detectorName)
    coverage = offsetCoverage(noisyRows, badPixels, irpN)
    byParity = offsetCoverageByParity(noisyRows, badPixels, irpN, nchan=raw.nchan)
    return IrpAlignment(detectorName, irpN, raw.irpOffset, threshold,
                        noisyRows, badPixels, coverage, byParity=byParity,
                        readOrders=raw.getH4channelReadOrder(), nchan=raw.nchan)


@dataclass
class CorrectionValidation:
    """Result of `validateCorrection`."""

    detectorName: str
    threshold: float
    noisyRows: np.ndarray       # flagged with the correction off
    repairedRows: np.ndarray    # dropped below threshold once corrected
    remainingRows: np.ndarray   # still above threshold -> candidate new bad pixels
    excessOff: np.ndarray       # per-row excess noise, correction off
    excessOn: np.ndarray        # per-row excess noise, correction on
    dataId: Optional[dict] = None

    def report(self) -> str:
        n = len(self.noisyRows)
        lines = [
            f"detector {self.detectorName}: threshold={self.threshold}",
            f"noisy rows (correction off): {n}",
            f"repaired by correction:      {len(self.repairedRows)}/{n}",
        ]
        if len(self.remainingRows):
            lines.append(
                f"still noisy ({len(self.remainingRows)}): {self.remainingRows.tolist()}"
                f"  -> candidate new bad pixels for badRefPixels.yaml"
            )
        else:
            lines.append("every noisy row was repaired.")
        return "\n".join(lines)


def validateCorrection(task, raw, threshold: float = 2.0,
                       r0: int = 0, r1: int = 1) -> CorrectionValidation:
    """Check that the IRP bad-pixel correction actually repairs the noisy rows.

    Builds the IRP1 CDS difference with the bad-pixel correction disabled and
    enabled, and compares the per-row noise above each channel's baseline (see
    `excessRowNoise`). Rows that fall from above to below ``threshold`` once
    corrected are repaired; rows that stay above are candidate new bad pixels
    not yet in ``badRefPixels.yaml``.

    Parameters
    ----------
    task : `lsst.obs.pfs.isrTask.PfsIsrTask`
        Task used to build the IRP cube (e.g. from `makeQuickIsrTask`).
    raw : `lsst.obs.pfs.PfsRaw`
        The raw ramp.
    threshold : `float`
        Excess (above channel-median iqrStd, ADU) for both flagging and the
        still-noisy test.
    r0, r1 : `int`
        Reads to difference for the CDS (default 0 and 1).

    Returns
    -------
    result : `CorrectionValidation`
    """
    from lsst.obs.pfs.nirBadRefPixels import NirBadRefPixels

    det = raw.detector.getName()
    cfg = task.config.h4
    prevDo = cfg.doIRPbadPixels
    prevCalib = getattr(task, "_badRefPixels", None)
    try:
        cfg.doIRPbadPixels = False
        off = task.makeRawIrpNcube(raw, r0=r0, r1=r1)
        excessOff = excessRowNoise(iqrStd(off[-1] - off[0], axis=1), raw.nchan)
        # The corrected pass reaches getBadIRPpixels(), which reads the calib the
        # butler would supply; provide it from the current file-based list.
        cfg.doIRPbadPixels = True
        task._badRefPixels = NirBadRefPixels.fromList(loadBadRefPixels(det), det)
        on = task.makeRawIrpNcube(raw, r0=r0, r1=r1)
        excessOn = excessRowNoise(iqrStd(on[-1] - on[0], axis=1), raw.nchan)
    finally:
        cfg.doIRPbadPixels = prevDo
        task._badRefPixels = prevCalib

    noisy = np.where(excessOff > threshold)[0]
    repaired = noisy[excessOn[noisy] <= threshold]
    remaining = noisy[excessOn[noisy] > threshold]
    return CorrectionValidation(raw.detector.getName(), threshold,
                                noisy, repaired, remaining, excessOff, excessOn)


def rampExcess(task, raw, nreads: Optional[int] = None,
               aggregate=np.mean, nchan: Optional[int] = None) -> np.ndarray:
    """Per-row IRP noise across a ramp's reads, aggregated over read-pair diffs.

    Streams the reference reads, forms consecutive-read CDS diffs, computes each
    diff's per-row excess over the channel baseline (`excessRowNoise`), and
    aggregates over the diffs. The default ``mean`` keeps a persistently-noisy
    reference pixel high while diluting single-read transients, and avoids the
    quantization of ``median`` (``iqrStd`` on integer-ADU data lands on ~0.741-ADU
    steps, so a median of reads stays quantized). ``max`` is the most sensitive
    but transient-prone. Intended for IRP1 darks, where every reference pixel is
    read and each flagged row is directly a bad-pixel ordinate.

    Parameters
    ----------
    task : `lsst.obs.pfs.isrTask.PfsIsrTask`
        Task whose ``makeRawIrpArray`` reads one IRP1 reference image.
    raw : `lsst.obs.pfs.PfsRaw`
        The raw dark ramp.
    nreads : `int`, optional
        Use this many evenly-spaced reads (minimum 2); ``None`` uses all reads.
    aggregate : callable
        Reduction over the per-diff excess stack along axis 0 (e.g. ``np.mean``,
        ``np.max``, ``np.median``), called as ``aggregate(stack, axis=0)``.
    nchan : `int`, optional
        Number of channels; defaults to ``raw.nchan``.

    Returns
    -------
    excess : `numpy.ndarray`
        Aggregated per-row excess noise, length ``H4_SIZE``.
    """
    return aggregate(_rampExcessStack(task, raw, nreads, nchan), axis=0)


def _rampExcessStack(task, raw, nreads=None, nchan=None, doCorrect=False) -> np.ndarray:
    """Per-diff, per-row excess-noise stack for a ramp, shape ``(ndiff, H4_SIZE)``.

    Streams the reference reads, forms consecutive-read CDS diffs, and returns
    each diff's `excessRowNoise`. The bad-pixel correction is disabled by default
    (``doCorrect=False``) so known-bad pixels stay visible; set ``doCorrect=True``
    to measure the residual after the current ``badRefPixels.yaml`` is applied.
    """
    nchan = raw.nchan if nchan is None else nchan
    numReads = raw.getNumReads()
    if nreads is None:
        reads = list(range(numReads))
    else:
        if nreads < 2:
            raise ValueError(f"nreads must be >= 2; got {nreads}")
        reads = sorted(set(np.linspace(0, numReads - 1, nreads).round().astype(int).tolist()))
    if len(reads) < 2:
        raise ValueError(f"need at least two reads; got {reads}")

    cfg = task.config.h4
    prev = cfg.doIRPbadPixels
    cfg.doIRPbadPixels = doCorrect
    try:
        prevImg = task.makeRawIrpArray(raw, reads[0])
        excesses = []
        for rn in reads[1:]:
            img = task.makeRawIrpArray(raw, rn)
            excesses.append(excessRowNoise(iqrStd(img - prevImg, axis=1), nchan))
            prevImg = img
    finally:
        cfg.doIRPbadPixels = prev
    return np.stack(excesses)


@dataclass
class BadRefPixelScan:
    """Proposed bad reference pixels for one detector (see `scanBadRefPixels`)."""

    detectorName: str
    threshold: float
    maxThreshold: float
    minVisits: int
    nVisits: int
    counts: Dict[int, int]      # row -> number of darks that flagged it
    pixels: np.ndarray          # rows flagged in >= minVisits darks (sorted)
    meanMetric: np.ndarray      # per-row mean excess, max over darks
    maxMetric: np.ndarray       # per-row max excess, max over darks

    def report(self) -> str:
        known = set(loadBadRefPixels(self.detectorName).tolist())
        proposed = self.pixels.tolist()
        new = [p for p in proposed if p not in known]
        missing = sorted(known - set(proposed))
        return "\n".join([
            f"detector {self.detectorName}: {len(proposed)} proposed bad ref pixels "
            f"from {self.nVisits} dark(s) (mean>{self.threshold} or max>{self.maxThreshold}, "
            f"minVisits={self.minVisits})",
            f"  {len(new)} not in current badRefPixels.yaml: {new[:20]}",
            f"  {len(missing)} current entries not re-flagged: {missing[:20]}",
        ])


def scanBadRefPixels(task, raws, threshold: float = 2.0, maxThreshold: float = 4.0,
                     nreads: Optional[int] = None, minVisits: int = 1) -> BadRefPixelScan:
    """Propose bad reference pixels for one detector from IRP1 dark ramps.

    For each dark, the per-read excess-noise stack gives two per-row metrics: the
    ``mean`` over reads (catches persistently-noisy pixels, de-quantized) and the
    ``max`` over reads (catches intermittent / telegraph pixels that the mean
    dilutes). A row is flagged when its mean exceeds ``threshold`` OR its max
    exceeds ``maxThreshold``. A pixel is proposed when flagged in at least
    ``minVisits`` of the supplied darks (raise it above 1 to suppress transient
    single-dark hits).

    Parameters
    ----------
    task : `lsst.obs.pfs.isrTask.PfsIsrTask`
        Task used to read the reference images (e.g. from `makeQuickIsrTask`).
    raws : sequence of `lsst.obs.pfs.PfsRaw`
        IRP1 dark ramps for the *same* detector.
    threshold : `float`
        Mean-excess threshold (ADU) for persistently-noisy pixels.
    maxThreshold : `float`
        Max-excess threshold (ADU) for intermittent pixels.
    nreads : `int`, optional
        Reads per ramp (see `rampExcess`).
    minVisits : `int`
        Minimum number of darks that must flag a pixel for it to be proposed.

    Returns
    -------
    scan : `BadRefPixelScan`
    """
    raws = list(raws)
    if not raws:
        raise ValueError("need at least one dark")
    detName = raws[0].detector.getName()

    counts: Dict[int, int] = {}
    meanMetric = maxMetric = None
    for raw in raws:
        if raw.irpN != 1:
            raise ValueError(
                f"scanBadRefPixels expects IRP1 darks; {detName} ramp is irpN={raw.irpN}. "
                f"For IRP4 data flag with noisyIrpRows then collapseToIrp1Pixels."
            )
        if raw.detector.getName() != detName:
            raise ValueError(f"mixed detectors: {raw.detector.getName()} != {detName}")
        stack = _rampExcessStack(task, raw, nreads)
        meanExc, maxExc = stack.mean(axis=0), stack.max(axis=0)
        meanMetric = meanExc if meanMetric is None else np.maximum(meanMetric, meanExc)
        maxMetric = maxExc if maxMetric is None else np.maximum(maxMetric, maxExc)
        flagged = (meanExc > threshold) | (maxExc > maxThreshold)
        for r in np.where(flagged)[0].tolist():
            counts[r] = counts.get(r, 0) + 1

    pixels = np.array(sorted(r for r, c in counts.items() if c >= minVisits), dtype=int)
    return BadRefPixelScan(detName, threshold, maxThreshold, minVisits, len(raws),
                           counts, pixels, meanMetric, maxMetric)


def scanBadRefPixelsAcrossVisits(butler, visits, arm: str = "n",
                                 spectrographs=(1, 2, 3, 4), threshold: float = 2.0,
                                 maxThreshold: float = 4.0, nreads: Optional[int] = None,
                                 minVisits: int = 1, task=None):
    """Scan many IRP1 dark visits and return a `BadRefPixelScan` per detector.

    Convenience wrapper over `scanBadRefPixels`: for each spectrograph it loads
    the ``raw`` of every visit it can find and scans them together. Visits with
    no dataset for a given detector are skipped.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler with collections set.
    visits : sequence of `int`
        Dark visit numbers.
    arm : `str`
        Arm (``"n"`` for the H4s).
    spectrographs : sequence of `int`
        Spectrograph numbers to scan (one detector each).
    threshold, maxThreshold, nreads, minVisits :
        Passed through to `scanBadRefPixels`.
    task : `lsst.obs.pfs.isrTask.PfsIsrTask`, optional
        Task to use; a `makeQuickIsrTask` is built if omitted.

    Returns
    -------
    scans : `list` of `BadRefPixelScan`
        One per spectrograph that had at least one dark.
    """
    task = task or makeQuickIsrTask()
    scans = []
    for spec in spectrographs:
        raws, missing = [], []
        for v in visits:
            try:
                raws.append(butler.get("raw", dataId=dict(visit=v, arm=arm, spectrograph=spec)))
            except Exception:
                missing.append(v)
        if missing:
            print(f"spectrograph {spec}: skipped {len(missing)} visits with no raw")
        if raws:
            scans.append(scanBadRefPixels(task, raws, threshold=threshold,
                                          maxThreshold=maxThreshold, nreads=nreads,
                                          minVisits=minVisits))
    return scans


def proposedBadRefPixelsYaml(scans, metadata=None) -> str:
    """Render `BadRefPixelScan` results as ``badRefPixels.yaml`` text.

    If ``metadata`` (a provenance dict: who/when/visits/dates/thresholds) is
    given, it is emitted under a top-level ``metadata`` key. The pixel-list
    loaders (`loadBadRefPixels`, ``PfsIsrTask.loadBadIRPpixels``) look up by
    detector name, so that key is inert for them.
    """
    data = {}
    if metadata is not None:
        data["metadata"] = metadata
    for s in scans:
        data[s.detectorName] = [int(p) for p in s.pixels]
    return yaml.safe_dump(data, default_flow_style=False, sort_keys=False)


def _plotExcessProfile(ax, metric, scan, rows, ylabel, hline):
    """Draw one excess-vs-row panel from a combined metric + kept/pruned/new marks."""
    det = scan.detectorName
    sl = rows if rows is not None else slice(None)
    lo = sl.start or 0
    hi = sl.stop if sl.stop is not None else len(metric)
    x = np.arange(len(metric))
    ax.plot(x[sl], metric[sl], color="0.5", lw=0.7, label="excess (max over darks)")
    ax.axhline(hline, color="C3", ls="--", lw=0.8, label=f"threshold={hline}")

    old = set(loadBadRefPixels(det).tolist())
    proposed = set(scan.pixels.tolist())

    def mark(idx, style, label, **kw):
        arr = np.array(sorted(idx), dtype=int)
        inr = arr[(arr >= lo) & (arr < hi)]
        ax.plot(inr, metric[inr], style, label=f"{label} ({len(arr)})", **kw)

    mark(old & proposed, "rx", "kept", ms=6)
    mark(old - proposed, "o", "pruned", ms=7, mfc="none", mec="C1", mew=1.3)
    mark(proposed - old, "g+", "new", ms=8, mew=1.6)
    ax.set_ylabel(f"iqrStd excess\n({ylabel})")
    ax.legend(loc="upper right", fontsize=8)


def plotScan(scan, rows: Optional[slice] = None):
    """Two-row plot of a scan's combined per-row excess: mean and worst read.

    Uses the scan's combined-over-darks metrics, so every flagged pixel shows its
    spike (unlike a single-dark trace, where a pixel that fired in a *different*
    dark would look flat). Current-list pixels are marked with red ×, newly
    proposed ones with green +.

    - Top (mean over reads): persistently-noisy pixels.
    - Bottom (max over reads = the worst read): pixels from *occasional* offsets,
      flat in the mean but spiking in a single bad read.

    Parameters
    ----------
    scan : `BadRefPixelScan`
        Result from `scanBadRefPixels`.
    rows : `slice`, optional
        Restrict to a row range (e.g. ``slice(1650, 1850)`` to zoom).

    Returns
    -------
    fig : `matplotlib.figure.Figure`
    """
    import matplotlib.pyplot as plt

    fig, (axMean, axMax) = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    _plotExcessProfile(axMean, scan.meanMetric, scan, rows, "mean over reads", scan.threshold)
    _plotExcessProfile(axMax, scan.maxMetric, scan, rows, "max over reads (worst read)",
                       scan.maxThreshold)
    axMax.set_xlabel("IRP row")
    fig.suptitle(f"{scan.detectorName} IRP1 dark bad-ref-pixel scan "
                 f"({scan.nVisits} darks, minVisits={scan.minVisits})")
    fig.tight_layout()
    return fig


def plotIrpRowNoise(metric: np.ndarray, threshold: float,
                    alignment: Optional[IrpAlignment] = None,
                    rows: Optional[slice] = None, ax=None,
                    showChannels: bool = False, nchan: int = 32,
                    readOrders: Optional[Sequence[int]] = None):
    """Plot the per-row IRP noise metric with the threshold and flagged rows.

    Parameters
    ----------
    metric : `numpy.ndarray`
        Per-row metric from `noisyIrpRows` (excess iqrStd above channel baseline).
    threshold : `float`
        The flagging threshold; drawn as a horizontal line.
    alignment : `IrpAlignment`, optional
        If given, shade the rows repaired at its best offset.
    rows : `slice`, optional
        Restrict the plot to a row range (e.g. ``slice(3800, 3950)`` to zoom on
        the 3896-3899 block).
    ax : `matplotlib.axes.Axes`, optional
        Axes to draw on; created if omitted.
    showChannels : `bool`
        Overlay shaded bands at the ``nchan`` channel boundaries. With
        ``readOrders`` the bands are coloured by temporal readout direction
        (forward vs reversed); otherwise alternate channels are shaded grey.
    nchan : `int`
        Number of readout channels (channel height = ``len(metric) // nchan``).
    readOrders : sequence of `int`, optional
        The 2-tuple from ``pfsRaw.getH4channelReadOrder()`` ``(even, odd)``;
        channel ``c`` is reversed if ``readOrders[c % 2]``. Used to colour the
        channel bands so the alternating readout direction is visible.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(12, 4))
    sl = rows if rows is not None else slice(None)
    lo = sl.start or 0
    hi = sl.stop if sl.stop is not None else len(metric)
    x = np.arange(len(metric))

    if showChannels:
        chanHeight = len(metric) // nchan
        seen = set()
        for c in range(nchan):
            clo, chi = c * chanHeight, (c + 1) * chanHeight
            if chi <= lo or clo >= hi:
                continue
            if readOrders is not None:
                reversed_ = bool(readOrders[c % 2])
                color = "C1" if reversed_ else "C2"
                label = "reversed channel" if reversed_ else "forward channel"
            else:
                if c % 2 == 0:
                    continue
                color, label = "0.5", "odd channel"
            ax.axvspan(clo, chi, color=color, alpha=0.12, zorder=0,
                       label=label if label not in seen else None)
            seen.add(label)

    ax.plot(x[sl], metric[sl], color="0.4", lw=0.8, label="iqrStd excess")
    ax.axhline(threshold, color="C3", ls="--", lw=1, label=f"threshold={threshold}")

    if alignment is not None:
        # Rows the actual per-channel correction repairs (forward/reversed phase
        # per channel), NOT a single global offset.
        fixed = alignment.correctedRows()
        first = True
        for r in fixed[(fixed >= lo) & (fixed < hi)]:
            ax.axvline(r, color="C0", alpha=0.4, lw=0.8,
                       label="repaired by correction" if first else None)
            first = False
        fwd, rev = alignment.headerOffset - 1, alignment.irpN - alignment.headerOffset
        ax.set_title(f"{alignment.detectorName} IRP{alignment.irpN} row noise "
                     f"(phase fwd={fwd}/rev={rev}, header W_H4IRPO={alignment.headerOffset})")
    ax.set_xlabel("IRP row")
    ax.set_ylabel("iqrStd excess over channel median (ADU)")
    ax.legend()
    return ax


def runBadRefPixelSurvey(butler, visits, outdir: str = ".", arm: str = "n",
                         spectrographs=(1, 2, 3, 4), threshold: float = 2.0,
                         maxThreshold: float = 4.0, minVisits: int = 3,
                         nreads: Optional[int] = 30, task=None):
    """Scan all detectors over dark visits and write plots + proposed YAML to ``outdir``.

    For each spectrograph it runs `scanBadRefPixels` over the visits, writes a
    ``scan_<det>.png`` (`plotScan`) and, once all detectors are done,
    ``badRefPixels_proposed.yaml`` (`proposedBadRefPixelsYaml`). Returns the list
    of `BadRefPixelScan`. ``outdir`` is created if needed; give each run its own
    directory. Defaults to ``nreads=30`` (evenly-spaced, to keep runtime down;
    use ``None`` for the whole ramp) and the ``mean>2 or max>4``, ``minVisits=3``
    operating point.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler with collections set.
    visits : sequence of `int`
        Dark visit numbers.
    outdir : `str`
        Directory to write plots and YAML into.
    (other parameters as in `scanBadRefPixels`.)

    Returns
    -------
    scans : `list` of `BadRefPixelScan`
    """
    import os
    import datetime
    import getpass
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    task = task or makeQuickIsrTask()
    os.makedirs(outdir, exist_ok=True)
    visits = list(visits)
    visitDates = {}
    scans = []
    for spec in spectrographs:
        raws = []
        for v in visits:
            try:
                raw = butler.get("raw", dataId=dict(visit=v, arm=arm, spectrograph=spec))
            except Exception:
                continue
            raws.append(raw)
            visitDates.setdefault(v, _obsDate(raw))
        if not raws:
            continue
        scan = scanBadRefPixels(task, raws, threshold=threshold, maxThreshold=maxThreshold,
                                nreads=nreads, minVisits=minVisits)
        scans.append(scan)
        fig = plotScan(scan)
        fig.savefig(os.path.join(outdir, f"scan_{scan.detectorName}.png"),
                    dpi=110, bbox_inches="tight")
        plt.close(fig)

    metadata = {
        "generatedBy": getpass.getuser(),
        "generatedAt": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "arm": arm,
        "visits": visits,
        "visitDates": {v: visitDates.get(v) for v in visits},
        "threshold": threshold,
        "maxThreshold": maxThreshold,
        "minVisits": minVisits,
        "nreads": nreads,
    }
    with open(os.path.join(outdir, "badRefPixels_proposed.yaml"), "w") as f:
        f.write(proposedBadRefPixelsYaml(scans, metadata=metadata))
    return scans


def _obsDate(raw):
    """Observation date (YYYY-MM-DD) of a raw, or None if unavailable."""
    try:
        return str(raw.obsInfo.datetime_begin.iso)[:10]
    except Exception:
        return None


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    from lsst.daf.butler import Butler

    parser = argparse.ArgumentParser(description="Determine the IRP4 bad-pixel offset")
    parser.add_argument("repo", help="Butler repository")
    parser.add_argument("--collections", required=True, nargs="+", help="Input collections")
    parser.add_argument("--visit", type=int, required=True, help="Visit of an IRP4 dark")
    parser.add_argument("--arm", default="n")
    parser.add_argument("--spectrograph", type=int, required=True)
    parser.add_argument("--threshold", type=float, default=2.0,
                        help="excess iqrStd over channel median (ADU) to flag a row")
    args = parser.parse_args(argv)

    butler = Butler(args.repo, collections=args.collections)
    raw = butler.get("raw", dataId=dict(visit=args.visit, arm=args.arm,
                                        spectrograph=args.spectrograph))
    task = makeQuickIsrTask()
    result = alignByCorrection(task, raw, threshold=args.threshold)
    result.dataId = dict(visit=args.visit, arm=args.arm, spectrograph=args.spectrograph)
    print(result.report())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
