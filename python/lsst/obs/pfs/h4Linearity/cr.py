"""Post-linearization cosmic-ray detection and repair on H4 ramps.

Detect single-read positive rate outliers in a linearized cumulative ramp,
repair them in place by interpolating the affected delta to the per-pixel
median rate, and report a per-pixel CR mask. Replaces the older
delta-space CR path on H4 ramps.

The pure-numpy algorithm is here; the wrapper that integrates with
``PfsIsrTask`` lives in ``lsst.obs.pfs.isrTask``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


DEFAULT_N_SIGMA = 4.0
DEFAULT_SIGMA_FLOOR_ADU = 8.0
DEFAULT_EXCESS_FLOOR_ADU = 15.0
DEFAULT_MIN_READS = 8


@dataclass
class CRRepairResult:
    """Outcome of a `detectAndRepair` call.

    Attributes
    ----------
    nFlagged : int
    flagMask : np.ndarray
        Boolean ``(H, W)``; True where a CR was flagged.
    crReadIdx : np.ndarray
        Int ``(H, W)``; per-pixel index ``k`` into ``deltas`` of the
        per-pixel max-delta read. Always populated.
    excess : np.ndarray
        Float ``(H, W)``; ``max_delta - median_delta`` (i.e., the amount
        subtracted from cumulative reads ``>= k+1``). Populated for every
        candidate (pixel that passed the pre-screen); zero elsewhere.
    candidateMask : np.ndarray, optional
        Diagnostics-only. Boolean ``(H, W)``; True for pixels that passed
        the cheap pre-screen and had per-pixel statistics computed.
        ``None`` unless ``returnDiagnostics=True``.
    medianDelta : np.ndarray, optional
        Diagnostics-only. Float ``(H, W)``; per-pixel median delta.
        Populated only at candidate pixels; zero elsewhere. ``None``
        unless ``returnDiagnostics=True``.
    sigma : np.ndarray, optional
        Diagnostics-only. Float ``(H, W)``; per-pixel
        ``max(1.4826 * MAD, sigmaFloorADU)``. Populated only at candidate
        pixels; zero elsewhere. ``None`` unless ``returnDiagnostics=True``.
    """

    nFlagged: int
    flagMask: np.ndarray
    crReadIdx: np.ndarray
    excess: np.ndarray
    candidateMask: Optional[np.ndarray] = None
    medianDelta: Optional[np.ndarray] = None
    sigma: Optional[np.ndarray] = None


def detectAndRepair(
    flux: np.ndarray,
    *,
    goodPixelMask: np.ndarray,
    nSigma: float = DEFAULT_N_SIGMA,
    sigmaFloorADU: float = DEFAULT_SIGMA_FLOOR_ADU,
    excessFloorADU: float = DEFAULT_EXCESS_FLOOR_ADU,
    minReads: int = DEFAULT_MIN_READS,
    returnDiagnostics: bool = False,
) -> CRRepairResult:
    """Detect and repair single-read positive rate outliers in-place.

    For each pixel where ``goodPixelMask`` is True, compute per-read deltas
    of the linearized cumulative cube and find the read ``k`` with the
    largest positive delta. Flag the pixel as a CR when both:

        (max_delta - median_delta) > nSigma * sigma
        (max_delta - median_delta) > excessFloorADU

    where ``sigma = max(1.4826 * MAD(deltas), sigmaFloorADU)``. The MAD
    floor keeps the threshold from collapsing on faint pixels whose
    empirical scatter is read-noise-limited.

    Pre-screen: pixels whose ``max_delta`` is at or below ``excessFloorADU``
    cannot satisfy the (excess > floor) criterion (since
    ``excess <= max_delta`` whenever ``median >= 0``), so the per-pixel
    median/MAD computation is skipped for them. On a typical dark frame
    this filters out ~99% of pixels, giving a large speedup on the
    median-along-axis-0 step that otherwise dominates.

    Repair: subtract ``(max_delta - median_delta)`` from cumulative reads
    ``>= k + 1``, equivalent to replacing the outlier delta with the
    per-pixel median.

    Parameters
    ----------
    flux : np.ndarray
        Linearized cumulative ramp, shape ``(N, H, W)``. Modified in place.
    goodPixelMask : np.ndarray
        Boolean ``(H, W)`` mask. Only pixels where this is True are
        candidates for CR flagging.
    nSigma, sigmaFloorADU, excessFloorADU, minReads
        See module-level defaults.
    returnDiagnostics : bool
        If True, populate ``candidateMask``, ``medianDelta``, and
        ``sigma`` on the returned result so the threshold behavior can be
        inspected (histograms etc.). Off by default to keep the
        production path light.

    Returns
    -------
    CRRepairResult
    """
    if flux.ndim != 3:
        raise ValueError(f"flux must be 3-D (N, H, W); got {flux.shape}")
    if goodPixelMask.shape != flux.shape[1:]:
        raise ValueError(
            f"goodPixelMask shape {goodPixelMask.shape} != flux H,W {flux.shape[1:]}"
        )

    nReads, H, W = flux.shape

    def _empty(candidateMask=None, medianDelta=None, sigma=None):
        return CRRepairResult(
            nFlagged=0,
            flagMask=np.zeros((H, W), dtype=bool),
            crReadIdx=np.zeros((H, W), dtype=np.int32),
            excess=np.zeros((H, W), dtype=np.float32),
            candidateMask=candidateMask,
            medianDelta=medianDelta,
            sigma=sigma,
        )

    if nReads < minReads:
        return _empty()

    deltas = np.diff(flux, axis=0)
    kCR = np.argmax(deltas, axis=0).astype(np.int32)
    maxDelta = np.take_along_axis(deltas, kCR[None], axis=0)[0]

    candidateMask = (maxDelta > excessFloorADU) & goodPixelMask

    medianFull = np.zeros((H, W), dtype=np.float32) if returnDiagnostics else None
    sigmaFull = np.zeros((H, W), dtype=np.float32) if returnDiagnostics else None
    excessFull = np.zeros((H, W), dtype=np.float32)
    flag = np.zeros((H, W), dtype=bool)

    nCand = int(candidateMask.sum())
    if nCand == 0:
        return _empty(
            candidateMask=candidateMask if returnDiagnostics else None,
            medianDelta=medianFull,
            sigma=sigmaFull,
        )

    ys, xs = np.where(candidateMask)
    candDeltas = deltas[:, ys, xs]
    candMedian = np.median(candDeltas, axis=0)
    candMAD = np.median(np.abs(candDeltas - candMedian), axis=0)
    candSigma = np.maximum(1.4826 * candMAD, np.float32(sigmaFloorADU))
    candExcess = (maxDelta[ys, xs] - candMedian).astype(np.float32)
    excessFull[ys, xs] = candExcess

    if returnDiagnostics:
        medianFull[ys, xs] = candMedian.astype(np.float32)
        sigmaFull[ys, xs] = candSigma.astype(np.float32)

    candFlag = (candExcess > nSigma * candSigma) & (candExcess > excessFloorADU)
    flag[ys[candFlag], xs[candFlag]] = True
    nFlagged = int(candFlag.sum())

    if nFlagged > 0:
        ysF = ys[candFlag]
        xsF = xs[candFlag]
        ksF = kCR[ysF, xsF]
        excF = candExcess[candFlag].astype(flux.dtype)
        for r in range(nReads):
            sel = ksF < r
            if sel.any():
                flux[r, ysF[sel], xsF[sel]] -= excF[sel]

    return CRRepairResult(
        nFlagged=nFlagged,
        flagMask=flag,
        crReadIdx=kCR,
        excess=excessFull,
        candidateMask=candidateMask if returnDiagnostics else None,
        medianDelta=medianFull,
        sigma=sigmaFull,
    )
