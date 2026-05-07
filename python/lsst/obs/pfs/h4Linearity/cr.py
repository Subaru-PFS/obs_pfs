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
        flagged read. Only meaningful where ``flagMask`` is True.
    excess : np.ndarray
        Float ``(H, W)``; ``max_delta - median_delta`` (i.e., the amount
        subtracted from cumulative reads ``>= k+1``). Only meaningful
        where ``flagMask`` is True.
    """

    nFlagged: int
    flagMask: np.ndarray
    crReadIdx: np.ndarray
    excess: np.ndarray


def detectAndRepair(
    flux: np.ndarray,
    *,
    goodPixelMask: np.ndarray,
    nSigma: float = DEFAULT_N_SIGMA,
    sigmaFloorADU: float = DEFAULT_SIGMA_FLOOR_ADU,
    excessFloorADU: float = DEFAULT_EXCESS_FLOOR_ADU,
    minReads: int = DEFAULT_MIN_READS,
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
    empty = CRRepairResult(
        nFlagged=0,
        flagMask=np.zeros((H, W), dtype=bool),
        crReadIdx=np.zeros((H, W), dtype=np.int32),
        excess=np.zeros((H, W), dtype=np.float32),
    )
    if nReads < minReads:
        return empty

    deltas = np.diff(flux, axis=0)
    median = np.median(deltas, axis=0)
    mad = np.median(np.abs(deltas - median), axis=0)
    sigma = np.maximum(1.4826 * mad, np.float32(sigmaFloorADU))

    kCR = np.argmax(deltas, axis=0).astype(np.int32)
    maxDelta = np.take_along_axis(deltas, kCR[None], axis=0)[0]
    excess = (maxDelta - median).astype(np.float32)

    flag = (excess > nSigma * sigma) & (excess > excessFloorADU) & goodPixelMask

    nFlagged = int(flag.sum())
    if nFlagged == 0:
        return empty

    ys, xs = np.where(flag)
    ks = kCR[ys, xs]
    exc = excess[ys, xs].astype(flux.dtype)
    for r in range(nReads):
        sel = ks < r
        if sel.any():
            flux[r, ys[sel], xs[sel]] -= exc[sel]

    return CRRepairResult(
        nFlagged=nFlagged,
        flagMask=flag,
        crReadIdx=kCR,
        excess=excess,
    )
