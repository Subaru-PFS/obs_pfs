"""Post-linearization CR + ASIC-glitch detection and repair on H4 ramps.

Iteratively detect single-read positive rate outliers (cosmic rays) and
matched up/down delta pairs (ASIC digital glitches) in a linearized
cumulative ramp. Repair them by replacing the affected deltas with the
per-pixel UTR rate and re-cumsumming (pair-aware so post-glitch reads
are preserved), and report per-delta CR / glitch masks. Replaces all
prior H4 CR-correction paths.

The pure-numpy algorithm is here; the wrapper that integrates with
``PfsIsrTask`` lives in ``lsst.obs.pfs.isrTask``.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np


DEFAULT_SIGMA_FLOOR_ADU = 8.0

# Defaults for iterativeUtrDetectAndRepair.
DEFAULT_ITER_N_SIGMA = 5.0
DEFAULT_MAX_ITERATIONS = 5


@dataclass
class IterativeRepairResult:
    """Outcome of `iterativeUtrDetectAndRepair`.

    Attributes
    ----------
    nIterations : int
        Number of iterations actually run.
    crFlagMask : np.ndarray
        Boolean ``(N-1, H, W)``; True at delta positions identified as
        single-read CRs. The "delta" axis corresponds to ``np.diff(cube,
        axis=0)``; a flag at index ``k`` repaired the delta from read
        ``k`` to read ``k+1``.
    glitchFlagMask : np.ndarray
        Boolean ``(N-1, H, W)``; True at delta positions identified as
        ASIC glitches. Glitches come in pairs (one positive, one negative)
        so this mask is True at BOTH delta indices of each pair.
    rate : np.ndarray
        ``(H, W)`` float32; the per-pixel UTR rate — the mean of the
        un-flagged deltas (CR/glitch-flagged deltas excluded).
    sigma : np.ndarray
        ``(H, W)`` float32; per-delta noise estimate used for thresholding.
    nCRs : int
        Count of delta positions flagged as CR.
    nGlitchPairs : int
        Count of distinct glitch pairs (= nGlitchFlags / 2).
    nByIteration : list
        Per-iteration ``(newCRs, newGlitchPairs)`` tuples for inspection.
    iterationTimings : list
        Wall-clock seconds per iteration. Useful for profiling on real
        (4096^2) data.
    """

    nIterations: int
    crFlagMask: np.ndarray
    glitchFlagMask: np.ndarray
    rate: np.ndarray
    sigma: np.ndarray
    nCRs: int
    nGlitchPairs: int
    nByIteration: list
    iterationTimings: list = field(default_factory=list)


def _rampQR(deltas: np.ndarray, hChunk: int = 512) -> tuple:
    """Per-pixel 25th/50th/75th percentile of ``deltas`` along axis 0.

    Equivalent to ``np.percentile(deltas, [25, 50, 75], axis=0)`` with
    method='linear'.

    Implementation note — why we transpose chunks to ``(h_chunk, W, N)``
    before partitioning:

    Naive ``np.partition(deltas, kth, axis=0)`` on a C-contiguous
    ``(N, H, W)`` cube sorts each pixel's N-tuple along the slowest-
    varying axis. Per-pixel access walks reads strided ``H*W*4`` bytes
    apart, which is cache-hostile. Worse, chunking along axis 1 didn't
    help on this layout — measured 3.4× slowdown vs unchunked, since
    each chunked sub-partition still had the same per-pixel stride.

    Copying each chunk to a contiguous ``(h_chunk, W, N)`` transposed
    buffer puts each pixel's N reads in contiguous memory; numpy's
    in-place ``ndarray.partition`` on axis=-1 then runs ~30% faster
    than the original axis=0 path AND chunks linearly. The transposed
    buffer is partitioned in place, so the per-chunk peak is one
    ``(h_chunk, W, N)`` allocation (~730 MB at hChunk=512) rather than
    the 5.85 GB unchunked ``dp``.

    The 50th percentile (median) is returned as float32 because the
    rest of the CR detector consumes it as a per-pixel rate in float32.
    ``p25``/``p75`` retain ``deltas``'s dtype.
    """
    N, H, W = deltas.shape

    def _ranks(q):
        r = (N - 1) * q
        lo = int(np.floor(r))
        hi = int(np.ceil(r))
        return lo, hi, r - lo

    lo25, hi25, f25 = _ranks(0.25)
    lo50, hi50, f50 = _ranks(0.50)
    lo75, hi75, f75 = _ranks(0.75)
    kth = sorted({lo25, hi25, lo50, hi50, lo75, hi75})

    p25 = np.empty((H, W), dtype=deltas.dtype)
    median = np.empty((H, W), dtype=np.float32)
    p75 = np.empty((H, W), dtype=deltas.dtype)

    # Subset / column-vector paths (H == 1) — one iteration, equivalent
    # to the unchunked behavior.
    for h0 in range(0, H, hChunk):
        h1 = min(h0 + hChunk, H)
        # Force a fresh contiguous copy of the transposed chunk. We
        # CANNOT use ``np.ascontiguousarray`` here: when the chunk has a
        # single-element axis (the subset path's ``(N-1, 1, M)`` becomes
        # ``(1, M, N-1)`` after transpose), numpy considers the view
        # "trivially C-contiguous" and ascontiguousarray returns the
        # input view unchanged. The in-place ``chunkT.partition`` would
        # then mutate the caller's ``deltas`` array. ``.copy()`` always
        # allocates a fresh buffer.
        chunkT = deltas[:, h0:h1, :].transpose(1, 2, 0).copy()
        # Partition IN PLACE along axis=-1; no separate dp buffer.
        chunkT.partition(kth, axis=-1)

        def _interp(lo, hi, frac):
            if lo == hi:
                return chunkT[..., lo]
            return (1.0 - frac) * chunkT[..., lo] + frac * chunkT[..., hi]

        p25[h0:h1, :] = _interp(lo25, hi25, f25)
        median[h0:h1, :] = _interp(lo50, hi50, f50).astype(
            np.float32, copy=False
        )
        p75[h0:h1, :] = _interp(lo75, hi75, f75)
        del chunkT  # free the ~730 MB transpose buffer before next iter

    return p25, median, p75


def _utrRateSimple(cube: np.ndarray) -> np.ndarray:
    """Robust per-pixel rate via median of cumulative deltas (ADU/read).

    Median has a 50% breakdown point — a few CR/glitch hits in a ramp
    don't bias the rate. LSQ slope (the previous implementation) is non-
    robust: a single outlier near the center of a short ramp can drag
    the slope by ~5% of the outlier amplitude, which feeds back into the
    iterative detection loop and causes the rate to diverge runaway
    (e.g. on a 22-read half-ramp with a 6000-ADU glitch, LSQ slope ~ -38
    in iter 1, propagating to ~ -1200 by iter 5; most normal deltas get
    false-flagged because residuals end up at or beyond threshold).
    For a clean ramp the median of deltas matches the mean (= LSQ slope)
    up to sampling noise.
    """
    deltas = np.diff(cube, axis=0)
    return np.median(deltas, axis=0).astype(np.float32, copy=False)


def _detectAndRepairOnce(deltas, goodPixelMask, glitchActive,
                         crAccum, glitchAccum, boundaryAccum,
                         sigmaFloorADU, nSigma, repair, correctGlitches,
                         glitchAmplitudeMinADU=0.0):
    """Run one IQR-sigma detection/repair iteration on a delta cube.

    Modifies ``deltas``, ``crAccum``, ``glitchAccum``, ``boundaryAccum``
    in place. ``crAccum`` collects CRs, ``glitchAccum`` interior glitch
    pairs, ``boundaryAccum`` end glitches (lone flagged delta at index 0
    or N-2). Accepts any ``(Nd, H', W')`` shape, so it works for both the
    full delta cube (H', W' = H, W) and a column-vector subset
    (H'=1, W'=nActive).

    The simple in-iter repair sets flagged deltas to the per-pixel rate.
    Pair-aware repair (which preserves the original cumulative value at
    pair-second positions) is the caller's responsibility — done once at
    the end of :func:`iterativeUtrDetectAndRepair` from the saved
    original deltas.

    ``glitchAmplitudeMinADU`` adds a minimum-amplitude floor for glitch pair
    classification: at least one of the two deltas must have ``|residual|``
    above this value. Use it to suppress faint-end deglitching where the
    classifier is less reliable. 0 disables the extra floor.

    Returns ``(rate, sigma, newCR, newGlitchPairs)`` — counts are *new this
    call* against ``crAccum`` / ``glitchAccum`` at entry.
    """
    # IQR percentiles + median in a single partition pass (see _rampQR).
    # Done BEFORE converting deltas to residual so we can still read the
    # raw delta values.
    p25, rate, p75 = _rampQR(deltas)

    iqrSigma = 0.741 * (p75 - p25)
    sigma = np.maximum(iqrSigma, sigmaFloorADU).astype(np.float32, copy=False)
    threshold = (nSigma * sigma).astype(np.float32, copy=False)

    # Convert deltas → residual IN PLACE. From here ``residual`` aliases
    # the deltas buffer; the original delta values are no longer reachable
    # except by adding ``rate[None]`` back. Saves the ~5.85 GB second
    # (N-1, H, W) cube the previous code held alive concurrently.
    deltas -= rate[None, :, :]
    residual = deltas

    # Pre-filter to narrow detection work to "candidate" pixels — those
    # with at least one delta exceeding threshold. The full pair-detect
    # bool/float expressions are then evaluated on the gathered subset
    # ((N-1, n_cand) instead of (N-1, H, W)), avoiding ~10 GB of
    # transient float buffers and ~30 s of CPU on a 4096²×88 ramp.
    # ``rate``, ``sigma``, and ``threshold`` are still per-pixel from
    # _rampQR, so the detection criterion is bit-identical to the
    # full-cube version — only the *where* changes, not the *what*.

    # Slice-wise candidate scan: build the (H, W) any-delta-exceeds-
    # threshold mask without materializing the full (N-1, H, W) flagged
    # cube (saves ~1.5 GB bool + ~5.85 GB float abs transient).
    cand2D = np.zeros(residual.shape[1:], dtype=bool)
    absSlice = np.empty(residual.shape[1:], dtype=residual.dtype)
    for k in range(residual.shape[0]):
        np.abs(residual[k], out=absSlice)
        cand2D |= absSlice > threshold
    cand2D &= goodPixelMask
    del absSlice

    candYs, candXs = np.where(cand2D)
    nCand = candYs.size

    newCR = 0
    newGlitch = 0

    if nCand > 0:
        # Gather candidate data into compact (N-1, nCand) arrays. These
        # are O(nCand) in size — for the production dark, ~30 K candidates
        # out of 16.7 M pixels.
        residC = residual[:, candYs, candXs]                  # (N-1, nCand)
        threshC = threshold[candYs, candXs]                   # (nCand,)
        gActiveC = glitchActive[candYs, candXs]               # (nCand,)
        flaggedC = np.abs(residC) > threshC[None]             # (N-1, nCand)

        # Pair criterion: at least ONE adjacent delta crosses threshold,
        # the two residuals have opposite signs, and they cancel within
        # threshold. The cancellation test is what discriminates glitch
        # from CR — a real CR leaves residual[k+1] ~ 0, so the sum
        # |resid[k] + resid[k+1]| ~ resid[k] is too large to pass.

        # Opposite-sign via bool XOR (avoids a float product temp).
        negPrev = residC[:-1] < 0
        negNext = residC[1:] < 0
        oppositeSign = negPrev != negNext
        del negPrev, negNext

        # Cancellation: |residual[:-1] + residual[1:]| < threshold.
        # abs in place on the sum.
        sumResid = residC[:-1] + residC[1:]
        np.abs(sumResid, out=sumResid)
        cancels = sumResid < threshC[None]
        del sumResid

        flaggedEither = flaggedC[:-1] | flaggedC[1:]

        pairMatch = oppositeSign & cancels & flaggedEither & gActiveC[None]
        del oppositeSign, cancels, flaggedEither

        if glitchAmplitudeMinADU > 0.0:
            # Extra floor: at least one delta in the pair clears the
            # specified amplitude. Suppresses faint-end glitch
            # classification.
            absT = np.abs(residC[:-1])
            bigEnough = absT > glitchAmplitudeMinADU
            del absT
            absT = np.abs(residC[1:])
            bigEnough |= absT > glitchAmplitudeMinADU
            del absT
            pairMatch &= bigEnough
            del bigEnough

        isPairC = np.zeros_like(flaggedC)
        isPairC[:-1] |= pairMatch
        isPairC[1:] |= pairMatch
        del pairMatch

        # Boundary ("end") glitch heuristic: a flagged delta at index 0
        # or N-2 can never pair because one of its neighbors lies off the
        # ramp. It is classified as an ASIC glitch (so it is not
        # misclassified as a CR) but tracked separately from interior
        # pairs — a lone end glitch has no partner to cancel against, so
        # it must always be corrected.
        isBoundaryC = np.zeros_like(flaggedC)
        isBoundaryC[0] = flaggedC[0] & ~isPairC[0] & gActiveC
        isBoundaryC[-1] = flaggedC[-1] & ~isPairC[-1] & gActiveC

        isCRC = flaggedC & ~isPairC & ~isBoundaryC & (residC > 0)
        del flaggedC

        # Read existing accumulator values at candidate positions; count
        # new flags; scatter the OR'd values back. (Fancy indexing on
        # the LHS produces a copy, so we can't use ``|=`` directly.)
        crAccumC = crAccum[:, candYs, candXs]
        glitchAccumC = glitchAccum[:, candYs, candXs]
        boundaryAccumC = boundaryAccum[:, candYs, candXs]
        newCR = int((isCRC & ~crAccumC).sum())
        newGlitch = int((isPairC & ~glitchAccumC).sum()) // 2
        crAccum[:, candYs, candXs] = crAccumC | isCRC
        glitchAccum[:, candYs, candXs] = glitchAccumC | isPairC
        boundaryAccum[:, candYs, candXs] = boundaryAccumC | isBoundaryC
        del isCRC, isPairC, isBoundaryC, crAccumC, glitchAccumC, boundaryAccumC

    if repair:
        # In-iter repair: flagged residuals → 0 (they become ``rate``
        # after add-back) so the next iteration sees a cleaner sample.
        # CRs and end glitches are always repaired — an end glitch has no
        # pair partner, so it cannot cancel itself out. Interior glitch
        # pairs are repaired only when ``correctGlitches`` is set;
        # otherwise the symmetric +A/-A pair is left in place to cancel
        # in the mean UTR rate on its own.
        allFlag = crAccum | boundaryAccum
        if correctGlitches:
            allFlag = allFlag | glitchAccum
        if allFlag.any():
            residual[allFlag] = 0.0
        del allFlag
    # Always restore delta space (even when nothing was flagged this
    # iter — keeps the in/out contract: input and output both live in
    # delta space).
    residual += rate[None, :, :]

    return rate, sigma, newCR, newGlitch


def iterativeUtrDetectAndRepair(
    deltas: np.ndarray,
    *,
    goodPixelMask: np.ndarray,
    glitchPixelMask: Optional[np.ndarray] = None,
    sigmaFloorADU: float = DEFAULT_SIGMA_FLOOR_ADU,
    nSigma: float = DEFAULT_ITER_N_SIGMA,
    maxIterations: int = DEFAULT_MAX_ITERATIONS,
    repair: bool = True,
    correctGlitches: bool = True,
    glitchAmplitudeMinADU: float = 0.0,
) -> IterativeRepairResult:
    """Iterative CR + ASIC-glitch detection on a linearized delta cube.

    Operates in delta-space throughout — no ``np.diff`` or ``np.cumsum``
    happens inside this function. The caller is responsible for:

      1. ``deltas = np.diff(flux, axis=0)`` ONCE before this call, on the
         linearized cumulative ramp.
      2. Cumulative reconstruction (``np.cumsum``) ONCE after, if needed.

    For linearized data, deltas have uncorrelated read noise and the
    per-pixel mean is the optimal (BLUE) UTR rate estimator; ``result.rate``
    is that mean taken over the un-flagged deltas (CR/glitch-flagged
    deltas excluded). Detection still uses the robust median + IQR
    (computed in a single ``np.partition`` pass via :func:`_rampQR`), so
    iter-1's threshold doesn't get inflated by contamination.

    Each iteration:

      1. Robust median rate + IQR sigma from current deltas (via _rampQR).
      2. Flag deltas with ``|delta - rate| > nSigma * max(sigma,
         sigmaFloorADU)``.
      3. Classify each flagged delta:

         - **Glitch pair (ASIC)**: at least one of two adjacent deltas
           is flagged, the two residuals have opposite signs, and
           ``|resid[k]+resid[k+1]| < nSigma*sigma`` (the recovery
           cancels). Both deltas are marked.
         - **CR (single read)**: a positive flagged delta whose neighbor
           is not flagged as a pair partner.

      4. Simple in-iter repair (flagged → rate) so the next iteration
         sees a cleaner sample. Pair-aware repair happens once at the
         end from the saved original deltas.
      5. Stop when no new flags or rate change is small.

    The detection threshold is the noise floor: glitches/CRs much
    larger than ``nSigma * sigma`` are reliably detected; those near or
    below the noise are not (and shouldn't be — they're indistinguishable
    from per-read noise).

    Parameters
    ----------
    deltas : np.ndarray
        ``(N-1, H, W)`` linearized delta cube; ``deltas[k] = flux[k+1]
        - flux[k]`` for the linearized cumulative cube. Modified in
        place if ``repair=True``: on return, holds the FINAL repaired
        deltas (pair-aware).
    goodPixelMask : np.ndarray
        ``(H, W)`` bool. Pixels where this is False are skipped (no
        flags raised).
    glitchPixelMask : np.ndarray, optional
        ``(H, W)`` bool. ASIC-glitch pair detection runs only where this
        is True. Default ``None`` disables glitch detection entirely.
        Production callers should build the mask from
        ``PfsIsrTask.asicBadChannelMask(...)``.
    sigmaFloorADU : float
        Lower bound on per-pixel IQR sigma.
    nSigma : float
        Detection threshold in sigma units. Defaults to 5.
    maxIterations : int
        Hard cap on iteration count.
    repair : bool
        If False, only flag; deltas are not modified.
    correctGlitches : bool
        Controls correction of *interior* ASIC-glitch pairs. When True
        (default), interior glitch pairs are repaired in the cube and
        excluded from ``result.rate`` alongside CRs. When False, interior
        pairs are detected (so a glitch up-spike is not misclassified as
        a CR) but left in place — kept in the rate mean, where a
        symmetric +A/-A pair cancels on its own. *End glitches* (a lone
        flagged delta at the first or last delta, with no pair partner)
        are always repaired regardless of this flag: with no partner they
        cannot self-cancel.
    glitchAmplitudeMinADU : float
        Minimum |residual| amplitude (ADU) for ASIC-glitch pair
        classification. 0 (default) uses only the CR threshold.

    Returns
    -------
    IterativeRepairResult
        ``rate`` is the per-pixel UTR rate — the mean of the un-flagged
        deltas (CR-flagged always excluded; glitch-flagged excluded only
        when ``correctGlitches``), ready for use as the science rate.
        ``sigma`` is the IQR-based per-pixel scatter from the last
        iteration.
    """
    if deltas.ndim != 3:
        raise ValueError(
            f"deltas must be 3-D (N-1, H, W); got {deltas.shape}"
        )
    nDeltas, H, W = deltas.shape
    if nDeltas < 2:
        raise ValueError(f"deltas must have at least 2 entries; got {nDeltas}.")

    if glitchPixelMask is None:
        # Off by default: glitch pair detection finds nothing.
        glitchActive2D = np.zeros((H, W), dtype=bool)
    else:
        if glitchPixelMask.shape != (H, W):
            raise ValueError(
                f"glitchPixelMask shape {glitchPixelMask.shape} != deltas H,W ({H}, {W})."
            )
        glitchActive2D = np.asarray(glitchPixelMask, dtype=bool)

    crFlagAccum = np.zeros((nDeltas, H, W), dtype=bool)
    glitchFlagAccum = np.zeros((nDeltas, H, W), dtype=bool)
    boundaryFlagAccum = np.zeros((nDeltas, H, W), dtype=bool)
    nByIter = []
    iterTimings = []
    iteration = 0
    rateTolerance = 0.05  # ADU/read; loop stops when rate change < this

    # Pair-aware repair (preserving the original cumulative value at the
    # post-glitch read) is unnecessary on real H4 data: ASIC glitches are
    # symmetric digital bit-flips (e.g., +2048 then -2048), so simple
    # "both deltas → rate" replacement gives the same cumulative ramp
    # downstream of the glitch. Median + IQR statistics are robust
    # enough that the few glitch outliers per pixel don't bias the rate
    # in iter 1 either. So we don't save origDeltas — we just iterate
    # the simple in-iter repair to convergence.

    # ---- Iter 1: process the full delta cube. ----
    tIter0 = time.time()
    rateFull, sigmaFull, newCR, newGlitch = _detectAndRepairOnce(
        deltas, goodPixelMask, glitchActive2D,
        crFlagAccum, glitchFlagAccum, boundaryFlagAccum,
        sigmaFloorADU, nSigma, repair, correctGlitches,
        glitchAmplitudeMinADU=glitchAmplitudeMinADU,
    )
    iterTimings.append(time.time() - tIter0)
    nByIter.append((newCR, newGlitch))

    # ---- Iter 2..N: restrict to the subset of pixels that got any flag in
    #      iter 1. Pixels with no iter-1 flag never see new flags in later
    #      iterations because their deltas are unchanged. ----
    activeMask = (crFlagAccum.any(axis=0) | glitchFlagAccum.any(axis=0)
                  | boundaryFlagAccum.any(axis=0))
    ratePrev = rateFull

    for iteration in range(1, maxIterations):
        ys, xs = np.where(activeMask)
        if ys.size == 0:
            break

        tIter0 = time.time()
        # Subset views are copies (fancy indexing); modify, then write back.
        deltasSub = deltas[:, ys, xs][:, np.newaxis, :].copy()
        crAccumSub = crFlagAccum[:, ys, xs][:, np.newaxis, :].copy()
        glitchAccumSub = glitchFlagAccum[:, ys, xs][:, np.newaxis, :].copy()
        boundaryAccumSub = boundaryFlagAccum[:, ys, xs][:, np.newaxis, :].copy()
        goodSub = goodPixelMask[ys, xs][np.newaxis, :]
        glitchSub = glitchActive2D[ys, xs][np.newaxis, :]

        rateSub, sigmaSub, newCR, newGlitch = _detectAndRepairOnce(
            deltasSub, goodSub, glitchSub,
            crAccumSub, glitchAccumSub, boundaryAccumSub,
            sigmaFloorADU, nSigma, repair, correctGlitches,
            glitchAmplitudeMinADU=glitchAmplitudeMinADU,
        )

        deltas[:, ys, xs] = deltasSub[:, 0, :]
        crFlagAccum[:, ys, xs] = crAccumSub[:, 0, :]
        glitchFlagAccum[:, ys, xs] = glitchAccumSub[:, 0, :]
        boundaryFlagAccum[:, ys, xs] = boundaryAccumSub[:, 0, :]

        rateFull = ratePrev.copy()
        rateFull[ys, xs] = rateSub[0, :]
        sigmaFull[ys, xs] = sigmaSub[0, :]

        iterTimings.append(time.time() - tIter0)
        nByIter.append((newCR, newGlitch))

        # Converge when the per-pixel rate stops moving on the active set.
        rateChange = float(np.max(np.abs(rateFull - ratePrev)))
        if rateChange < rateTolerance:
            break
        ratePrev = rateFull

    # result.rate is the per-pixel UTR rate: the mean of the un-flagged
    # deltas. CR-flagged and end-glitch-flagged deltas are always
    # excluded (both are repaired in the cube). Interior glitch pairs are
    # excluded only when ``correctGlitches`` is set; otherwise they stay
    # in the mean, where a symmetric +A/-A ASIC pair cancels on its own.
    # Excluding (rather than the median or a truncated mean) gives the
    # optimal uniform-weight UTR estimator. The median ``rateFull`` is
    # kept only for detection; a pixel with every delta flagged falls
    # back to it. Accumulated slice-wise to avoid an (N-1, H, W) transient.
    sumUnflagged = np.zeros(deltas.shape[1:], dtype=np.float32)
    nUnflagged = np.zeros(deltas.shape[1:], dtype=np.int32)
    for k in range(nDeltas):
        flagK = crFlagAccum[k] | boundaryFlagAccum[k]
        if correctGlitches:
            flagK = flagK | glitchFlagAccum[k]
        unflagK = ~flagK
        sumUnflagged += deltas[k] * unflagK
        nUnflagged += unflagK
    with np.errstate(invalid="ignore", divide="ignore"):
        rateFinal = (sumUnflagged / nUnflagged.astype(np.float32)).astype(
            np.float32, copy=False)
    allFlagged = nUnflagged == 0
    if allFlagged.any():
        rateFinal[allFlagged] = rateFull[allFlagged]
    sigmaFinal = sigmaFull.astype(np.float32, copy=False)

    # glitchFlagMask reports every ASIC glitch — interior pairs and end
    # glitches — for the ASIC_GLITCH mask plane. nGlitchPairs counts only
    # true pairs (from glitchFlagAccum); end glitches are singletons.
    return IterativeRepairResult(
        nIterations=iteration + 1,
        crFlagMask=crFlagAccum,
        glitchFlagMask=glitchFlagAccum | boundaryFlagAccum,
        rate=rateFinal,
        sigma=sigmaFinal,
        nCRs=int(crFlagAccum.sum()),
        nGlitchPairs=int(glitchFlagAccum.sum()) // 2,
        nByIteration=nByIter,
        iterationTimings=iterTimings,
    )

