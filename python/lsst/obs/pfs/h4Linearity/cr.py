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
        ``(H, W)`` float32; the final per-pixel UTR rate after repair.
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


def _detectAndRepairOnce(cube, goodPixelMask, glitchActive,
                         crAccum, glitchAccum,
                         sigmaFloorADU, nSigma, repair,
                         glitchAmplitudeMinADU=0.0):
    """Run one IQR-sigma detection/repair iteration.

    Modifies ``cube``, ``crAccum``, ``glitchAccum`` in place. Accepts any
    ``(N, H', W')`` shape, so it works for both the full cube (H', W'=H, W)
    and a column-vector subset (H'=1, W'=nActive).

    ``glitchAmplitudeMinADU`` adds a minimum-amplitude floor for glitch pair
    classification: at least one of the two deltas must have ``|residual|``
    above this value. Use it to suppress faint-end deglitching where the
    classifier is less reliable. 0 disables the extra floor.

    Returns ``(rate, sigma, newCR, newGlitchPairs)`` — counts are *new this
    call* against ``crAccum`` / ``glitchAccum`` at entry.
    """
    # Compute deltas once and reuse for the rate, sigma, and threshold
    # passes. Avoid calling `_utrRateSimple(cube)` here — it would do its
    # own np.diff and discard it, doubling the largest allocation in the
    # loop.
    deltas = np.diff(cube, axis=0)
    rate = np.median(deltas, axis=0).astype(np.float32, copy=False)
    residual = deltas - rate[None, :, :]

    p25, p75 = np.percentile(deltas, [25, 75], axis=0)
    iqrSigma = 0.741 * (p75 - p25)
    sigma = np.maximum(iqrSigma, sigmaFloorADU).astype(np.float32, copy=False)
    threshold = (nSigma * sigma).astype(np.float32, copy=False)

    flagged = np.abs(residual) > threshold[None, :, :]
    flagged &= goodPixelMask[None, :, :]

    # Pair criterion: at least ONE adjacent delta crosses threshold, the two
    # residuals have opposite signs, and they cancel within threshold. The
    # cancellation test is what discriminates glitch from CR — a real CR
    # leaves residual[k+1] ~ 0, so |resid[k] + resid[k+1]| ~ resid[k] is too
    # large to pass. Requiring both deltas independently flagged misses
    # asymmetric glitches where the recovery sits just under threshold.
    pairMatch = (
        (flagged[:-1] | flagged[1:])
        & (residual[:-1] * residual[1:] < 0)
        & (np.abs(residual[:-1] + residual[1:]) < threshold[None])
        & glitchActive[None]
    )
    if glitchAmplitudeMinADU > 0.0:
        # Extra floor: require at least one delta in the pair to clear the
        # specified amplitude. Suppresses faint-end glitch classification.
        bigEnough = (np.abs(residual[:-1]) > glitchAmplitudeMinADU) | (
            np.abs(residual[1:]) > glitchAmplitudeMinADU
        )
        pairMatch &= bigEnough
    isPair = np.zeros_like(flagged)
    isPair[:-1] |= pairMatch
    isPair[1:] |= pairMatch

    # Boundary heuristic: a flagged delta at index 0 or N-2 can never pair
    # because one of its neighbors lies off the ramp. Empirically these
    # produce digital-glitch-shaped residuals (clusters at ±4096, ±8192...)
    # — mark the affected pixel as ASIC_GLITCH (joins glitchAccum) so it
    # gets the same repair path as a normal pair and isn't misclassified
    # as a CR. Only fires where glitch detection is enabled.
    isBoundary = np.zeros_like(flagged)
    isBoundary[0] = flagged[0] & ~isPair[0] & glitchActive[None]
    isBoundary[-1] = flagged[-1] & ~isPair[-1] & glitchActive[None]
    isPair |= isBoundary

    isCR = flagged & ~isPair & (residual > 0)

    newCR = int((isCR & ~crAccum).sum())
    newGlitch = int((isPair & ~glitchAccum).sum()) // 2

    crAccum |= isCR
    glitchAccum |= isPair

    if repair:
        allFlag = crAccum | glitchAccum
        if allFlag.any():
            rate3D = np.broadcast_to(rate[None, :, :], deltas.shape)
            deltasRepaired = np.where(
                allFlag, rate3D, deltas
            ).astype(cube.dtype, copy=False)
            cube[1:] = cube[0:1] + np.cumsum(deltasRepaired, axis=0)

    return rate, sigma, newCR, newGlitch


def iterativeUtrDetectAndRepair(
    cube: np.ndarray,
    *,
    goodPixelMask: np.ndarray,
    glitchPixelMask: Optional[np.ndarray] = None,
    sigmaFloorADU: float = DEFAULT_SIGMA_FLOOR_ADU,
    nSigma: float = DEFAULT_ITER_N_SIGMA,
    maxIterations: int = DEFAULT_MAX_ITERATIONS,
    repair: bool = True,
    glitchAmplitudeMinADU: float = 0.0,
) -> IterativeRepairResult:
    """Iterative UTR-rate-based detection of CRs and ASIC glitches.

    Each iteration:

      1. Estimate per-pixel rate via the UTR (LSQ slope) of the
         current cube.
      2. Compute per-read residual ``delta - rate``.
      3. Estimate per-pixel sigma as ``0.741 * (p75 − p25)`` of the
         per-read deltas (axis=0 IQR), floored at ``sigmaFloorADU``.
         This adapts to each pixel's actual scatter — hot / RTS pixels
         get a higher sigma than the noise-model assumption would give
         — and is robust to ~25% outlier contamination per pixel.
      4. Flag deltas with ``|residual| > nSigma * sigma``.
      5. Classify each flagged delta:

         - **Glitch pair (ASIC)**: at least one of two adjacent deltas
           is flagged, the two residuals have opposite signs, and
           ``|resid[k]+resid[k+1]| < nSigma*sigma`` (the recovery
           cancels). Both deltas are marked. Requiring only one
           independently flagged catches asymmetric pairs where the
           recovery delta sits just under threshold; the sum-cancellation
           is what discriminates glitch from CR.
         - **CR (single read)**: a positive flagged delta whose neighbor
           is not flagged as a pair partner.

      6. Repair by replacing flagged deltas with ``rate`` and
         re-cumsumming. The cube is modified in place.
      7. Stop when no new flags are produced, or after ``maxIterations``.

    The detection threshold is the noise floor of the data: glitches or
    CRs much larger than ``nSigma * sigma`` are reliably detected; those
    near or below the noise are not (and shouldn't be — they're
    indistinguishable from per-read noise).

    Parameters
    ----------
    cube : np.ndarray
        ``(N, H, W)`` linearized cumulative cube. Modified in place if
        ``repair=True``.
    goodPixelMask : np.ndarray
        ``(H, W)`` bool. Pixels where this is False are skipped (no
        flags raised).
    glitchPixelMask : np.ndarray, optional
        ``(H, W)`` bool. ASIC-glitch pair detection runs only where this
        is True. Default ``None`` disables glitch detection entirely
        (the result has ``nGlitchPairs == 0``), so flagged deltas
        outside an opt-in region are classified as single-read CRs
        regardless of their neighbor's residual. Production callers
        should build the mask from ``PfsIsrTask.asicBadChannelMask(...)``
        (which uses ``loadBadAsicChannels`` — ``n3`` channel 24 today,
        128 rows along y per channel).
    sigmaFloorADU : float
        Lower bound on per-pixel IQR sigma. Protects pixels whose IQR
        collapses (e.g. perfectly-quiet dark pixels) from getting an
        unrealistically small threshold.
    nSigma : float
        Detection threshold in sigma units. Defaults to 5.
    maxIterations : int
        Hard cap on iteration count.
    repair : bool
        If False, only flag; cube is not modified. The rate estimate
        will not converge as well without repair, so the flag counts
        from ``repair=False`` are a lower bound on what ``repair=True``
        finds.
    glitchAmplitudeMinADU : float
        Minimum |residual| amplitude (ADU) required for an ASIC-glitch
        pair to be classified. 0 (default) uses only the CR threshold
        ``nSigma * sigma``. Set higher to suppress faint-end deglitching
        where the classifier is less reliable; the condition is that at
        least one delta in the pair has |residual| above this value.

    Returns
    -------
    IterativeRepairResult
    """
    nReads, H, W = cube.shape
    if nReads < 3:
        raise ValueError(f"cube must have at least 3 reads; got {nReads}.")

    if glitchPixelMask is None:
        # Off by default: glitch pair detection finds nothing.
        glitchActive2D = np.zeros((H, W), dtype=bool)
    else:
        if glitchPixelMask.shape != (H, W):
            raise ValueError(
                f"glitchPixelMask shape {glitchPixelMask.shape} != cube H,W ({H}, {W})."
            )
        glitchActive2D = np.asarray(glitchPixelMask, dtype=bool)

    crFlagAccum = np.zeros((nReads - 1, H, W), dtype=bool)
    glitchFlagAccum = np.zeros((nReads - 1, H, W), dtype=bool)
    nByIter = []
    iterTimings = []
    iteration = 0
    rateTolerance = 0.05  # ADU/read; loop stops when rate change < this

    # Save the unmodified cube so the final repair pass can use ORIGINAL
    # deltas — necessary for glitch pairs, where preserving the original
    # cube[k+2] requires knowing the original delta[k+1] amplitude.
    # The naive in-place repair inside _detectAndRepairOnce is still useful
    # for rate convergence across iterations; we just override the final
    # cube state at the end.
    if repair:
        cubeOriginal = cube.copy()

    # ---- Iter 1: process the full cube. ----
    tIter0 = time.time()
    rateFull, sigmaFull, newCR, newGlitch = _detectAndRepairOnce(
        cube, goodPixelMask, glitchActive2D,
        crFlagAccum, glitchFlagAccum,
        sigmaFloorADU, nSigma, repair,
        glitchAmplitudeMinADU=glitchAmplitudeMinADU,
    )
    iterTimings.append(time.time() - tIter0)
    nByIter.append((newCR, newGlitch))

    # ---- Iter 2..N: restrict to the subset of pixels that got any flag in
    #      iter 1. Pixels with no iter-1 flag never see new flags in later
    #      iterations because their cube is unchanged, so their rate and
    #      threshold are unchanged. ----
    activeMask = crFlagAccum.any(axis=0) | glitchFlagAccum.any(axis=0)
    ratePrev = rateFull

    for iteration in range(1, maxIterations):
        ys, xs = np.where(activeMask)
        if ys.size == 0:
            break

        tIter0 = time.time()
        # Subset views are copies (fancy indexing); modify, then write back.
        cubeSub = cube[:, ys, xs][:, np.newaxis, :].copy()
        crAccumSub = crFlagAccum[:, ys, xs][:, np.newaxis, :].copy()
        glitchAccumSub = glitchFlagAccum[:, ys, xs][:, np.newaxis, :].copy()
        goodSub = goodPixelMask[ys, xs][np.newaxis, :]
        glitchSub = glitchActive2D[ys, xs][np.newaxis, :]

        rateSub, sigmaSub, newCR, newGlitch = _detectAndRepairOnce(
            cubeSub, goodSub, glitchSub,
            crAccumSub, glitchAccumSub,
            sigmaFloorADU, nSigma, repair,
            glitchAmplitudeMinADU=glitchAmplitudeMinADU,
        )

        cube[:, ys, xs] = cubeSub[:, 0, :]
        crFlagAccum[:, ys, xs] = crAccumSub[:, 0, :]
        glitchFlagAccum[:, ys, xs] = glitchAccumSub[:, 0, :]

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

    # Final rate / sigma: composed from iter-1's full-cube values (for
    # inactive pixels) and the last subset iteration (for active pixels).
    # No need to recompute IQR on the full cube — that's the most expensive
    # operation in the function (~150 s on the dark).
    rateFinal = rateFull
    sigmaFinal = sigmaFull.astype(np.float32, copy=False)

    # ---- Final pair-aware repair pass ----
    # Rebuild ``cube`` from the saved original using all accumulated flags.
    # For glitch pairs (consecutive glitchFlagAccum flags at k, k+1):
    #   delta[k]_repaired   = rate         (zeros out the spike at read k+1)
    #   delta[k+1]_repaired = origDelta[k] + origDelta[k+1] - rate
    # This makes cube[k+2] match the original (no shift propagating to
    # subsequent reads), while cube[k+1] becomes the smooth interpolation
    # cube[k] + rate. For CRs, delta becomes rate as before.
    if repair:
        origDeltas = np.diff(cubeOriginal, axis=0)
        allFlag = crFlagAccum | glitchFlagAccum
        rate3D = np.broadcast_to(rateFinal[None, :, :], origDeltas.shape)
        deltasRepaired = np.where(
            allFlag, rate3D, origDeltas
        ).astype(cube.dtype, copy=False)

        # pairSecond[k]=True means glitchFlagAccum[k-1] AND glitchFlagAccum[k]
        # are both True — k is the second delta of a pair.
        pairSecond = np.zeros_like(glitchFlagAccum)
        pairSecond[1:] = glitchFlagAccum[:-1] & glitchFlagAccum[1:]
        if pairSecond.any():
            origSum = np.zeros_like(origDeltas)
            origSum[1:] = origDeltas[:-1] + origDeltas[1:]
            deltasRepaired = np.where(
                pairSecond, origSum - rate3D, deltasRepaired
            ).astype(cube.dtype, copy=False)

        cube[1:] = cubeOriginal[0:1] + np.cumsum(deltasRepaired, axis=0)
        del cubeOriginal  # release the ~2.5 GB copy as soon as possible

    return IterativeRepairResult(
        nIterations=iteration + 1,
        crFlagMask=crFlagAccum,
        glitchFlagMask=glitchFlagAccum,
        rate=rateFinal,
        sigma=sigmaFinal,
        nCRs=int(crFlagAccum.sum()),
        nGlitchPairs=int(glitchFlagAccum.sum()) // 2,
        nByIteration=nByIter,
        iterationTimings=iterTimings,
    )

