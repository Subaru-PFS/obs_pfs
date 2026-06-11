"""Post-linearization CR + ASIC-glitch detection and repair on H4 ramps.

Iteratively detect single-read positive rate outliers (cosmic rays) and
matched up/down delta pairs (ASIC digital glitches) in a linearized
cumulative ramp. Repair them by replacing the affected deltas with the
per-pixel UTR rate and re-cumsumming (pair-aware so post-glitch reads
are preserved), and report per-delta CR / glitch masks.

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
        Boolean ``(H, W, N-1)``; True at delta positions identified as
        single-read CRs. The last axis is the delta axis: a flag at
        ``[y, x, k]`` repaired the delta from read ``k`` to read
        ``k+1`` at pixel ``(y, x)``.
    glitchFlagMask : np.ndarray
        Boolean ``(H, W, N-1)``; True at delta positions identified as
        ASIC glitches. Glitches come in pairs (one positive, one negative)
        so this mask is True at BOTH delta indices of each pair.
    unclassifiedFlagMask : np.ndarray
        Boolean ``(H, W, N-1)``; True at delta positions whose absolute
        residual exceeded the detection threshold but did not fall into
        the CR / glitch pair / boundary classes — negative-residual
        outliers and CR candidates demoted by the cumulative-drop check.
        Diagnostic only; not subtracted from ``rate``.
    badPixelMask : np.ndarray
        Boolean ``(H, W)``; True at pixels whose ramp is dominated by
        RTS / telegraph-noise behaviour rather than a single CR event
        (many 3σ delta excursions AND a significantly negative ``rate``).
        Downstream consumers should OR this into the ``BAD`` mask plane
        and ignore any CR/glitch flags at these pixels — their
        classifications are unreliable.
    rate : np.ndarray
        ``(H, W)`` float32; the per-pixel UTR-weighted rate
        (= optimal least-squares slope) computed on the
        CR/glitch-corrected ramp. Equivalent to
        :meth:`PfsIsrTask.calcUTRrates` applied to the reconstructed
        cumulative flux ``read0 + cumsum(deltas)``.
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
    unclassifiedFlagMask: np.ndarray
    badPixelMask: np.ndarray
    rate: np.ndarray
    sigma: np.ndarray
    nCRs: int
    nGlitchPairs: int
    nByIteration: list
    iterationTimings: list = field(default_factory=list)


def _rampQR(deltas: np.ndarray) -> tuple:
    """Per-pixel 25th/50th/75th percentile of ``deltas`` along the time axis.

    Equivalent to ``np.percentile(deltas, [25, 50, 75], axis=-1)`` with
    method='linear', but uses a single ``ndarray.partition`` pass to
    extract the three percentiles together (~3x faster than three
    separate percentile calls).

    ``deltas`` is ``(H, W, N-1)`` — the time axis is last and
    contiguous, so per-pixel partitioning is cache-friendly and works
    directly on the input without copying or transposing. The partition
    is performed on a one-time scratch buffer so the caller's ``deltas``
    is not reordered.

    The 50th percentile (median) is returned as float32 because the
    rest of the CR detector consumes it as a per-pixel rate in float32.
    ``p25``/``p75`` retain ``deltas``'s dtype.
    """
    N = deltas.shape[-1]

    def _ranks(q):
        r = (N - 1) * q
        lo = int(np.floor(r))
        hi = int(np.ceil(r))
        return lo, hi, r - lo

    lo25, hi25, f25 = _ranks(0.25)
    lo50, hi50, f50 = _ranks(0.50)
    lo75, hi75, f75 = _ranks(0.75)
    kth = sorted({lo25, hi25, lo50, hi50, lo75, hi75})

    def _interp(arr, lo, hi, frac):
        """Linear interpolation between ranks ``lo`` and ``hi`` on axis=-1."""
        if lo == hi:
            return arr[..., lo]
        return (1.0 - frac) * arr[..., lo] + frac * arr[..., hi]

    # Partition IN PLACE on a scratch copy so the caller's ``deltas``
    # array is not reordered. The transient copy is the same size as
    # ``deltas`` (~5.85 GB on a 4096²×88 ramp); this is the tradeoff
    # for keeping ``deltas`` immutable here.
    scratch = deltas.copy()
    scratch.partition(kth, axis=-1)

    p25 = _interp(scratch, lo25, hi25, f25)
    median = _interp(scratch, lo50, hi50, f50).astype(np.float32, copy=False)
    p75 = _interp(scratch, lo75, hi75, f75)
    return p25, median, p75


def _utrRateSimple(cube: np.ndarray) -> np.ndarray:
    """Robust per-pixel rate via median of cumulative deltas (ADU/read).

    Median has a 50% breakdown point — a few CR/glitch hits in a ramp
    don't bias the rate. LSQ slope (the obvious alternative) is non-
    robust: a single outlier near the center of a short ramp can drag
    the slope by ~5% of the outlier amplitude, which feeds back into
    the iterative detection loop and causes the rate to diverge. For a
    clean ramp the median of deltas matches the mean (= LSQ slope) up
    to sampling noise.

    ``cube`` is ``(H, W, N)`` cumulative ADU; the time axis is last.
    """
    deltas = np.diff(cube, axis=-1)
    return np.median(deltas, axis=-1).astype(np.float32, copy=False)


def _detectAndRepairOnce(deltas, goodPixelMask, glitchActive,
                         crAccum, glitchAccum, boundaryAccum, unclassAccum,
                         sigmaFloorADU, nSigma, repair, correctGlitches,
                         glitchAmplitudeMinADU=0.0,
                         maxDropFraction=0.5,
                         nDropSigma=3.0):
    """Run one IQR-sigma detection/repair iteration on a delta cube.

    Modifies ``deltas``, ``crAccum``, ``glitchAccum``, ``boundaryAccum``
    in place. ``crAccum`` collects CRs, ``glitchAccum`` interior glitch
    pairs, ``boundaryAccum`` end glitches (lone flagged delta at index 0
    or N-2). Accepts any ``(H', W', Nd)`` shape, so it works for both
    the full delta cube (H', W' = H, W) and a column-vector subset
    (H'=1, W'=nActive). The time axis is always the last axis.

    The simple in-iter repair sets flagged deltas to the per-pixel rate.

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

    # Convert deltas → residual IN PLACE — ``residual`` aliases the
    # deltas buffer; the raw delta values can only be recovered by
    # adding ``rate[..., None]`` back. Done in place to avoid a second
    # (H, W, N-1) cube — ~5.85 GB at 4096²×88.
    deltas -= rate[..., None]
    residual = deltas

    # Pre-filter to narrow detection work to "candidate" pixels — those
    # with at least one delta exceeding threshold. The full pair-detect
    # bool/float expressions are then evaluated on the gathered subset
    # ((nCand, N-1) instead of (H, W, N-1)), avoiding ~10 GB of transient
    # float buffers and ~30 s of CPU on a 4096²×88 ramp. ``rate``,
    # ``sigma``, and ``threshold`` are still per-pixel from _rampQR, so
    # the detection criterion is bit-identical to the full-cube version —
    # only the *where* changes, not the *what*.

    # Slice-wise candidate scan: build the (H, W) any-delta-exceeds-
    # threshold mask without materializing the full (H, W, N-1) flagged
    # cube (saves ~1.5 GB bool + ~5.85 GB float abs transient).
    cand2D = np.zeros(residual.shape[:-1], dtype=bool)
    absSlice = np.empty(residual.shape[:-1], dtype=residual.dtype)
    for k in range(residual.shape[-1]):
        np.abs(residual[..., k], out=absSlice)
        cand2D |= absSlice > threshold
    cand2D &= goodPixelMask
    del absSlice

    candYs, candXs = np.where(cand2D)
    nCand = candYs.size

    newCR = 0
    newGlitch = 0

    if nCand > 0:
        # Gather candidate data into compact (nCand, N-1) arrays. These
        # are O(nCand) in size — for the production dark, ~30 K candidates
        # out of 16.7 M pixels. Gathered along axes 0/1 (spatial) so the
        # last axis (time) stays contiguous per candidate.
        residC = residual[candYs, candXs, :]                  # (nCand, N-1)
        threshC = threshold[candYs, candXs]                   # (nCand,)
        gActiveC = glitchActive[candYs, candXs]               # (nCand,)
        flaggedC = np.abs(residC) > threshC[:, None]          # (nCand, N-1)

        # Pair criterion: at least ONE adjacent delta crosses threshold,
        # the two residuals have opposite signs, and they cancel within
        # threshold. The cancellation test is what discriminates glitch
        # from CR — a real CR leaves residual[k+1] ~ 0, so the sum
        # |resid[k] + resid[k+1]| ~ resid[k] is too large to pass.

        # Opposite-sign via bool XOR (avoids a float product temp).
        negPrev = residC[:, :-1] < 0
        negNext = residC[:, 1:] < 0
        oppositeSign = negPrev != negNext
        del negPrev, negNext

        # Cancellation: |residual[:-1] + residual[1:]| < threshold.
        # abs in place on the sum.
        sumResid = residC[:, :-1] + residC[:, 1:]
        np.abs(sumResid, out=sumResid)
        cancels = sumResid < threshC[:, None]
        del sumResid

        flaggedEither = flaggedC[:, :-1] | flaggedC[:, 1:]

        pairMatch = oppositeSign & cancels & flaggedEither & gActiveC[:, None]
        del oppositeSign, cancels, flaggedEither

        if glitchAmplitudeMinADU > 0.0:
            # Extra floor: at least one delta in the pair clears the
            # specified amplitude. Suppresses faint-end glitch
            # classification.
            absT = np.abs(residC[:, :-1])
            bigEnough = absT > glitchAmplitudeMinADU
            del absT
            absT = np.abs(residC[:, 1:])
            bigEnough |= absT > glitchAmplitudeMinADU
            del absT
            pairMatch &= bigEnough
            del bigEnough

        isPairC = np.zeros_like(flaggedC)
        isPairC[:, :-1] |= pairMatch
        isPairC[:, 1:] |= pairMatch
        del pairMatch

        # Boundary ("end") glitch heuristic: a flagged delta at index 0
        # or N-2 can never pair because one of its neighbors lies off the
        # ramp. It is classified as an ASIC glitch (so it is not
        # misclassified as a CR) but tracked separately from interior
        # pairs — a lone end glitch has no partner to cancel against, so
        # it must always be corrected.
        isBoundaryC = np.zeros_like(flaggedC)
        isBoundaryC[:, 0] = flaggedC[:, 0] & ~isPairC[:, 0] & gActiveC
        isBoundaryC[:, -1] = flaggedC[:, -1] & ~isPairC[:, -1] & gActiveC

        isCRC = flaggedC & ~isPairC & ~isBoundaryC & (residC > 0)

        # Cumulative-drop check: a real cosmic ray deposits charge that
        # stays for the rest of the ramp, so the per-pixel running
        # cumulative residual from the CR delta onward should hover near
        # zero (modulo random-walk noise). A transient up-spike (RTS,
        # decaying glitch, etc.) leaks back and drives the cumulative
        # residual significantly negative. For each candidate CR delta k
        # with amplitude ``A = residC[k] > 0``, compute
        # ``minDropAfter[k] = min over m > k of sum_{j=k+1..m} residC[j]``.
        # We REJECT only when BOTH:
        #
        #   (a) the drop exceeds a fraction of the CR amplitude
        #       (``minDropAfter < -maxDropFraction * A``), and
        #   (b) the drop exceeds the noise-driven cumulative
        #       random-walk excursion
        #       (``minDropAfter < -nDropSigma * sigma * sqrt(N-1-k)``).
        #
        # Both criteria together: real persistent CRs survive because
        # their cumulative drift is noise-limited (small in sigma units);
        # genuine transients fail both because their cumulative drops
        # back ~A which is both >> 0.5*A and >> a few-sigma random walk.
        # (Note: not to be confused with H4RG image persistence — the
        # latent-charge phenomenon — which is a different effect.) CRs
        # at the very last delta have no "after" to check and are
        # accepted. Skip the check entirely when ``maxDropFraction`` is
        # infinite.
        if isCRC.any() and np.isfinite(maxDropFraction):
            forwardCS = np.cumsum(residC, axis=-1)  # (nCand, Nd)
            # runningMinFwd[i, m] = min over j >= m of forwardCS[i, j]
            runningMinFwd = np.minimum.accumulate(
                forwardCS[:, ::-1], axis=-1
            )[:, ::-1]
            minDropAfter = np.zeros_like(forwardCS)
            minDropAfter[:, :-1] = runningMinFwd[:, 1:] - forwardCS[:, :-1]
            del forwardCS, runningMinFwd
            # Criterion (a): relative-to-amplitude.
            relativeOk = minDropAfter > -maxDropFraction * residC
            if np.isfinite(nDropSigma):
                # Criterion (b): noise-aware. nAfter[k] = N-1-k is the
                # number of deltas summed in the cumulative residual at
                # the worst-drop position; its std grows as σ√nAfter.
                nDeltasC = residC.shape[-1]
                nAfter = np.arange(nDeltasC - 1, -1, -1, dtype=np.float32)
                sigmaC = sigma[candYs, candXs].astype(np.float32, copy=False)
                sigmaScale = sigmaC[:, None] * np.sqrt(nAfter)[None, :]
                sigmaOk = minDropAfter > -nDropSigma * sigmaScale
                dropOk = relativeOk | sigmaOk
                del sigmaC, sigmaScale, sigmaOk
            else:
                dropOk = relativeOk
            isCRC &= dropOk
            del minDropAfter, relativeOk, dropOk

        # Outliers above |residual| threshold that the classifier did not
        # assign to any class (not CR, not glitch pair, not boundary).
        # These include negative-residual flagged deltas (RTS bursts,
        # persistence echoes, defects) and CR candidates demoted by the
        # cumulative-drop check. They aren't physically interpretable as
        # ramp signal, so the final rate excludes them.
        isUnclassC = flaggedC & ~isPairC & ~isBoundaryC & ~isCRC
        del flaggedC

        # Read existing accumulator values at candidate positions; count
        # new flags; scatter the OR'd values back. (Fancy indexing on
        # the LHS produces a copy, so we can't use ``|=`` directly.)
        crAccumC = crAccum[candYs, candXs, :]
        glitchAccumC = glitchAccum[candYs, candXs, :]
        boundaryAccumC = boundaryAccum[candYs, candXs, :]
        unclassAccumC = unclassAccum[candYs, candXs, :]
        newCR = int((isCRC & ~crAccumC).sum())
        newGlitch = int((isPairC & ~glitchAccumC).sum()) // 2
        crAccum[candYs, candXs, :] = crAccumC | isCRC
        glitchAccum[candYs, candXs, :] = glitchAccumC | isPairC
        boundaryAccum[candYs, candXs, :] = boundaryAccumC | isBoundaryC
        unclassAccum[candYs, candXs, :] = unclassAccumC | isUnclassC
        del isCRC, isPairC, isBoundaryC, isUnclassC
        del crAccumC, glitchAccumC, boundaryAccumC, unclassAccumC

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
    residual += rate[..., None]

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
    correctGlitches: bool = False,
    glitchAmplitudeMinADU: float = 0.0,
    maxDropFraction: float = 0.5,
    nDropSigma: float = 3.0,
    badPixelMinOutliers: int = 4,
    badPixelOutlierSigma: float = 4.0,
) -> IterativeRepairResult:
    """Iterative CR + ASIC-glitch detection on a linearized delta cube.

    Operates in delta-space throughout — no ``np.diff`` or ``np.cumsum``
    happens inside this function. The caller is responsible for:

      1. ``deltas = np.diff(flux, axis=-1)`` ONCE before this call, on
         a ``(H, W, N)`` linearized cumulative ramp.
      2. Cumulative reconstruction (``np.cumsum(axis=-1)``) ONCE after,
         if needed.

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
         sees a cleaner sample.
      5. Stop when no new flags or rate change is small.

    The detection threshold is the noise floor: glitches/CRs much
    larger than ``nSigma * sigma`` are reliably detected; those near or
    below the noise are not (and shouldn't be — they're indistinguishable
    from per-read noise).

    Parameters
    ----------
    deltas : np.ndarray
        ``(H, W, N-1)`` linearized delta cube with the time axis last;
        ``deltas[y, x, k] = flux[y, x, k+1] - flux[y, x, k]`` for the
        linearized cumulative cube. Modified in place if ``repair=True``:
        on return, holds the FINAL repaired deltas.
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
        Controls correction of *interior* ASIC-glitch pairs. Default
        False — the policy is to detect ASIC glitches but not repair
        them: interior pairs are flagged (so a glitch up-spike is not
        misclassified as a CR) and left in place, where the symmetric
        +A/-A pair cancels on its own in the rate mean. Set True to
        also repair interior pairs in the cube and exclude their deltas
        from ``result.rate`` alongside CRs. *End glitches* (a lone
        flagged delta at the first or last delta, with no pair partner)
        are always repaired regardless of this flag: with no partner
        they cannot self-cancel.
    glitchAmplitudeMinADU : float
        Minimum |residual| amplitude (ADU) for ASIC-glitch pair
        classification. 0 (default) uses only the CR threshold.
    maxDropFraction : float
        Cumulative-drop check for the CR classifier — amplitude
        criterion. A real CR deposits charge that stays for the rest
        of the ramp, so the per-pixel running cumulative residual from
        the CR delta onward should hover near zero. A candidate CR
        with amplitude ``A`` fails this criterion when the cumulative
        residual drops by more than ``maxDropFraction * A`` at any
        later read. Default 0.5. Set to ``float('inf')`` to disable
        the check (alone). (This is unrelated to H4RG image
        persistence, which is a different detector-level phenomenon.)
    nDropSigma : float
        Cumulative-drop check for the CR classifier — noise criterion.
        The cumulative residual after the CR delta is a random walk
        whose std grows as ``σ × √(N-1-k)``; even a true CR can drift
        a couple of sigma negative just from noise. A candidate fails
        this criterion when the cumulative drop exceeds
        ``nDropSigma × σ × √(N-1-k)``. The CR is **rejected** only
        when BOTH criteria fail (drop exceeds the amplitude fraction
        AND exceeds the noise floor), which keeps marginal CRs whose
        noise-driven cumulative drift would otherwise be misread as
        decay. Default 3.0. Set to ``float('inf')`` to disable the
        noise criterion (then the amplitude criterion alone gates
        rejection).
    badPixelMinOutliers : int
        BAD-pixel gate (count criterion). A pixel is marked BAD in
        ``result.badPixelMask`` when it has at least this many delta
        residuals exceeding ``badPixelOutlierSigma × σ_IQR`` from the
        per-pixel median. Counts are taken on the pristine input
        deltas (before in-place repair). Default 4. Set to ``0`` to
        disable the BAD-pixel pass entirely.
    badPixelOutlierSigma : float
        BAD-pixel gate (sigma criterion). Per-delta outlier threshold,
        in units of ``σ_IQR``. Default 4.0 — for clean Gaussian noise
        the expected outlier count per N=80 ramp is ~0.005, so the
        combined ``count ≥ 4 at ≥ 4σ`` gate is essentially
        false-positive-free; real H4 deltas have mild non-Gaussian
        tails (shot/read-noise mixture + linearity residuals) that
        ``σ_IQR`` underestimates by ~30 %, and 3σ at the count=4
        threshold pulls in many merely-noisy pixels.

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
            f"deltas must be 3-D (H, W, N-1); got {deltas.shape}"
        )
    H, W, nDeltas = deltas.shape
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

    crFlagAccum = np.zeros((H, W, nDeltas), dtype=bool)
    glitchFlagAccum = np.zeros((H, W, nDeltas), dtype=bool)
    boundaryFlagAccum = np.zeros((H, W, nDeltas), dtype=bool)
    # Outlier deltas flagged by threshold but not assigned a class
    # (CR / glitch pair / end glitch). Tracked so the final ``rate``
    # excludes them — leaving them in biases the mean toward the
    # outlier amplitude even though we can't physically classify them.
    unclassFlagAccum = np.zeros((H, W, nDeltas), dtype=bool)

    # BAD-pixel pass (count criterion): tally 3σ delta excursions on the
    # PRISTINE input deltas. Must run before iter 1 because
    # ``_detectAndRepairOnce`` mutates ``deltas`` in place when
    # ``repair=True``. A single ``_rampQR`` pass gives the per-pixel
    # median + IQR-σ used to define an outlier; the rate criterion is
    # applied at the end against ``rateFinal``.
    if badPixelMinOutliers > 0:
        p25Init, p50Init, p75Init = _rampQR(deltas)
        sigmaInit = np.maximum(
            0.741 * (p75Init - p25Init).astype(np.float32, copy=False),
            sigmaFloorADU,
        )
        threshInit = (badPixelOutlierSigma * sigmaInit).astype(
            np.float32, copy=False
        )
        nLargeOutliers = np.zeros((H, W), dtype=np.int32)
        absSlice = np.empty((H, W), dtype=deltas.dtype)
        for k in range(nDeltas):
            np.abs(deltas[..., k] - p50Init, out=absSlice)
            nLargeOutliers += absSlice > threshInit
        del absSlice, threshInit, sigmaInit, p25Init, p50Init, p75Init
    else:
        nLargeOutliers = None

    nByIter = []
    iterTimings = []
    iteration = 0
    rateTolerance = 0.05  # ADU/read; loop stops when rate change < this

    # ---- Iter 1: process the full delta cube. ----
    tIter0 = time.time()
    rateFull, sigmaFull, newCR, newGlitch = _detectAndRepairOnce(
        deltas, goodPixelMask, glitchActive2D,
        crFlagAccum, glitchFlagAccum, boundaryFlagAccum, unclassFlagAccum,
        sigmaFloorADU, nSigma, repair, correctGlitches,
        glitchAmplitudeMinADU=glitchAmplitudeMinADU,
        maxDropFraction=maxDropFraction,
        nDropSigma=nDropSigma,
    )
    iterTimings.append(time.time() - tIter0)
    nByIter.append((newCR, newGlitch))

    # ---- Iter 2..N: restrict to the subset of pixels that got any flag in
    #      iter 1. Pixels with no iter-1 flag never see new flags in later
    #      iterations because their deltas are unchanged. ----
    activeMask = (crFlagAccum.any(axis=-1) | glitchFlagAccum.any(axis=-1)
                  | boundaryFlagAccum.any(axis=-1)
                  | unclassFlagAccum.any(axis=-1))
    ratePrev = rateFull

    for iteration in range(1, maxIterations):
        ys, xs = np.where(activeMask)
        if ys.size == 0:
            break

        tIter0 = time.time()
        # Subset views are copies (fancy indexing); modify, then write back.
        # Subset shape is (1, nActive, N-1) — a single-row column-vector
        # H'=1 keeps the (H', W', Nd) convention _detectAndRepairOnce
        # expects.
        deltasSub = deltas[ys, xs, :][np.newaxis, :, :].copy()
        crAccumSub = crFlagAccum[ys, xs, :][np.newaxis, :, :].copy()
        glitchAccumSub = glitchFlagAccum[ys, xs, :][np.newaxis, :, :].copy()
        boundaryAccumSub = boundaryFlagAccum[ys, xs, :][np.newaxis, :, :].copy()
        unclassAccumSub = unclassFlagAccum[ys, xs, :][np.newaxis, :, :].copy()
        goodSub = goodPixelMask[ys, xs][np.newaxis, :]
        glitchSub = glitchActive2D[ys, xs][np.newaxis, :]

        rateSub, sigmaSub, newCR, newGlitch = _detectAndRepairOnce(
            deltasSub, goodSub, glitchSub,
            crAccumSub, glitchAccumSub, boundaryAccumSub, unclassAccumSub,
            sigmaFloorADU, nSigma, repair, correctGlitches,
            glitchAmplitudeMinADU=glitchAmplitudeMinADU,
            maxDropFraction=maxDropFraction,
            nDropSigma=nDropSigma,
        )

        deltas[ys, xs, :] = deltasSub[0, :, :]
        crFlagAccum[ys, xs, :] = crAccumSub[0, :, :]
        glitchFlagAccum[ys, xs, :] = glitchAccumSub[0, :, :]
        boundaryFlagAccum[ys, xs, :] = boundaryAccumSub[0, :, :]
        unclassFlagAccum[ys, xs, :] = unclassAccumSub[0, :, :]

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

    # result.rate is the per-pixel UTR-weighted rate — the optimal
    # least-squares slope on the (CR-/glitch-corrected) ramp. Equivalent
    # to ``PfsIsrTask.calcUTRrates`` applied to ``read0 + cumsum(deltas)``,
    # derived directly in delta space via the closed-form delta weights
    # ``u[j] = 6(j+1)(N-1-j) / (N(N-1)(N+1))``, j = 0..N-2 — these sum
    # to 1 and reproduce the read-space UTR weights w[i] = (12i - 6(N-1))
    # / (N(N²-1)) after the cumsum/diff change of variables.
    #
    # Flagged-delta handling: CR + interior glitch pairs + end glitches
    # + unclassified outliers are all excluded from the rate.
    # Unclassified deltas are above-threshold residuals the classifier
    # could not assign (mostly negative outliers from RTS/persistence/
    # defects + CR candidates demoted by the cumulative-drop check) —
    # they are real ramp aberrations, not signal, and leaving them in
    # pulls the slope toward the outlier amplitude. Glitch pairs are
    # excluded regardless of ``correctGlitches`` because the UTR
    # weights are asymmetric across the pair — the +A/-A cancellation
    # that held for the unweighted mean does not extend here.
    # ``correctGlitches`` now controls only the in-cube repair
    # (whether the spikes survive in the returned ``deltas``), not
    # the rate.
    #
    # Exclusion is implemented as weight re-normalization: rate = sum
    # over unflagged k of u[k]·δ[k], divided by 1 − sum over flagged k
    # of u[k]. Equivalently, fill flagged deltas with the rate itself
    # and apply the closed-form UTR weights — i.e., the fixed point of
    # ``calcUTRrates(read0 + cumsum(deltas_with_flagged_filled_by_rate))``.
    # This makes ``result.rate`` agree exactly with calcUTRrates on the
    # reconstructed-and-repaired ramp; the median ``rateFull`` is kept
    # only as the all-flagged fallback. Accumulated slice-wise to avoid
    # an ``(H, W, N-1)`` transient.
    nReads = nDeltas + 1
    ks = np.arange(nDeltas, dtype=np.float32)
    utrW = np.float32(6.0) * (ks + np.float32(1.0)) * (
        np.float32(nReads - 1) - ks
    ) / np.float32(nReads * (nReads - 1) * (nReads + 1))
    sumUnflaggedW = np.zeros(deltas.shape[:-1], dtype=np.float32)
    sumFlaggedW = np.zeros(deltas.shape[:-1], dtype=np.float32)
    for k in range(nDeltas):
        flagK = (crFlagAccum[..., k] | boundaryFlagAccum[..., k]
                 | glitchFlagAccum[..., k] | unclassFlagAccum[..., k])
        unflagK = ~flagK
        sumUnflaggedW += utrW[k] * deltas[..., k] * unflagK
        sumFlaggedW += utrW[k] * flagK
    denom = np.float32(1.0) - sumFlaggedW
    with np.errstate(divide="ignore", invalid="ignore"):
        rateFinal = (sumUnflaggedW / denom).astype(np.float32, copy=False)
    allFlagged = denom < np.float32(1e-6)
    if allFlagged.any():
        rateFinal[allFlagged] = rateFull[allFlagged]
    sigmaFinal = sigmaFull.astype(np.float32, copy=False)

    # BAD-pixel pass: a pixel with ≥ badPixelMinOutliers delta
    # excursions above ``badPixelOutlierSigma × σ_IQR`` (counted on the
    # pristine input deltas before in-place repair) is RTS /
    # telegraph-noise. Restricted to ``goodPixelMask`` so already-masked
    # pixels don't appear in this output independently.
    if nLargeOutliers is not None:
        badPixelMask = (
            (nLargeOutliers >= badPixelMinOutliers) & goodPixelMask
        )
    else:
        badPixelMask = np.zeros((H, W), dtype=bool)

    # glitchFlagMask reports every ASIC glitch — interior pairs and end
    # glitches — for the ASIC_GLITCH mask plane. nGlitchPairs counts only
    # true pairs (from glitchFlagAccum); end glitches are singletons.
    return IterativeRepairResult(
        nIterations=iteration + 1,
        crFlagMask=crFlagAccum,
        glitchFlagMask=glitchFlagAccum | boundaryFlagAccum,
        unclassifiedFlagMask=unclassFlagAccum,
        badPixelMask=badPixelMask,
        rate=rateFinal,
        sigma=sigmaFinal,
        nCRs=int(crFlagAccum.sum()),
        nGlitchPairs=int(glitchFlagAccum.sum()) // 2,
        nByIteration=nByIter,
        iterationTimings=iterTimings,
    )
