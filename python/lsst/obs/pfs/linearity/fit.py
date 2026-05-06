"""Top-level fit(): tile-iterate over (H, W) and delegate to model.fitBlock."""

from __future__ import annotations

import os
from collections.abc import Sequence
from concurrent.futures import Future, ThreadPoolExecutor, as_completed

import numpy as np

from .models import Model, PolynomialModel
from .types import (
    BORDER_PIX,
    FIT_FAILED,
    INSUFFICIENT_POINTS,
    MASKED_BY_INPUT,
    NON_MONOTONIC,
    Diagnostics,
    LinearityCorrection,
    Ramp,
)

# Worker-count resolution constants. Tunable at module level; the tests
# monkeypatch `os.cpu_count` rather than these, so changing them does not
# break tests but will change the default behavior for small/large frames.
_SMALL_FRAME_PIXEL_LIMIT = 1_000_000   # H*W below this → sequential default
_DEFAULT_WORKER_CAP = 8                # auto-detected cpu_count is capped here

# Override point for tests. Default is the real ThreadPoolExecutor; a test
# can `monkeypatch.setattr("nirLinearity.fit._executorFactory", ...)` to observe
# construction or to inject a recording executor.
_executorFactory = ThreadPoolExecutor


def _resolveWorkerCount(workers: int | None, H: int, W: int) -> int:
    """Resolve the effective worker count for a `fit()` call.

    - If ``workers`` is an ``int``: returned as-is; must be >= 1.
    - If ``workers`` is ``None``:
        - H*W < ``_SMALL_FRAME_PIXEL_LIMIT`` → 1 (sequential default).
        - Otherwise → ``min(os.cpu_count() or 1, _DEFAULT_WORKER_CAP)``.

    Raises:
        ValueError: if ``workers`` is an int less than 1.
    """
    if workers is None:
        if H * W < _SMALL_FRAME_PIXEL_LIMIT:
            return 1
        return min(os.cpu_count() or 1, _DEFAULT_WORKER_CAP)
    if workers < 1:
        raise ValueError(f"workers must be >= 1, got {workers}")
    return workers


def fit(
    ramps: Sequence[Ramp],
    model: Model | None = None,
    blockSize: tuple[int, int] = (512, 512),
    workers: int | None = None,
    conditionNumberLimit: float = 1e12,
    deviationLimit: float | None = None,
    deviationStart: float = 0.5,
    nRefReads: int = 5,
    saturationLevel: float | None = None,
    lowFluxFraction: float = 0.5,
    borderWidth: int = 4,
) -> LinearityCorrection:
    """Fit a per-pixel nonlinearity correction from one or more ramps.

    See ``docs/superpowers/specs/2026-04-16-relin-package-design.md`` for
    the full algorithm description.

    Parameters
    ----------
    ramps : sequence of Ramp
        One or more ramps to fit jointly. All ramps must share the same
        ``(H, W)`` frame shape.
    model : Model, optional
        Model to fit. Defaults to ``PolynomialModel(order=4)``.
    blockSize : (int, int), optional
        Tile size in pixels for the per-tile normal-equations fit.
        Default is ``(512, 512)``. Smaller tiles reduce peak memory;
        larger tiles reduce per-tile overhead.
    workers : int or None, optional
        Number of worker threads for the tile loop.

        - ``1`` (explicit): sequential — no thread pool is constructed.
        - ``N > 1`` (explicit): run the tile loop on a
          ``ThreadPoolExecutor`` with ``max_workers=N``. No upper cap is
          applied to explicit values.
        - ``None`` (default): heuristic. If ``H * W < 1_000_000``, use
          ``1`` worker (sequential). Otherwise use
          ``min(os.cpu_count() or 1, 8)``. The cap at 8 applies only
          to the auto-detected default.

        Output is deterministic and worker-count-independent: every tile
        writes to a disjoint slice of the preallocated output arrays
        on the main thread, so the fit result is byte-identical
        regardless of ``workers``.

        Note on BLAS/LAPACK: numpy's linear-algebra routines may spawn
        additional internal threads. On multi-core machines, combining
        ``workers > 1`` with an uncontrolled BLAS thread count can lead
        to oversubscription and diminishing returns. If threaded
        speedup plateaus below expectations, try setting
        ``OMP_NUM_THREADS=1`` and/or ``MKL_NUM_THREADS=1`` in the
        process environment.
    conditionNumberLimit : float, optional
        Pixels whose normal-equations matrix has a condition number
        above this threshold are flagged as ``FIT_FAILED`` and left
        with zeroed coefficients. Default ``1e12``.
    deviationLimit : float or None, optional
        When set, reads where the per-read delta deviates from the
        per-pixel reference rate by more than this fraction are excluded
        from fitting. The reference rate is the median of the first
        ``nRefReads`` deltas for each pixel. For each pixel, the first
        read exceeding the threshold causes all subsequent reads to be
        masked as well. This clips the fitting range to the linear
        regime, excluding e.g. the saturated tail. Default ``None``
        (disabled — all valid reads are used).
    nRefReads : int, optional
        Number of early reads to median when computing the per-pixel
        reference rate for ``deviationLimit``. Default ``5``.
    saturationLevel : float or None, optional
        When set, reads where the cumulative signal ``m`` exceeds this
        value are excluded from fitting. For each pixel, the first read
        exceeding the threshold causes all subsequent reads to be masked.
        Default ``None`` (disabled).
    lowFluxFraction : float, optional
        Pixels whose median delta over the first ``nRefReads`` reads is
        below ``lowFluxFraction * globalMedianRate`` are masked entirely
        (all reads invalidated). This rejects dead or very dim pixels
        that would otherwise pass per-pixel deviation clipping.
        Default ``0.5``.
    borderWidth : int, optional
        Number of border pixels to exclude on each edge of the frame.
        These are hardware reference pixels and are flagged with
        ``BORDER_PIX``. Default ``4``.

    Returns
    -------
    LinearityCorrection
        Fitted coefficients, range bounds, bad-pixel mask, and
        diagnostics.
    """
    if model is None:
        model = PolynomialModel(order=4)

    if len(ramps) == 0:
        raise ValueError("fit() requires at least one ramp")

    # Validate shapes.
    H, W = ramps[0].reads.shape[1:]
    for k, ramp in enumerate(ramps):
        if ramp.reads.ndim != 3:
            raise ValueError(
                f"ramps[{k}].reads must be 3-D (N, H, W); got {ramp.reads.shape}"
            )
        if ramp.reads.shape[1:] != (H, W):
            raise ValueError(
                f"ramps[{k}].reads H,W = {ramp.reads.shape[1:]} "
                f"does not match ramps[0] H,W = {(H, W)}"
            )
        if ramp.validMask is not None and ramp.validMask.shape != (H, W):
            raise ValueError(
                f"ramps[{k}].validMask shape {ramp.validMask.shape} != {(H, W)}"
            )

    effectiveWorkers = _resolveWorkerCount(workers, H, W)

    # Per-ramp precomputation.
    cumulatives: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    for ramp in ramps:
        m = ramp.reads.astype(np.float32)
        cumulatives.append(m)
        # Rate R_k: median of first read (= first delta) over allowed pixels.
        firstRead = m[0]
        if ramp.validMask is not None:
            allowed = ramp.validMask == 0
            if allowed.any():
                rate = float(np.median(firstRead[allowed]))
            else:
                rate = float(np.median(firstRead))
        else:
            rate = float(np.median(firstRead))
        Nk = ramp.reads.shape[0]
        targets.append(rate * np.arange(1, Nk + 1, dtype=np.float32))

    # Deviation-based clipping: for each ramp, mask reads where the measured
    # cumulative deviates from the target by more than deviationLimit.
    # Store per-ramp (Nk, H, W) validity arrays for use in tile assembly.
    rampValidity: list[np.ndarray] = []
    for k, ramp in enumerate(ramps):
        Nk = ramp.reads.shape[0]
        if ramp.validMask is not None:
            v = np.broadcast_to(
                (ramp.validMask == 0)[None], (Nk, H, W)
            ).copy()
        else:
            v = np.ones((Nk, H, W), dtype=bool)

        # Mask border pixels (skip if frame is too small).
        if borderWidth > 0 and H > 2 * borderWidth and W > 2 * borderWidth:
            v[:, :borderWidth, :] = False
            v[:, -borderWidth:, :] = False
            v[:, :, :borderWidth] = False
            v[:, :, -borderWidth:] = False

        # Recover per-read deltas for deviation and low-flux checks.
        m = cumulatives[k]
        deltas = np.empty_like(m)
        deltas[0] = m[0]
        deltas[1:] = np.diff(m, axis=0)
        nRef = min(nRefReads, Nk)
        refDelta = np.median(deltas[:nRef], axis=0)  # (H, W)

        # Reject low-flux pixels: mask all reads for pixels whose
        # reference delta is below lowFluxFraction * global median rate.
        rate = float(targets[k][0])  # targets[k][0] = rate * 1
        lowFluxThreshold = lowFluxFraction * rate
        lowFlux = refDelta < lowFluxThreshold  # (H, W)
        v[:, lowFlux] = False

        if deviationLimit is not None:
            startRead = int(Nk * deviationStart)
            with np.errstate(divide="ignore", invalid="ignore"):
                frac = np.abs(deltas - refDelta[None]) / np.abs(refDelta[None])
            frac = np.where(np.isfinite(frac), frac, 0.0)
            exceeds = frac > deviationLimit  # (Nk, H, W)
            # Only apply from startRead onward (early reads are not
            # near saturation and shouldn't trigger the limit).
            exceeds[:startRead] = False
            exceeds = np.maximum.accumulate(exceeds, axis=0)
            v[exceeds] = False

        if saturationLevel is not None:
            m = cumulatives[k]  # (Nk, H, W)
            saturated = m > saturationLevel
            saturated = np.maximum.accumulate(saturated, axis=0)
            v[saturated] = False

        rampValidity.append(v)

    # Concatenated targets across ramps — used per tile.
    tConcat = np.concatenate(targets)

    # Preallocate full-frame outputs.
    # The shape of coefficients is determined by the model:
    # - PolynomialModel: (order+1, H, W)
    # We discover the coefficient shape by running a 1x1 dummy block first.
    coefShape = _peekCoefShape(model)
    coefficients = np.zeros((coefShape, H, W), dtype=np.float32)
    fitMin = np.zeros((H, W), dtype=np.float32)
    fitMax = np.zeros((H, W), dtype=np.float32)
    residualRms = np.zeros((H, W), dtype=np.float32)
    maxAbsResidual = np.zeros((H, W), dtype=np.float32)
    nPointsUsed = np.zeros((H, W), dtype=np.int32)
    conditionNumber = np.zeros((H, W), dtype=np.float32)
    monotonic = np.zeros((H, W), dtype=bool)
    badPixelMask = np.zeros((H, W), dtype=np.uint8)

    # Iterate over tiles. Tile-assembly (mTile, validTile) is identical
    # for sequential and threaded paths; factor it into a closure so both
    # paths call model.fitBlock with exactly the same inputs.
    bH, bW = blockSize

    def _assembleTile(
        rowStart: int, rowEnd: int, colStart: int, colEnd: int
    ) -> tuple[np.ndarray, np.ndarray]:
        mSegments: list[np.ndarray] = []
        validSegments: list[np.ndarray] = []
        for k in range(len(ramps)):
            mSegments.append(
                cumulatives[k][:, rowStart:rowEnd, colStart:colEnd]
            )
            validSegments.append(
                rampValidity[k][:, rowStart:rowEnd, colStart:colEnd]
            )
        mTile = np.concatenate(mSegments, axis=0)
        validTile = np.concatenate(validSegments, axis=0)
        return mTile, validTile

    def _storeResult(
        rowStart: int, rowEnd: int, colStart: int, colEnd: int, result
    ) -> None:
        coefficients[:, rowStart:rowEnd, colStart:colEnd] = result.coefficients
        fitMin[rowStart:rowEnd, colStart:colEnd] = result.fitMin
        fitMax[rowStart:rowEnd, colStart:colEnd] = result.fitMax
        residualRms[rowStart:rowEnd, colStart:colEnd] = result.residualRms
        maxAbsResidual[rowStart:rowEnd, colStart:colEnd] = result.maxAbsResidual
        nPointsUsed[rowStart:rowEnd, colStart:colEnd] = result.nPointsUsed
        conditionNumber[rowStart:rowEnd, colStart:colEnd] = result.conditionNumber
        monotonic[rowStart:rowEnd, colStart:colEnd] = result.monotonic
        badPixelMask[rowStart:rowEnd, colStart:colEnd] = result.badPixelMask

    if effectiveWorkers == 1:
        # Sequential fast path — no executor involvement.
        for rowStart in range(0, H, bH):
            rowEnd = min(rowStart + bH, H)
            for colStart in range(0, W, bW):
                colEnd = min(colStart + bW, W)
                mTile, validTile = _assembleTile(
                    rowStart, rowEnd, colStart, colEnd
                )
                result = model.fitBlock(
                    m=mTile, t=tConcat, valid=validTile,
                    conditionNumberLimit=conditionNumberLimit,
                )
                _storeResult(rowStart, rowEnd, colStart, colEnd, result)
    else:
        # Threaded path. Submit each tile as a future; consume completed
        # futures on the main thread and stitch into disjoint output slices.
        # Tile-assembly runs on the submitting thread so workers do pure
        # compute on independent numpy arrays (no shared mutable state).
        # Note: ThreadPoolExecutor.submit has no back-pressure, so all
        # tiles' assembled (mTile, validTile) arrays coexist in memory
        # until their futures complete. On 4096x4096 with blockSize=(512,
        # 512), that is 64 tiles of roughly tile-sized float32 plus a
        # small bool array each — manageable on the reference workload
        # but worth keeping in mind if tile size grows.
        with _executorFactory(max_workers=effectiveWorkers) as executor:
            futures: dict[Future, tuple[int, int, int, int]] = {}
            for rowStart in range(0, H, bH):
                rowEnd = min(rowStart + bH, H)
                for colStart in range(0, W, bW):
                    colEnd = min(colStart + bW, W)
                    mTile, validTile = _assembleTile(
                        rowStart, rowEnd, colStart, colEnd
                    )
                    fut = executor.submit(
                        model.fitBlock,
                        m=mTile, t=tConcat, valid=validTile,
                        conditionNumberLimit=conditionNumberLimit,
                    )
                    futures[fut] = (rowStart, rowEnd, colStart, colEnd)

            for fut in as_completed(futures):
                rs, re, cs, ce = futures[fut]
                try:
                    result = fut.result()
                except Exception as e:
                    # Cancel any futures that haven't started; in-flight
                    # tasks still run to completion but their results are
                    # discarded when the `with` block shuts down.
                    for other in futures:
                        if other is not fut:
                            other.cancel()
                    raise RuntimeError(
                        f"fitBlock failed at tile "
                        f"[rows {rs}:{re}, cols {cs}:{ce}]"
                    ) from e
                _storeResult(rs, re, cs, ce, result)

    # Flag border pixels (skip if frame is too small).
    if borderWidth > 0 and H > 2 * borderWidth and W > 2 * borderWidth:
        badPixelMask[:borderWidth, :] |= BORDER_PIX
        badPixelMask[-borderWidth:, :] |= BORDER_PIX
        badPixelMask[:, :borderWidth] |= BORDER_PIX
        badPixelMask[:, -borderWidth:] |= BORDER_PIX

    # Propagate caller-supplied input masks.
    for ramp in ramps:
        if ramp.validMask is not None:
            inputBad = (ramp.validMask != 0)
            badPixelMask[inputBad] |= MASKED_BY_INPUT

    # Dataset-wide summary.
    goodPixels = (badPixelMask == 0)
    totalPixels = int(H * W)
    summary: dict = {
        "totalPixels": totalPixels,
        "goodPixelFraction": float(goodPixels.sum()) / totalPixels,
        "badPixelFraction_borderPix": float((badPixelMask & BORDER_PIX > 0).sum()) / totalPixels,
        "badPixelFraction_maskedByInput": float((badPixelMask & MASKED_BY_INPUT > 0).sum()) / totalPixels,
        "badPixelFraction_insufficientPoints": float((badPixelMask & INSUFFICIENT_POINTS > 0).sum()) / totalPixels,
        "badPixelFraction_fitFailed": float((badPixelMask & FIT_FAILED > 0).sum()) / totalPixels,
        "badPixelFraction_nonMonotonic": float((badPixelMask & NON_MONOTONIC > 0).sum()) / totalPixels,
        "modelName": model.modelName,
        "nRamps": len(ramps),
        "order": model.order,
        "blockSize": f"{blockSize[0]}x{blockSize[1]}",
        "borderWidth": borderWidth,
        "conditionNumberLimit": conditionNumberLimit,
        "nRefReads": nRefReads,
        "lowFluxFraction": lowFluxFraction,
    }
    if deviationLimit is not None:
        summary["deviationLimit"] = deviationLimit
        summary["deviationStart"] = deviationStart
    if saturationLevel is not None:
        summary["saturationLevel"] = saturationLevel
    if goodPixels.any():
        goodRms = residualRms[goodPixels]
        summary["residualRmsP50"] = float(np.percentile(goodRms, 50))
        summary["residualRmsP95"] = float(np.percentile(goodRms, 95))
        summary["residualRmsP99"] = float(np.percentile(goodRms, 99))
    else:
        summary["residualRmsP50"] = float("nan")
        summary["residualRmsP95"] = float("nan")
        summary["residualRmsP99"] = float("nan")

    diagnostics = Diagnostics(
        residualRms=residualRms,
        maxAbsResidual=maxAbsResidual,
        nPointsUsed=nPointsUsed,
        monotonic=monotonic,
        conditionNumber=conditionNumber,
        summary=summary,
    )

    return LinearityCorrection(
        model=model,
        coefficients=coefficients,
        fitMin=fitMin,
        fitMax=fitMax,
        badPixelMask=badPixelMask,
        diagnostics=diagnostics,
    )


def _peekCoefShape(model: Model) -> int:
    """Return the first-axis size of coefficients the model will produce.

    For ``PolynomialModel``, this is ``order + 1``.
    Models that don't expose ``order`` must run a throwaway 1x1 fit.
    """
    if isinstance(model, PolynomialModel):
        return model.order + 1
    # Fallback: run a minimal 1x1 block fit with 2*(order+1) dummy points.
    nPoints = 8
    m = np.linspace(0.0, 1.0, nPoints, dtype=np.float32)[:, None, None]
    t = m[:, 0, 0].copy()
    valid = np.ones((nPoints, 1, 1), dtype=bool)
    result = model.fitBlock(
        m=m, t=t, valid=valid, conditionNumberLimit=1e12
    )
    return int(result.coefficients.shape[0])
