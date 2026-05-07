"""Validate H4 linearity corrections by comparing first/second-half rates.

Run a NIR ramp through the full ISR (dark + linearization + defects), split
the resulting cumulative-flux cube into halves, fit a per-pixel rate to each
half, and compare. A correctly-tuned linearity model flattens the half-vs-
half rate; residual structure points at a poor or wrong correction.

Run twice (linearize on / off) so the *improvement* is visible, not just the
absolute residual.

Usage from JupyterLab::

    from lsst.obs.pfs.h4Linearity import validate
    cmp = validate.runComparison(butler, dataId)
    validate.printComparison(cmp)
    validate.plotComparison(cmp)

Notes
-----
- ``PfsIsrTask.makeNirExposure`` mutates the underlying ramp arrays in
  place (in-place ``cumsum`` and dark subtraction). Each call therefore
  needs a fresh ``raw``; ``runComparison`` re-fetches between passes.
- The returned cube is *pre-gain* ADU. ``minRate`` and ``fluxBins``
  defaults are in pre-gain ADU.
- Linearity coefficients exist only for the ``n3`` detector today. Calling
  with ``doLinearize=True`` against any other detector raises before any
  ISR work is done.
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from lsst.obs.pfs import isrTask as pfsIsrTask


SUPPORTED_LINEARITY_CAMS = ("n3",)

DEFAULT_MIN_RATE = 3.0
# "auto" → bins are derived per-call from minRate and the ramp midpoint:
#   [(0, minRate*midRead), (minRate*midRead, 100), (100, 1000),
#    (1000, 5000), (5000, 10000), (10000, 50000)]
# Pass an explicit sequence of (lo, hi) pairs to override.
DEFAULT_FLUX_BINS = "auto"
DEFAULT_AGREE_PCT = (0.5, 1.0, 2.0)


def _resolveFluxBins(fluxBins, minRate, midRead):
    if fluxBins == "auto":
        thr = float(minRate) * float(midRead)
        return (
            (0.0, thr),
            (thr, 100.0),
            (100.0, 1_000.0),
            (1_000.0, 5_000.0),
            (5_000.0, 10_000.0),
            (10_000.0, 50_000.0),
        )
    return tuple((float(lo), float(hi)) for lo, hi in fluxBins)


@dataclass
class HalvesResult:
    slope1: np.ndarray
    slope2: np.ndarray
    meanFlux: np.ndarray
    relDiff: np.ndarray
    goodMask: np.ndarray         # unmasked pixels with finite slopes
    statsMask: np.ndarray        # goodMask AND avgRate > minRate
    exposureMaskBits: dict
    nReads: int
    midRead: int
    cam: str
    doLinearize: bool
    minRate: float
    summary: dict
    cube: Optional[np.ndarray] = None


@dataclass
class ComparisonResult:
    on: HalvesResult
    off: HalvesResult
    delta: dict = field(default_factory=dict)


def _makeDefaultIsrTask(doLinearize: bool):
    """Build a `PfsIsrTask` matching the notebook recipe for this test."""
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doFlat = False
    config.doDark = True
    config.doDefect = True
    config.doSaturationInterpolation = False

    config.h4.quickCDS = False
    config.h4.doCR = False
    config.h4.doIPC = False
    config.h4.doWriteRawCube = True
    config.h4.doLinearize = doLinearize
    # Rate-based CR rejection requires linearization to have run.
    config.h4.doRateCR = doLinearize

    config.validate()
    return pfsIsrTask.PfsIsrTask(config=config)


def _ensureIsrTask(isrTask, doLinearize: bool):
    """Return an isrTask whose configuration matches the requested mode."""
    if isrTask is None:
        return _makeDefaultIsrTask(doLinearize)

    needsRebuild = (
        bool(isrTask.config.h4.doLinearize) != bool(doLinearize)
        or bool(isrTask.config.h4.doRateCR) != bool(doLinearize)
        or not bool(isrTask.config.h4.doWriteRawCube)
        or bool(isrTask.config.h4.quickCDS)
    )
    if needsRebuild:
        return _makeDefaultIsrTask(doLinearize)
    return isrTask


def _resolveCam(cam, raw, butler, dataId):
    if cam is not None:
        return cam, raw
    if raw is None:
        raw = butler.get("raw", dataId)
    return raw.detector.getName(), raw


def _uniformSlope(cube: np.ndarray) -> np.ndarray:
    """Per-pixel least-squares slope for uniform x = arange(k).

    Parameters
    ----------
    cube : np.ndarray
        Shape ``(k, H, W)``. Treated as ``y[k]`` per pixel.

    Returns
    -------
    slope : np.ndarray
        Shape ``(H, W)``, dtype float32. Slope of ``cube`` against
        ``arange(k)``. Constant offsets in ``cube`` cancel.
    """
    k = cube.shape[0]
    if k < 2:
        raise ValueError(f"Need at least 2 reads for a slope; got {k}.")
    x = np.arange(k, dtype=np.float64)
    Sx = float(x.sum())
    Sxx = float((x * x).sum())
    denom = k * Sxx - Sx * Sx

    Sc = cube.sum(axis=0, dtype=np.float64)
    Sxc = np.einsum("k,khw->hw", x, cube, dtype=np.float64)

    slope = (k * Sxc - Sx * Sc) / denom
    return slope.astype(np.float32)


def _mad(values: np.ndarray) -> float:
    if values.size == 0:
        return float("nan")
    med = float(np.median(values))
    return float(1.4826 * np.median(np.abs(values - med)))


def _binStats(values: np.ndarray, agreePct):
    if values.size == 0:
        return {
            "median": float("nan"),
            "mad": float("nan"),
            "p1": float("nan"),
            "p16": float("nan"),
            "p84": float("nan"),
            "p99": float("nan"),
            "fracWithin": {pct: float("nan") for pct in agreePct},
            "nPix": 0,
        }
    p1, p16, med, p84, p99 = np.percentile(values, [1, 16, 50, 84, 99])
    return {
        "median": float(med),
        "mad": _mad(values),
        "p1": float(p1),
        "p16": float(p16),
        "p84": float(p84),
        "p99": float(p99),
        "fracWithin": {
            pct: float(np.mean(np.abs(values) < pct / 100.0)) for pct in agreePct
        },
        "nPix": int(values.size),
    }


def _summarize(relDiff, meanFlux, goodMask, statsMask, fluxBins, agreePct):
    nTotal = int(goodMask.size)
    nGood = int(goodMask.sum())
    nStats = int(statsMask.sum())

    if nStats == 0:
        empty = _binStats(np.array([], dtype=np.float32), agreePct)
        return {
            "nGood": nGood,
            "nStats": 0,
            "nTotal": nTotal,
            "overall": empty,
            "byBin": {(lo, hi): empty for (lo, hi) in fluxBins},
        }

    flat = relDiff[statsMask]
    flatMean = meanFlux[statsMask]

    byBin = {}
    for lo, hi in fluxBins:
        sel = (flatMean >= lo) & (flatMean < hi)
        byBin[(lo, hi)] = _binStats(flat[sel], agreePct)

    overall = _binStats(flat, agreePct)
    overall.pop("nPix", None)

    return {
        "nGood": nGood,
        "nStats": nStats,
        "nTotal": nTotal,
        "overall": overall,
        "byBin": byBin,
    }


def _exposureMaskBitCounts(exp):
    bits = {}
    arr = exp.mask.array
    mask = exp.mask
    for name in ("BAD", "SAT", "INTRP", "CR", "EDGE", "DETECTED", "NO_DATA"):
        try:
            bit = mask.getPlaneBitMask(name)
        except Exception:
            continue
        bits[name] = int(((arr & bit) != 0).sum())
    return bits


def runHalvesTest(
    butler,
    dataId,
    *,
    doLinearize: bool = True,
    cam: Optional[str] = None,
    isrTask=None,
    raw=None,
    nirDark=None,
    defects=None,
    minRate: float = DEFAULT_MIN_RATE,
    fluxBins=DEFAULT_FLUX_BINS,
    agreePct=DEFAULT_AGREE_PCT,
    returnCube: bool = False,
    log=None,
) -> HalvesResult:
    """Run NIR ISR once and report half-vs-half rate statistics.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
    dataId : mapping
    doLinearize : bool
        Whether to apply the new H4 per-read linearity correction.
    cam : str, optional
        Detector name (e.g. ``"n3"``). If ``None``, inferred from the raw.
        Required to be ``"n3"`` when ``doLinearize=True``.
    isrTask : `PfsIsrTask`, optional
        Reuse a prebuilt task. If its config disagrees with the requested
        mode, a fresh task is built (the supplied task is not mutated).
    raw, nirDark, defects : optional
        Pre-fetched butler objects. Fetched on demand if not given.
    minRate : float
        Pixels with mean rate < ``minRate`` ADU/read are excluded from
        statistics (avoids dividing by ~zero in relDiff).
    fluxBins : sequence of (lo, hi)
        Cumulative-flux bins (in pre-gain ADU at the ramp midpoint) used
        for per-bin statistics.
    agreePct : sequence of float
        Percent thresholds for the ``fracWithin`` summary.
    returnCube : bool
        If True, the returned ``HalvesResult`` carries the full flux cube.
        Off by default to keep peak memory at one cube.
    log : logger, optional

    Returns
    -------
    HalvesResult
    """
    cam, raw = _resolveCam(cam, raw, butler, dataId)

    if doLinearize and cam not in SUPPORTED_LINEARITY_CAMS:
        raise RuntimeError(
            f"Linearity corrections are only available for {SUPPORTED_LINEARITY_CAMS}; "
            f"got cam={cam!r}. Re-run with doLinearize=False to test this detector."
        )

    isrTask = _ensureIsrTask(isrTask, doLinearize)
    if log is None:
        log = isrTask.log

    if raw is None:
        raw = butler.get("raw", dataId)
    if nirDark is None:
        nirDark = butler.get("nirDark", dataId)
    if defects is None:
        defects = butler.get("defects", dataId)

    linearity = isrTask.resolveNirLinearity(cam) if doLinearize else None
    if doLinearize and linearity is None:
        raise RuntimeError(
            f"resolveNirLinearity({cam!r}) returned None; cannot linearize."
        )

    log.info(
        f"validate.runHalvesTest: cam={cam} doLinearize={doLinearize} "
        f"visit={dataId.get('visit', '?')}"
    )

    exp, cube = isrTask.makeNirExposure(
        raw,
        nirDark=nirDark,
        defects=defects,
        linearity=linearity,
        doReturnRawCube=True,
    )
    if cube is None:
        raise RuntimeError(
            "makeNirExposure returned no flux cube; check h4.doWriteRawCube and h4.quickCDS."
        )

    # makeNirExposure() only flags defects in the mask plane when
    # linearization is enabled (via the MASKED_BY_INPUT path). PfsIsrTask.run()
    # normally calls maskDefect() afterwards; we are bypassing that, so apply
    # defects ourselves. Idempotent when defects are already flagged.
    if defects is not None:
        defects.maskPixels(exp.mask, "BAD")

    n = cube.shape[0]
    mid = n // 2
    if mid < 2 or n - mid < 2:
        raise RuntimeError(
            f"Ramp too short to split into halves with >= 2 reads each (n={n})."
        )

    slope1 = _uniformSlope(cube[:mid])
    slope2 = _uniformSlope(cube[mid:])
    meanFlux = cube.mean(axis=0, dtype=np.float64).astype(np.float32)

    unmasked = exp.mask.array == 0
    goodMask = unmasked & np.isfinite(slope1) & np.isfinite(slope2)

    with np.errstate(divide="ignore", invalid="ignore"):
        denom = slope1 + slope2
        relDiff = np.where(denom != 0, 2.0 * (slope1 - slope2) / denom, np.nan).astype(
            np.float32
        )

    avgSlope = 0.5 * (slope1 + slope2)
    statsMask = goodMask & (avgSlope > minRate)

    fluxBins = _resolveFluxBins(fluxBins, minRate, mid)
    summary = _summarize(relDiff, meanFlux, goodMask, statsMask, fluxBins, agreePct)
    exposureMaskBits = _exposureMaskBitCounts(exp)

    result = HalvesResult(
        slope1=slope1,
        slope2=slope2,
        meanFlux=meanFlux,
        relDiff=relDiff,
        goodMask=goodMask,
        statsMask=statsMask,
        exposureMaskBits=exposureMaskBits,
        nReads=n,
        midRead=mid,
        cam=cam,
        doLinearize=doLinearize,
        minRate=float(minRate),
        summary=summary,
        cube=cube if returnCube else None,
    )

    if not returnCube:
        del cube
    return result


def runComparison(
    butler,
    dataId,
    *,
    cam: Optional[str] = None,
    isrTask=None,
    nirDark=None,
    defects=None,
    minRate: float = DEFAULT_MIN_RATE,
    fluxBins=DEFAULT_FLUX_BINS,
    agreePct=DEFAULT_AGREE_PCT,
    log=None,
) -> ComparisonResult:
    """Run the halves test twice (linearize on, then off) and compare.

    Re-fetches ``raw`` between the two passes since ``makeNirExposure``
    mutates the ramp data in place. ``defects`` and ``nirDark`` are
    read-only consumers and are reused.

    Returns
    -------
    ComparisonResult
    """
    cam, _raw = _resolveCam(cam, None, butler, dataId)

    if nirDark is None:
        nirDark = butler.get("nirDark", dataId)
    if defects is None:
        defects = butler.get("defects", dataId)

    common = dict(
        cam=cam,
        isrTask=isrTask,
        nirDark=nirDark,
        defects=defects,
        minRate=minRate,
        fluxBins=fluxBins,
        agreePct=agreePct,
        returnCube=False,
        log=log,
    )

    onResult = runHalvesTest(
        butler,
        dataId,
        doLinearize=True,
        raw=butler.get("raw", dataId),
        **common,
    )
    offResult = runHalvesTest(
        butler,
        dataId,
        doLinearize=False,
        raw=butler.get("raw", dataId),
        **common,
    )

    medOn = onResult.summary["overall"]["median"]
    medOff = offResult.summary["overall"]["median"]
    madOn = onResult.summary["overall"]["mad"]
    madOff = offResult.summary["overall"]["mad"]
    delta = {
        "median_on": medOn,
        "median_off": medOff,
        "mad_on": madOn,
        "mad_off": madOff,
        "madImprovementFactor": (madOff / madOn) if madOn and np.isfinite(madOn) else float("nan"),
    }
    return ComparisonResult(on=onResult, off=offResult, delta=delta)


def _formatBinTable(byBin, agreePct):
    lines = []
    headers = ["fluxBin", "nPix", "median", "MAD"]
    headers += [f"|<{p}%" for p in agreePct]
    lines.append("  ".join(f"{h:>12s}" for h in headers))
    for (lo, hi), s in byBin.items():
        row = [f"[{lo},{hi})", f"{s['nPix']}", f"{s['median']:+.4f}", f"{s['mad']:.4f}"]
        row += [f"{s['fracWithin'][p]:.3f}" for p in agreePct]
        lines.append("  ".join(f"{c:>12s}" for c in row))
    return "\n".join(lines)


def summarize(result: HalvesResult) -> str:
    """Return a multi-line text summary of a single halves-test result."""
    s = result.summary
    o = s["overall"]
    agreePct = sorted(o["fracWithin"].keys())
    nTotal = max(s["nTotal"], 1)
    lines = [
        f"cam={result.cam}  doLinearize={result.doLinearize}  "
        f"nReads={result.nReads}  midRead={result.midRead}",
        f"good (unmasked) pixels: {s['nGood']}/{s['nTotal']} "
        f"({100.0 * s['nGood'] / nTotal:.2f}%)",
        f"in stats (avgRate > minRate): {s['nStats']}/{s['nTotal']} "
        f"({100.0 * s['nStats'] / nTotal:.2f}%)",
        f"mask bits: " + ", ".join(f"{k}={v}" for k, v in result.exposureMaskBits.items()),
        f"overall relDiff: median={o['median']:+.4f}  MAD={o['mad']:.4f}  "
        f"p16/p84=[{o['p16']:+.4f},{o['p84']:+.4f}]  "
        f"p1/p99=[{o['p1']:+.4f},{o['p99']:+.4f}]",
        "fraction within: "
        + "  ".join(f"|relDiff|<{p}%: {o['fracWithin'][p]:.3f}" for p in agreePct),
        "by flux bin:",
        _formatBinTable(s["byBin"], agreePct),
    ]
    return "\n".join(lines)


def printComparison(cmp: ComparisonResult) -> None:
    """Print both summaries plus an improvement-factor line."""
    print("=== linearize=OFF ===")
    print(summarize(cmp.off))
    print()
    print("=== linearize=ON ===")
    print(summarize(cmp.on))
    print()
    d = cmp.delta
    print(
        f"Δ median: off={d['median_off']:+.4f}  on={d['median_on']:+.4f}\n"
        f"Δ MAD:    off={d['mad_off']:.4f}   on={d['mad_on']:.4f}   "
        f"improvement (off/on) = {d['madImprovementFactor']:.2f}x"
    )


def plotHalves(
    result: HalvesResult,
    *,
    fig=None,
    sample: int = 200_000,
    vmin: float = -0.05,
    vmax: float = 0.05,
    title: Optional[str] = None,
):
    """Three-panel diagnostic figure for a single halves-test result."""
    import matplotlib.pyplot as plt

    ownFigure = fig is None
    if ownFigure:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    else:
        axes = fig.subplots(1, 3)

    relDiff = result.relDiff
    good = result.goodMask
    stats = result.statsMask
    flat = relDiff[stats]
    fluxFlat = result.meanFlux[stats]

    ax = axes[0]
    if flat.size:
        clipped = flat[np.abs(flat) <= max(abs(vmin), abs(vmax))]
        ax.hist(clipped, bins=120, color="0.3")
        med = float(np.median(flat))
        mad = _mad(flat)
        ax.axvline(med, color="C3", lw=1.5, label=f"median={med:+.4f}")
        ax.axvline(med - mad, color="C3", lw=1.0, ls="--", label=f"±MAD={mad:.4f}")
        ax.axvline(med + mad, color="C3", lw=1.0, ls="--")
        ax.legend(loc="upper right", fontsize=9)
    ax.set_xlabel("relDiff = 2(s1-s2)/(s1+s2)")
    ax.set_ylabel(f"pixels with avgRate > {result.minRate:g}")
    ax.set_title("histogram (clipped)")

    ax = axes[1]
    img = np.where(good, relDiff, np.nan)
    im = ax.imshow(img, origin="lower", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    ax.set_title("relDiff map (unmasked pixels)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax = axes[2]
    if flat.size:
        if flat.size > sample:
            idx = np.random.default_rng(0).choice(flat.size, size=sample, replace=False)
            xs = fluxFlat[idx]
            ys = flat[idx]
        else:
            xs = fluxFlat
            ys = flat
        ax.hexbin(xs, ys, gridsize=80, cmap="viridis", mincnt=1, bins="log")
        ax.set_ylim(vmin, vmax)
        for (lo, hi), s in result.summary["byBin"].items():
            xc = 0.5 * (lo + hi)
            ax.errorbar(xc, s["median"], yerr=s["mad"], fmt="o", color="C3", capsize=4)
            ax.axvline(lo, color="0.7", lw=0.5)
            ax.axvline(hi, color="0.7", lw=0.5)
    ax.set_xlabel("mean cumulative flux (ADU)")
    ax.set_ylabel("relDiff")
    ax.set_title("relDiff vs mean flux")

    suptitle = title or (
        f"{result.cam}  doLinearize={result.doLinearize}  "
        f"nReads={result.nReads}  midRead={result.midRead}"
    )
    fig.suptitle(suptitle)
    if ownFigure:
        fig.tight_layout()
    return fig


def plotStatsMask(result: HalvesResult, *, ax=None):
    """Show which pixels are in statsMask (used for relDiff statistics).

    Three-class map:
      0 = masked (goodMask=False)
      1 = unmasked but below minRate (goodMask=True, statsMask=False)
      2 = in stats (statsMask=True)
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))
    else:
        fig = ax.figure

    classes = np.zeros(result.goodMask.shape, dtype=np.uint8)
    classes[result.goodMask & ~result.statsMask] = 1
    classes[result.statsMask] = 2

    cmap = ListedColormap(["#222222", "#e8a838", "#3a8fdc"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    im = ax.imshow(classes, origin="lower", cmap=cmap, norm=norm,
                   interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2], fraction=0.046, pad=0.02)
    cbar.ax.set_yticklabels([
        f"masked ({(~result.goodMask).sum()})",
        f"unmasked, rate≤{result.minRate:g}"
        f" ({int((result.goodMask & ~result.statsMask).sum())})",
        f"in stats ({int(result.statsMask.sum())})",
    ])
    ax.set_title(
        f"statsMask  cam={result.cam}  doLinearize={result.doLinearize}  "
        f"minRate={result.minRate:g}"
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return fig


def outlierMask(result: HalvesResult, *, thresh: float = 0.5) -> np.ndarray:
    """Return the boolean 2-D mask of statsMask pixels with |relDiff| > thresh."""
    return result.statsMask & (np.abs(result.relDiff) > thresh)


def plotOutliers(
    result: HalvesResult,
    *,
    thresh: float = 0.5,
    ax=None,
):
    """2-D map of statsMask pixels with |relDiff| > thresh.

    Background: gray = masked, light = in stats and within threshold.
    Foreground: red = positive outlier (s1 > s2), blue = negative outlier.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))
    else:
        fig = ax.figure

    bad = outlierMask(result, thresh=thresh)
    posOut = bad & (result.relDiff > 0)
    negOut = bad & (result.relDiff < 0)

    # 0 = masked, 1 = in-stats good, 2 = neg outlier, 3 = pos outlier
    classes = np.zeros(result.goodMask.shape, dtype=np.uint8)
    classes[result.goodMask & ~result.statsMask] = 1
    classes[result.statsMask] = 1
    classes[negOut] = 2
    classes[posOut] = 3

    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(["#222222", "#dddddd", "#3a6fff", "#dc3a3a"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    im = ax.imshow(classes, origin="lower", cmap=cmap, norm=norm,
                   interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2, 3], fraction=0.046, pad=0.02)
    cbar.ax.set_yticklabels([
        f"masked",
        f"in-bounds",
        f"relDiff < -{thresh:g} ({int(negOut.sum())})",
        f"relDiff > +{thresh:g} ({int(posOut.sum())})",
    ])
    ax.set_title(
        f"|relDiff| > {thresh:g} outliers  cam={result.cam}  "
        f"doLinearize={result.doLinearize}  total={int(bad.sum())}"
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return fig


def plotOutlierTraces(
    result: HalvesResult,
    cube: np.ndarray,
    *,
    thresh: float = 0.5,
    n: int = 20,
    normalize: str = "mean",
    coords=None,
    ax=None,
    rng=None,
    alpha: float = 0.6,
):
    """Plot per-read flux traces for outlier pixels.

    Parameters
    ----------
    result : HalvesResult
    cube : np.ndarray
        The (nreads, H, W) flux cube. Pass the cube returned from
        ``runHalvesTest(..., returnCube=True)``.
    thresh : float
        |relDiff| threshold defining outliers.
    n : int
        Number of outlier pixels to draw (random sample). Ignored if
        ``coords`` is given.
    normalize : {"mean", "max", "none"}
        Per-pixel normalization of the trace.
    coords : sequence of (y, x), optional
        Specific pixel coordinates to plot. Overrides ``thresh``/``n``.
    ax : matplotlib axis, optional
    rng : np.random.Generator, optional
        For reproducible random subsets. Defaults to seed=0.
    alpha : float
        Line alpha.

    Returns
    -------
    fig : matplotlib.figure.Figure
    coords : np.ndarray
        Shape ``(npicked, 2)`` of (y, x) pairs that were plotted.
    """
    import matplotlib.pyplot as plt

    if cube is None:
        raise ValueError(
            "cube is None. Re-run runHalvesTest(..., returnCube=True) and pass "
            "the resulting result.cube here."
        )
    if cube.shape[1:] != result.relDiff.shape:
        raise ValueError(
            f"cube spatial shape {cube.shape[1:]} != relDiff shape {result.relDiff.shape}"
        )

    if coords is None:
        bad = outlierMask(result, thresh=thresh)
        ys, xs = np.where(bad)
        if ys.size == 0:
            raise RuntimeError(
                f"No statsMask pixels with |relDiff|>{thresh}. Lower the threshold "
                "or check that statsMask has the pixels you expect."
            )
        if rng is None:
            rng = np.random.default_rng(0)
        if ys.size > n:
            pick = rng.choice(ys.size, size=n, replace=False)
            ys, xs = ys[pick], xs[pick]
        coords = np.column_stack([ys, xs])
    else:
        coords = np.asarray(coords)
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("coords must be (N, 2) of (y, x) pairs.")
        ys = coords[:, 0]
        xs = coords[:, 1]

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.figure

    nReads = cube.shape[0]
    t = np.arange(nReads)
    traces = cube[:, ys, xs].astype(np.float64)  # shape (nReads, npicked)

    if normalize == "mean":
        denom = traces.mean(axis=0, keepdims=True)
        ylabel = "flux / mean(flux)"
    elif normalize == "max":
        denom = np.abs(traces).max(axis=0, keepdims=True)
        ylabel = "flux / max(|flux|)"
    elif normalize == "none":
        denom = 1.0
        ylabel = "cumulative flux (ADU)"
    else:
        raise ValueError(f"normalize must be 'mean', 'max', or 'none'; got {normalize!r}")
    denom = np.where(denom == 0, 1.0, denom)
    norm_traces = traces / denom

    relDiffVals = result.relDiff[ys, xs]
    order = np.argsort(relDiffVals)
    cmap = plt.get_cmap("RdBu_r")
    vmax = max(1e-6, float(np.abs(relDiffVals).max()))

    for k in order:
        c = cmap(0.5 + 0.5 * relDiffVals[k] / vmax)
        ax.plot(t, norm_traces[:, k], color=c, alpha=alpha, lw=1.0)

    ax.axvline(result.midRead - 0.5, color="k", lw=0.8, ls="--",
               label=f"split at read {result.midRead}")
    ax.set_xlabel("read index")
    ax.set_ylabel(ylabel)
    ax.set_title(
        f"outlier traces  cam={result.cam}  doLinearize={result.doLinearize}  "
        f"n={len(ys)}  |relDiff|>{thresh:g}"
    )
    ax.legend(loc="best", fontsize=9)

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=-vmax, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("relDiff")

    return fig, coords


def plotComparison(
    cmp: ComparisonResult,
    *,
    fig=None,
    sample: int = 200_000,
    vmin: float = -0.05,
    vmax: float = 0.05,
):
    """2x3 grid: linearize-off on top, linearize-on on bottom."""
    import matplotlib.pyplot as plt

    if fig is None:
        fig = plt.figure(figsize=(16, 10))
    subfigs = fig.subfigures(2, 1)

    plotHalves(
        cmp.off, fig=subfigs[0], sample=sample, vmin=vmin, vmax=vmax,
        title=f"linearize=OFF  cam={cmp.off.cam}  nReads={cmp.off.nReads}",
    )
    plotHalves(
        cmp.on, fig=subfigs[1], sample=sample, vmin=vmin, vmax=vmax,
        title=f"linearize=ON   cam={cmp.on.cam}  nReads={cmp.on.nReads}",
    )

    d = cmp.delta
    fig.suptitle(
        f"H4 linearity halves test — MAD off/on = {d['madImprovementFactor']:.2f}×  "
        f"(off={d['mad_off']:.4f}, on={d['mad_on']:.4f})",
        fontsize=12,
    )
    return fig
