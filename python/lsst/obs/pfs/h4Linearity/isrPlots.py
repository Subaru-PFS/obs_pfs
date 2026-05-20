"""Plotting helpers for H4 ISR validation over arbitrary read spans.

Currently covers:

- Half-vs-half ISR validation (data shape produced by ``/tmp/run_halves_isr.py``).
  Functions: `load`, `printSummary`, `plotAdditivityMap`, `plotRateComparisonMap`,
  `plotAmpBias`, `plotResidualHist`, `randomPixels`.
- Per-pixel ramp inspection across one or more read spans (single full ramp, both
  halves overlaid, or any other split). Functions: `plotPixelRamp`, with data
  collected by `validate.collectPixelRampData`.

Pure numpy + matplotlib — no LSST stack imports, so this module can be loaded
in any Python environment that has access to the saved arrays.

Typical use from JupyterLab::

    from lsst.obs.pfs.h4Linearity import isrPlots

    data = isrPlots.load('/work/cloomis/outputs/halves_isr_v142109_n3.npz')

    isrPlots.printSummary(data)
    isrPlots.plotAdditivityMap(data, corner='top-left')
    isrPlots.plotRateComparisonMap(data, rateMin=5.0)
    isrPlots.plotAmpBias(data)

Each plot function returns the ``Figure`` for further customization.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Union

import numpy as np


# H4 readout has 32 horizontal channels (slow direction), each 128 pixels wide.
H4_AMP_WIDTH = 128
H4_N_AMPS = 32


@dataclass
class PixelRampData:
    """Per-read cubes aligned for a per-pixel ramp inspection.

    All cubes share shape ``(N, H, W)`` and are indexed by the same
    ``readIndices`` along axis 0. ``readIndices[i]`` is the absolute read
    index in the original ramp; ``cube*[i, y, x]`` is the cumulative
    value at that read for pixel (y, x).

    cubeRaw : pre-linearization, post-dark-subtraction cumulative (the
        input to ``h4Linearity.apply``).
    cubeLin : post-linearization cumulative (the output of
        ``h4Linearity.apply``, after re-anchoring at firstRead).
    cubeDark : dark-cube cumulative aligned with the same reads.
    cubePreDark : pre-dark-subtraction cumulative = cubeRaw + cubeDark.
        Useful for seeing the raw detector reads before any subtraction.

    For ``firstRead=0`` (full ramp) and the default ``applyUTRWeights=True``,
    ``cubeRaw[i]`` corresponds to absolute read ``i + 1`` of the original
    ramp (the first delta is read 0 -> read 1). For ``firstRead > 0`` the
    absolute read at index ``i`` is ``firstRead + i + 1``.
    """

    cubeRaw: np.ndarray
    cubeLin: np.ndarray
    cubeDark: np.ndarray
    cubePreDark: np.ndarray
    readIndices: np.ndarray
    firstRead: int
    lastRead: int
    visit: int = -1
    cam: str = "?"
    fitMin: Optional[np.ndarray] = None
    fitMax: Optional[np.ndarray] = None
    # Post-CR-repair linearized cumulative. None when CR detection didn't run.
    # Differs from cubeLin at pixels repaired by cr.iterativeUtrDetectAndRepair.
    cubeCR: Optional[np.ndarray] = None
    # Full ``exposure.mask.array`` after ISR (incl. CR + ASIC_GLITCH stamps).
    mask: Optional[np.ndarray] = None
    # Mask plane dictionary {name: bit_index} captured at construction time.
    # Used by `randomPixels(..., plane="CR")` etc. to look up plane bits
    # without importing afw at plot time.
    maskPlaneDict: Optional[dict] = None
    # Per-pixel UTR rate (ADU/read). From the iterative CR detector
    # (post-repair) when doCR=True, else computed from cubeLin.
    avgRate: Optional[np.ndarray] = None
    # Per-delta flag arrays from the iterative CR detector, shape
    # ``(N-1, H, W)`` bool. ``crFlagMask[k, y, x]`` True means the delta
    # from read k to read k+1 was flagged as a CR at pixel (x, y). Used
    # by `plotPixelRamp` to mark flagged reads along the top of each row.
    crFlagMask: Optional[np.ndarray] = None
    glitchFlagMask: Optional[np.ndarray] = None


@dataclass
class HalvesIsrData:
    """Per-pixel arrays from a half-vs-half ISR run.

    All arrays are ``(H, W)`` float32 unless noted. ``mask*`` are integer
    bitfields from `lsst.afw.image.Mask`.
    """

    img1: np.ndarray
    img2: np.ndarray
    imgF: np.ndarray
    rate1: np.ndarray
    rate2: np.ndarray
    addResid: np.ndarray
    addResidRel: np.ndarray
    relDiff: np.ndarray
    avgRate: np.ndarray
    mask_first: np.ndarray
    mask_second: np.ndarray
    mask_full: np.ndarray
    maskUnion: np.ndarray
    visit: int
    cam: str
    midRead: int
    nReads: int


def load(path: str) -> HalvesIsrData:
    """Load an .npz produced by the halves-ISR validation script."""
    with np.load(path, allow_pickle=False) as f:
        # Scalars come back as 0-d arrays; unwrap to Python ints/strings.
        return HalvesIsrData(
            img1=f['img1'], img2=f['img2'], imgF=f['imgF'],
            rate1=f['rate1'], rate2=f['rate2'],
            addResid=f['addResid'], addResidRel=f['addResidRel'],
            relDiff=f['relDiff'], avgRate=f['avgRate'],
            mask_first=f['mask_first'], mask_second=f['mask_second'],
            mask_full=f['mask_full'], maskUnion=f['maskUnion'],
            visit=int(f['visit']), cam=str(f['cam']),
            midRead=int(f['midRead']), nReads=int(f['nReads']),
        )


_RATE_BINS = (
    (-1e9, 0.0), (0.0, 1.0), (1.0, 5.0), (5.0, 20.0),
    (20.0, 100.0), (100.0, 1000.0), (1000.0, 5000.0), (5000.0, 50000.0),
)


def _asData(data: Union[str, HalvesIsrData]) -> HalvesIsrData:
    return load(data) if isinstance(data, str) else data


def printSummary(data: Union[str, HalvesIsrData]) -> None:
    """Print the same kind of breakdown the validation script prints."""
    d = _asData(data)
    unmasked = (d.maskUnion == 0) & np.isfinite(d.relDiff)
    nU = int(unmasked.sum())
    total = d.maskUnion.size
    print(f"visit={d.visit} cam={d.cam}  nReads={d.nReads}  midRead={d.midRead}")
    print(f"unmasked: {nU:,}/{total:,} ({100 * nU / total:.2f}%)")

    print("\n=== |addResidRel| over unmasked, by avgRate ===")
    print(f"{'rate range':>14s}  {'nPix':>10s}  {'pct50':>9s}  {'pct95':>9s}  "
          f"{'pct99.9':>9s}  {'max':>9s}")
    for lo, hi in _RATE_BINS:
        sel = unmasked & (d.avgRate > lo) & (d.avgRate <= hi)
        n = int(sel.sum())
        if n == 0:
            continue
        rr = np.abs(d.addResidRel[sel])
        p50, p95, p999 = np.percentile(rr, [50, 95, 99.9])
        print(f"({lo:>5.0f},{hi:>6.0f}]  {n:>10,d}  {p50:>9.4f}  {p95:>9.4f}  "
              f"{p999:>9.4f}  {float(rr.max()):>9.4f}")

    print("\n=== |relDiff| over unmasked, by avgRate (rate-space test) ===")
    print(f"{'rate range':>14s}  {'nPix':>10s}  {'median':>9s}  {'MAD':>9s}  "
          f"{'|<1%':>7s}  {'|<5%':>7s}  {'|<20%':>7s}")
    for lo, hi in _RATE_BINS:
        sel = unmasked & (d.avgRate > lo) & (d.avgRate <= hi)
        n = int(sel.sum())
        if n == 0:
            continue
        rd = d.relDiff[sel]
        med = float(np.median(rd))
        mad = float(1.4826 * np.median(np.abs(rd - med)))
        f1 = float(np.mean(np.abs(rd) < 0.01))
        f5 = float(np.mean(np.abs(rd) < 0.05))
        f20 = float(np.mean(np.abs(rd) < 0.20))
        print(f"({lo:>5.0f},{hi:>6.0f}]  {n:>10,d}  {med:+8.4f}  {mad:8.4f}  "
              f"{f1:7.4f}  {f5:7.4f}  {f20:7.4f}")


def _cornerSlice(corner: str, height: int, width: int, size: int):
    """Return (yslice, xslice) for one of: top-left, top-right, bottom-left, bottom-right.

    Corner names follow the y-up display convention: "top" means high y,
    "bottom" means low y. Plots use ``origin='lower'``, so this matches the
    visual orientation in matplotlib.
    """
    if corner == 'top-left':
        return slice(height - size, height), slice(0, size)
    if corner == 'top-right':
        return slice(height - size, height), slice(width - size, width)
    if corner == 'bottom-left':
        return slice(0, size), slice(0, size)
    if corner == 'bottom-right':
        return slice(0, size), slice(width - size, width)
    raise ValueError(f"unknown corner {corner!r}; expected top-left / top-right / "
                     "bottom-left / bottom-right.")


def plotAdditivityMap(
    data: Union[str, HalvesIsrData],
    *,
    corner: Optional[str] = None,
    cornerSize: int = 600,
    vlim: float = 0.5,
    rateMin: float = 0.0,
    fig=None,
):
    """imshow of the additivity residual map.

    Parameters
    ----------
    data
        Loaded `HalvesIsrData` or path to the .npz.
    corner : {None, 'top-left', 'top-right', 'bottom-left', 'bottom-right'}
        If given, zoom to a ``cornerSize x cornerSize`` corner of the detector.
    cornerSize : int
        Side length of the corner zoom in pixels.
    vlim : float
        Color-scale limit; map is clipped to ``[-vlim, +vlim]``.
    rateMin : float
        Set pixels with ``avgRate <= rateMin`` to NaN so they show as gray
        (their relative residual is dominated by noise).
    fig : matplotlib.figure.Figure, optional
        Reuse an existing figure.
    """
    import matplotlib.pyplot as plt

    d = _asData(data)
    H, W = d.addResidRel.shape
    arr = d.addResidRel.copy()
    if rateMin > 0:
        arr[d.avgRate <= rateMin] = np.nan
    arr[d.maskUnion != 0] = np.nan

    if corner is not None:
        ys, xs = _cornerSlice(corner, H, W, cornerSize)
        arr = arr[ys, xs]
        extent = (xs.start, xs.stop, ys.start, ys.stop)
        title = (f"addResidRel — {corner}  (cornerSize={cornerSize}, "
                 f"rateMin={rateMin})")
    else:
        extent = None
        title = f"addResidRel — full detector  (rateMin={rateMin})"

    if fig is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        ax = fig.add_subplot(111)
    im = ax.imshow(arr, origin='lower', cmap='RdBu_r',
                   vmin=-vlim, vmax=vlim,
                   extent=extent, interpolation='nearest')
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label='resid / |full|')
    ax.set_title(f"visit={d.visit} {d.cam}  {title}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.tight_layout()
    return fig


def plotRateComparisonMap(
    data: Union[str, HalvesIsrData],
    *,
    rateMin: float = 5.0,
    vlim: float = 0.2,
    corner: Optional[str] = None,
    cornerSize: int = 600,
    fig=None,
):
    """imshow of ``relDiff = 2(rate1-rate2)/(rate1+rate2)``, masked to bright pixels.

    Parameters
    ----------
    data
        Loaded `HalvesIsrData` or path to the .npz.
    rateMin : float
        Pixels with ``avgRate <= rateMin`` are masked (set to NaN). Default
        5 ADU/read selects pixels with reasonable per-half S/N.
    vlim : float
        Color-scale clip ``[-vlim, +vlim]``.
    corner, cornerSize
        Same semantics as `plotAdditivityMap`.
    """
    import matplotlib.pyplot as plt

    d = _asData(data)
    H, W = d.relDiff.shape
    arr = d.relDiff.copy()
    bad = (d.maskUnion != 0) | (d.avgRate <= rateMin) | ~np.isfinite(arr)
    arr[bad] = np.nan

    if corner is not None:
        ys, xs = _cornerSlice(corner, H, W, cornerSize)
        arr = arr[ys, xs]
        extent = (xs.start, xs.stop, ys.start, ys.stop)
        zoom = f"  ({corner}, cornerSize={cornerSize})"
    else:
        extent = None
        zoom = ""

    if fig is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        ax = fig.add_subplot(111)
    im = ax.imshow(arr, origin='lower', cmap='RdBu_r',
                   vmin=-vlim, vmax=vlim,
                   extent=extent, interpolation='nearest')
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02,
                 label='2*(r1-r2)/(r1+r2)')
    ax.set_title(f"visit={d.visit} {d.cam}  relDiff  rateMin={rateMin}{zoom}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.tight_layout()
    return fig


def plotAmpBias(
    data: Union[str, HalvesIsrData],
    *,
    rateMin: float = 5.0,
    ampWidth: int = H4_AMP_WIDTH,
    fig=None,
):
    """Per-amplifier breakdown of the rate-comparison bias.

    Splits the detector into vertical stripes of width ``ampWidth`` (default
    128 for H4, giving 32 amps along x) and reports per-amp:

    - median(relDiff) over unmasked bright pixels
    - MAD(relDiff)
    - npix in the comparison

    Useful for spotting whether a global +N% bias is detector-wide or
    concentrated in particular amps.
    """
    import matplotlib.pyplot as plt

    d = _asData(data)
    H, W = d.relDiff.shape
    nAmps = W // ampWidth
    if nAmps * ampWidth != W:
        raise ValueError(
            f"detector width {W} not divisible by ampWidth {ampWidth}"
        )

    unmasked = (d.maskUnion == 0) & np.isfinite(d.relDiff) & (d.avgRate > rateMin)

    ampMedians = np.zeros(nAmps, dtype=np.float64)
    ampMads = np.zeros(nAmps, dtype=np.float64)
    ampN = np.zeros(nAmps, dtype=np.int64)
    for k in range(nAmps):
        xs = slice(k * ampWidth, (k + 1) * ampWidth)
        sel = unmasked[:, xs]
        rd = d.relDiff[:, xs][sel]
        ampN[k] = rd.size
        if rd.size == 0:
            ampMedians[k] = np.nan
            ampMads[k] = np.nan
            continue
        med = float(np.median(rd))
        ampMedians[k] = med
        ampMads[k] = float(1.4826 * np.median(np.abs(rd - med)))

    if fig is None:
        fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    else:
        axes = fig.subplots(2, 1, sharex=True)

    ampIdx = np.arange(nAmps)
    axes[0].errorbar(ampIdx, ampMedians, yerr=ampMads, fmt='o', capsize=3,
                     color='C0', label='median ± MAD')
    axes[0].axhline(0.0, color='k', lw=0.5, ls=':')
    overallMed = float(np.median(d.relDiff[unmasked]))
    axes[0].axhline(overallMed, color='C3', lw=0.8, ls='--',
                    label=f'global median = {overallMed:+.4f}')
    axes[0].set_ylabel('relDiff (rate1 vs rate2)')
    axes[0].set_title(f"visit={d.visit} {d.cam}  per-amp rate bias  "
                      f"(rateMin={rateMin}, ampWidth={ampWidth})")
    axes[0].legend(loc='best', fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].bar(ampIdx, ampN, color='0.5')
    axes[1].set_xlabel(f'amplifier index (x = [k*{ampWidth}, (k+1)*{ampWidth}))')
    axes[1].set_ylabel('nPix in comparison')
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    return fig


_PIXEL_RAMP_LINESTYLES = ('-', '--', ':', '-.')


# Mask planes to highlight with colored badges on per-pixel ramp plots.
# Order is priority left-to-right; any plane set on the pixel but not in
# this dict is rendered with the fallback grey color and trailing.
_MASK_PLANE_BADGES = (
    ('CR', '#d62728'),  # red
    ('ASIC_GLITCH', '#ff7f0e'),  # orange
    ('SAT', '#bcbd22'),  # ochre
    ('BAD', '#444444'),
    ('INTRP', '#888888'),
    ('NO_DATA', '#bbbbbb'),
    ('SUSPECT', '#9467bd'),  # purple
    ('EDGE', '#17becf'),  # cyan
)


def _activeMaskPlanes(maskPixel: int, planeDict: dict) -> list:
    """Return the names of mask planes set in ``maskPixel``, in badge order.

    Known planes appear in ``_MASK_PLANE_BADGES`` priority order;
    unrecognized set bits are appended in alphabetical order.
    """
    if maskPixel == 0:
        return []
    knownNames = {name for name, _ in _MASK_PLANE_BADGES}
    active = [
        name for name, _ in _MASK_PLANE_BADGES
        if name in planeDict and (maskPixel >> planeDict[name]) & 1
    ]
    for name in sorted(planeDict):
        if name in knownNames:
            continue
        if (maskPixel >> planeDict[name]) & 1:
            active.append(name)
    return active


def _drawMaskBadges(ax, planeNames: list, yPos: float = 0.97) -> None:
    """Render a row of colored mask-plane badges along the top of ``ax``."""
    if not planeNames:
        return
    colorMap = dict(_MASK_PLANE_BADGES)
    # Place left-to-right just inside the top edge of the axes.
    xPos = 0.01
    for name in planeNames:
        color = colorMap.get(name, '#666666')
        ax.text(
            xPos, yPos, name,
            transform=ax.transAxes, ha='left', va='top',
            fontsize=7, color='white', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.18',
                      facecolor=color, edgecolor='none', alpha=0.92),
        )
        # Advance xPos by the rendered width of the badge (in axes coords).
        # matplotlib doesn't give us the width until draw time; approximate
        # by character count + a small pad. Each char ≈ 0.012 axes units at
        # fontsize=7 on a typical figure width.
        xPos += 0.014 * (len(name) + 2)
        if xPos > 0.85:   # wrap is unlikely; just stop drawing instead
            break


def topDiscrepantPixels(
    data: "HalvesIsrData",
    n: int = 20,
    *,
    metric: str = "relDiffWeighted",
    rateMin: float = 5.0,
    rateMax: Optional[float] = None,
    unmasked: bool = True,
):
    """Return the top-N most discrepant pixels by the chosen metric.

    Avoids the noise-amplification problem of sorting on
    ``|addResid|/|imgF|``: that ratio blows up when ``|imgF|`` is small,
    so the top pixels tend to be low-S/N rather than truly anomalous.

    Parameters
    ----------
    data : HalvesIsrData
    n : int
        Number of pixels to return.
    metric : {'relDiff', 'relDiffWeighted', 'addResid'}
        Sort key:

        - ``relDiff`` — pure rate-comparison: ``|2(r1-r2)/(r1+r2)|``.
          Sensitive to any half-vs-half mismatch but still noisy for
          low-rate pixels (Poisson + read noise dominate).
        - ``relDiffWeighted`` (default) — ``|relDiff| * sqrt(avgRate)``.
          Approximate S/N weighting that prefers high-rate pixels with
          large rate disagreement.
        - ``addResid`` — raw additivity residual in ADU:
          ``|imgF - (img1 + img2)|``. Picks the largest absolute
          discrepancies. Best when you want bright outliers regardless
          of relative scale.
    rateMin : float
        Drop pixels with ``avgRate <= rateMin`` (default 5 ADU/read).
        Faint pixels are excluded from the ranking because their metric
        values are dominated by statistical noise.
    rateMax : float, optional
        Drop pixels with ``avgRate > rateMax``. Use this to focus on a
        specific flux regime.
    unmasked : bool
        Restrict to pixels with no mask bits set on any of the three
        ISR passes (default True).

    Returns
    -------
    coords : list[tuple[int, int]]
        Pixel coordinates ordered by *descending* metric value.
    """
    sel = np.ones(data.relDiff.shape, dtype=bool)
    if unmasked:
        sel &= (data.maskUnion == 0)
    sel &= np.isfinite(data.relDiff) & np.isfinite(data.avgRate)
    sel &= data.avgRate > rateMin
    if rateMax is not None:
        sel &= data.avgRate <= rateMax

    if metric == "relDiff":
        score = np.abs(data.relDiff)
    elif metric == "relDiffWeighted":
        score = np.abs(data.relDiff) * np.sqrt(np.maximum(data.avgRate, 0.0))
    elif metric == "addResid":
        score = np.abs(data.addResid)
    else:
        raise ValueError(
            f"unknown metric {metric!r}; expected "
            "'relDiff', 'relDiffWeighted', or 'addResid'."
        )
    score = np.where(sel, score, -np.inf)

    flat = score.ravel()
    if not np.isfinite(flat).any():
        raise RuntimeError(
            f"No pixels match the selection (rateMin={rateMin}, "
            f"rateMax={rateMax}, unmasked={unmasked})."
        )
    nReturn = min(int(n), int(np.isfinite(flat).sum()))
    # argpartition gets the top-nReturn without sorting the rest; then
    # sort that small slice in descending order for a clean ranking.
    top_idx = np.argpartition(flat, -nReturn)[-nReturn:]
    top_idx = top_idx[np.argsort(flat[top_idx])[::-1]]
    ys, xs = np.unravel_index(top_idx, score.shape)
    return list(zip(xs.tolist(), ys.tolist()))


_DEFAULT_EXCLUDE_PLANES = ('BAD', 'SAT', 'BORDER')


def randomPixels(
    data,
    n: int = 20,
    *,
    plane: Optional[str] = None,
    excludePlanes=_DEFAULT_EXCLUDE_PLANES,
    rateMin: Optional[float] = None,
    rateMax: Optional[float] = None,
    rng=None,
):
    """Return up to ``n`` random ``(x, y)`` pixel coordinates passing the selection.

    Accepts either a `HalvesIsrData` (uses ``avgRate`` + ``maskUnion``) or
    a `PixelRampData` (uses ``mask`` + ``maskPlaneDict``).

    Parameters
    ----------
    data : HalvesIsrData or PixelRampData
    n : int
        Maximum number of pixels to return. Fewer are returned if the
        selection has fewer than ``n`` candidates.
    plane : str, optional
        Restrict to pixels where this mask plane bit is set in the mask.
        Useful for sampling e.g. ``"CR"`` or ``"ASIC_GLITCH"`` pixels.
        Requires ``maskPlaneDict`` on the data (carried by
        ``PixelRampData`` from ``validate.collectPixelRampData``).
        Composes with ``excludePlanes``, ``rateMin``, ``rateMax``.
    excludePlanes : sequence of str
        Drop pixels with any of these mask planes set. Default
        ``('BAD', 'SAT', 'BORDER')`` skips obviously unusable pixels but
        keeps CR / ASIC_GLITCH ones. Pass ``()`` to disable. Names that
        aren't in ``data.maskPlaneDict`` are silently ignored, so this
        default works across detectors. If the data carries a mask but
        no plane dict (e.g. ``HalvesIsrData``), falls back to requiring
        all mask bits to be zero.
    rateMin, rateMax : float, optional
        Restrict to ``rateMin <= avgRate <= rateMax``. Requires the data
        to carry ``avgRate``.

        Common shortcuts:
          - faint pixels: ``rateMin=1.0, rateMax=10.0``
          - bright pixels: ``rateMin=50.0`` (no upper bound)
          - dark pixels: ``rateMax=1.0``
    rng : np.random.Generator or int, optional
        Reproducible random state. ``int`` → seeded generator; ``None`` →
        fresh generator each call.

    Returns
    -------
    coords : list[tuple[int, int]]
    """
    if isinstance(rng, (int, np.integer)):
        rng = np.random.default_rng(int(rng))
    elif rng is None:
        rng = np.random.default_rng()

    # Mask source: PixelRampData.mask (preferred) or HalvesIsrData.maskUnion.
    mask = getattr(data, 'mask', None)
    if mask is None:
        mask = getattr(data, 'maskUnion', None)
    planeDict = getattr(data, 'maskPlaneDict', None)
    avgRate = getattr(data, 'avgRate', None)

    if mask is not None:
        shape = mask.shape
    elif avgRate is not None:
        shape = avgRate.shape
    else:
        raise ValueError(
            "data has neither a `mask`/`maskUnion` nor `avgRate`; nothing "
            "to filter on."
        )
    sel = np.ones(shape, dtype=bool)

    if plane is not None:
        if mask is None or planeDict is None:
            raise ValueError(
                "plane= filtering requires `mask` + `maskPlaneDict` on the data "
                "(use a PixelRampData from validate.collectPixelRampData)."
            )
        if plane not in planeDict:
            raise ValueError(
                f"plane {plane!r} not in maskPlaneDict; available planes: "
                f"{sorted(planeDict.keys())}"
            )
        bit = mask.dtype.type(1) << mask.dtype.type(int(planeDict[plane]))
        sel &= (mask & bit) != 0

    # Drop "really bad" pixels by default (BAD / SAT / BORDER). When the
    # data has a plane dict, only the named bits are checked, so CR /
    # ASIC_GLITCH pixels still pass through.
    if excludePlanes:
        if mask is not None and planeDict is not None:
            excludeBit = mask.dtype.type(0)
            for name in excludePlanes:
                if name in planeDict:
                    excludeBit = excludeBit | (
                        mask.dtype.type(1) << mask.dtype.type(int(planeDict[name]))
                    )
            if excludeBit != 0:
                sel &= (mask & excludeBit) == 0
        elif mask is not None:
            # No plane dict (e.g. HalvesIsrData) — fall back to "no mask
            # bits at all". Preserves the previous default behavior on
            # that data type.
            sel &= (mask == 0)
    if avgRate is not None:
        sel &= np.isfinite(avgRate)

    if rateMin is not None or rateMax is not None:
        if avgRate is None:
            raise ValueError(
                "rateMin/rateMax require `avgRate` on the data."
            )
        if rateMin is not None:
            sel &= avgRate >= rateMin
        if rateMax is not None:
            sel &= avgRate <= rateMax

    ys, xs = np.where(sel)
    if ys.size == 0:
        raise RuntimeError(
            f"No pixels match the selection "
            f"(plane={plane!r}, excludePlanes={tuple(excludePlanes)}, "
            f"rateMin={rateMin}, rateMax={rateMax})."
        )
    n = min(int(n), ys.size)
    pick = rng.choice(ys.size, size=n, replace=False)
    return list(zip(xs[pick].tolist(), ys[pick].tolist()))


def plotPixelRamp(
    data,
    coords,
    *,
    includePreDark: bool = True,
    includeDark: bool = True,
    includeCR: bool = False,
    showFitRange: bool = False,
    addResidualColumn: bool = False,
    rowHeight: float = 2.0,
    width: float = 10.0,
    fig=None,
):
    """Plot per-read cumulative values for one or more pixels, one per row.

    Up to four curves per pixel per data source (read-by-read cumulative ADU):

    - ``cubePreDark`` — raw cumulative before dark subtraction (blue).
      Skipped if ``includePreDark=False``.
    - ``cubeDark`` — dark-cube cumulative (green). Skipped if ``includeDark=False``.
    - ``cubeRaw`` — post-dark, pre-linearization cumulative (orange). The input
      to ``h4Linearity.apply``.
    - ``cubeLin`` — linearized cumulative (red).

    Multi-source view: pass a list of `PixelRampData` (e.g. one for the first
    half of the ramp and one for the second) to see how the pixels evolve
    across both halves on the same axes. Each data source gets its own line
    style (solid, dashed, dotted, dash-dot); colors stay tied to curve type.
    ``cubePreDark`` / ``cubeDark`` / ``cubeRaw`` are continuous across the
    half boundary because they are not re-anchored; ``cubeLin`` IS re-anchored
    per source, so each source's linearized curve starts near zero at its own
    ``firstRead``.

    Parameters
    ----------
    data : `PixelRampData` or sequence of `PixelRampData`
        Single ramp (one set of curves) or multi-source list (one set of
        curves per source, different line styles).
    coords : (x, y) tuple, or sequence of (x, y) tuples
    includePreDark, includeDark : bool
        Toggle the pre-dark and dark curves.
    includeCR : bool
        If True and ``data.cubeCR`` is populated, overlay the
        CR-repaired linearized cumulative as a fifth curve (purple).
        Lets you see the difference CR repair makes per pixel.
    showFitRange : bool
        If True (and ``data.fitMin`` / ``data.fitMax`` are populated),
        draw horizontal dotted lines at the linearity fit limits and
        annotate the values. Default False because fitMax is usually
        far above the data range, which forces the y-axis to expand
        and shrinks the actual ramp curves.
    addResidualColumn : bool
        If True, add a second column showing each source's ``cubeLin``
        minus its own linear fit, plotted against the within-half read
        offset (``reads - reads[0]``). With multi-source data
        (first/second halves), the two residual curves overlay and
        should be ~flat and overlapping where linearization is working;
        systematic curvature signals a linearization failure in the
        relevant flux regime.
    rowHeight, width : float
        Figure sizing. ``rowHeight=2.0`` keeps a 20-pixel stack to 40 in tall.
    fig : matplotlib.figure.Figure, optional

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # Normalize `data` to a list of sources.
    if isinstance(data, PixelRampData):
        dataList = [data]
    else:
        dataList = list(data)
    if not dataList:
        raise ValueError("data is empty.")
    H, W = dataList[0].cubeLin.shape[1:]
    for d in dataList[1:]:
        if d.cubeLin.shape[1:] != (H, W):
            raise ValueError("all PixelRampData must share the same (H, W).")

    # Normalize `coords` to a list of (x, y) tuples.
    if (
        isinstance(coords, tuple)
        and len(coords) == 2
        and not isinstance(coords[0], (tuple, list))
    ):
        coords = [coords]
    coords = list(coords)
    nPx = len(coords)
    if nPx == 0:
        raise ValueError("coords is empty.")

    ncols = 2 if addResidualColumn else 1
    # Add a fixed header band on top of the row stack for the suptitle
    # plus a 2-row fontsize=8 legend (measured ~0.34 in tall with default
    # padding), plus gaps for the per-read markers that sit at the very
    # top edge of the first axes.
    headerHeightIn = 0.75
    figHeight = rowHeight * nPx + headerHeightIn
    # Figure width is fixed at ``width`` regardless of the column count;
    # an enabled residual column splits that width equally with the main
    # column. Tight inter-column spacing keeps each plot as wide as
    # possible at the cost of letting the columns share a y-axis edge.
    if fig is None:
        fig, axesGrid = plt.subplots(
            nPx, ncols, figsize=(width, figHeight),
            squeeze=False, sharex='col',
            gridspec_kw={'wspace': 0.06} if ncols > 1 else None,
        )
    else:
        axesGrid = fig.subplots(
            nPx, ncols, squeeze=False, sharex='col',
            gridspec_kw={'wspace': 0.06} if ncols > 1 else None,
        )
    axes = axesGrid[:, 0]
    residAxes = axesGrid[:, 1] if addResidualColumn else None

    # Use the first source for fitMin/fitMax (the linearity calibration is
    # per-pixel, not per-read, so it's the same regardless of source).
    fitMin_arr = dataList[0].fitMin
    fitMax_arr = dataList[0].fitMax
    showFitRange = (
        showFitRange and fitMin_arr is not None and fitMax_arr is not None
    )

    for i, (x, y) in enumerate(coords):
        ax = axes[i]
        for s, d in enumerate(dataList):
            ls = _PIXEL_RAMP_LINESTYLES[s % len(_PIXEL_RAMP_LINESTYLES)]
            reads = d.readIndices
            if includePreDark:
                ax.plot(reads, d.cubePreDark[:, y, x], 'o', linestyle=ls,
                        color='C0', lw=1.0, markersize=3)
            if includeDark:
                ax.plot(reads, d.cubeDark[:, y, x], 'd', linestyle=ls,
                        color='C2', lw=1.0, markersize=3)
            ax.plot(reads, d.cubeRaw[:, y, x], 's', linestyle=ls,
                    color='C1', lw=1.0, markersize=3)
            ax.plot(reads, d.cubeLin[:, y, x], '^', linestyle=ls,
                    color='C3', lw=1.0, markersize=3)
            if includeCR and d.cubeCR is not None:
                ax.plot(reads, d.cubeCR[:, y, x], 'v', linestyle=ls,
                        color='C4', lw=1.0, markersize=3)
            # Faint UTR-rate reference line: cubeLin[0] + rate * (k - k0).
            # Lets the eye spot deviations from the per-pixel linear
            # extrapolation that the detector uses internally.
            if d.avgRate is not None:
                rate = float(d.avgRate[y, x])
                anchor = float(d.cubeLin[0, y, x])
                utrLine = anchor + rate * (reads - reads[0]).astype(np.float64)
                ax.plot(reads, utrLine, linestyle=ls, color='0.4',
                        lw=0.8, alpha=0.45)
        if showFitRange:
            fmax = float(fitMax_arr[y, x])
            fmin = float(fitMin_arr[y, x])
            ax.axhline(fmax, color='0.4', lw=0.8, ls=':', alpha=0.8)
            ax.axhline(fmin, color='0.4', lw=0.8, ls=':', alpha=0.4)
            # Compact annotation in the upper-right corner.
            ax.text(0.995, 0.95,
                    f'fitMax={fmax:.0f}  fitMin={fmin:.0f}',
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=7, color='0.3',
                    bbox=dict(boxstyle='round,pad=0.2',
                              facecolor='white', edgecolor='none',
                              alpha=0.7))
        ax.set_ylabel('ADU', fontsize=9)
        # Pixel coordinate in the top-left of the plot (paste-ready Python
        # tuple-unpack form).
        ax.text(
            0.01, 0.97, f'(x, y) = ({x}, {y})',
            transform=ax.transAxes, ha='left', va='top',
            fontsize=8, fontweight='bold', color='black',
            bbox=dict(boxstyle='round,pad=0.2',
                      facecolor='white', edgecolor='none', alpha=0.85),
        )
        ax.grid(True, alpha=0.3)

        # Per-read CR / ASIC_GLITCH markers along the top of the axes.
        # A flag at delta index k → marker at read k+1 (the cumulative
        # sample where the jump first appears). Drawn per source so
        # multi-half data shows each half's flags at the matching reads.
        for d in dataList:
            reads = d.readIndices
            if d.crFlagMask is not None:
                flags = np.asarray(d.crFlagMask[:, y, x], dtype=bool)
                if flags.any():
                    ax.plot(
                        reads[1:][flags], np.full(int(flags.sum()), 1.02),
                        marker='v', color='#d62728', linestyle='none',
                        markersize=7, mew=0, clip_on=False,
                        transform=ax.get_xaxis_transform(),
                    )
            if d.glitchFlagMask is not None:
                flags = np.asarray(d.glitchFlagMask[:, y, x], dtype=bool)
                # Glitches come in delta-pairs (the up + down). The
                # physically corrupted read is the cumulative sample
                # between them — for a pair at delta indices (k, k+1)
                # that's read k+1. Mark each pair-start once.
                pairStart = flags[:-1] & flags[1:]
                if pairStart.any():
                    pairReads = reads[1:-1][pairStart]
                    ax.plot(
                        pairReads, np.full(len(pairReads), 1.02),
                        marker='v', mfc='none', color='#ff7f0e',
                        linestyle='none', markersize=8, mew=1.4,
                        clip_on=False,
                        transform=ax.get_xaxis_transform(),
                    )

        # Per-pixel badges for the non-temporal planes (BAD, SAT, INTRP,
        # NO_DATA, SUSPECT, EDGE, ...) — CR and ASIC_GLITCH are already
        # surfaced by the per-read markers above.
        unionMask = 0
        planeDict = None
        for d in dataList:
            if d.mask is not None and d.maskPlaneDict:
                unionMask |= int(d.mask[y, x])
                planeDict = d.maskPlaneDict
        if planeDict is not None and unionMask != 0:
            planes = [
                p for p in _activeMaskPlanes(unionMask, planeDict)
                if p not in ('CR', 'ASIC_GLITCH')
            ]
            # Place badges below the coord text, not on top of it.
            _drawMaskBadges(ax, planes, yPos=0.86)

        if residAxes is not None:
            axR = residAxes[i]
            # Joint linear fit across all sources: a single common slope
            # is subtracted from each source's cubeLin (no intercept). That
            # way endpoint *offsets* between halves remain visible -- they
            # would be exactly zero if linearization were perfectly
            # consistent across the ramp.
            xAll, yAll = [], []
            for d in dataList:
                xAll.append(
                    (d.readIndices - d.readIndices[0]).astype(np.float64)
                )
                yAll.append(d.cubeLin[:, y, x].astype(np.float64, copy=False))
            xAll = np.concatenate(xAll)
            yAll = np.concatenate(yAll)
            slope = float(np.polyfit(xAll, yAll, 1)[0]) if xAll.size >= 2 else 0.0
            for s, d in enumerate(dataList):
                ls = _PIXEL_RAMP_LINESTYLES[s % len(_PIXEL_RAMP_LINESTYLES)]
                xRel = (d.readIndices - d.readIndices[0]).astype(np.float64)
                lin = d.cubeLin[:, y, x].astype(np.float64, copy=False)
                resid = (lin - slope * xRel).astype(np.float32, copy=False)
                axR.plot(xRel.astype(np.int32), resid, '^', linestyle=ls,
                         color='C3', lw=1.0, markersize=3)
            axR.axhline(0.0, color='k', lw=0.5, ls=':', alpha=0.7)
            axR.grid(True, alpha=0.3)
            # Slope annotation inside the axes rather than as an ylabel,
            # so the residual plot's data area gets all available width.
            axR.text(
                0.01, 0.97, f'residual:  cubeLin − {slope:.1f}·x',
                transform=axR.transAxes, ha='left', va='top',
                fontsize=8, color='0.3',
                bbox=dict(boxstyle='round,pad=0.2',
                          facecolor='white', edgecolor='none', alpha=0.85),
            )
            # Residual y-ticks on the right so the inter-column gap
            # doesn't have to reserve space for tick labels.
            axR.yaxis.tick_right()
            axR.yaxis.set_label_position('right')

    # Combined legend on the top subplot: curve types (by color) and (only
    # when there are >1 sources) line styles per source.
    legendItems = []
    if includePreDark:
        legendItems.append(Line2D([0], [0], marker='o', color='C0',
                                  lw=1.0, label='pre-dark (raw)'))
    if includeDark:
        legendItems.append(Line2D([0], [0], marker='d', color='C2',
                                  lw=1.0, label='dark'))
    legendItems.append(Line2D([0], [0], marker='s', color='C1',
                              lw=1.0, label='post-dark (pre-lin input)'))
    legendItems.append(Line2D([0], [0], marker='^', color='C3',
                              lw=1.0, label='linearized'))
    if includeCR and any(d.cubeCR is not None for d in dataList):
        legendItems.append(Line2D([0], [0], marker='v', color='C4',
                                  lw=1.0, label='linearized + CR repaired'))
    if any(d.avgRate is not None for d in dataList):
        legendItems.append(Line2D([0], [0], color='0.4', lw=0.8, alpha=0.6,
                                  label='UTR rate × k'))
    if showFitRange:
        legendItems.append(Line2D([0], [0], color='0.4', lw=0.8, ls=':',
                                  label='fitMax / fitMin'))
    if any(d.crFlagMask is not None for d in dataList):
        legendItems.append(Line2D([0], [0], marker='v', color='#d62728',
                                  lw=0, markersize=7, mew=0,
                                  label='CR flag (read)'))
    if any(d.glitchFlagMask is not None for d in dataList):
        legendItems.append(Line2D([0], [0], marker='v', color='#ff7f0e',
                                  mfc='none', lw=0, markersize=8, mew=1.4,
                                  label='ASIC_GLITCH flag (read)'))
    if len(dataList) > 1:
        for s, d in enumerate(dataList):
            ls = _PIXEL_RAMP_LINESTYLES[s % len(_PIXEL_RAMP_LINESTYLES)]
            legendItems.append(Line2D(
                [0], [0], color='k', linestyle=ls, lw=1.2,
                label=f'src {s}: reads {d.firstRead + 1}..{d.lastRead}',
            ))
    # Figure-level legend in the header band — placed in figure coords so
    # it sits below the suptitle at a fixed offset regardless of row count.
    # Top edge ~0.26 in below the figure top leaves ~0.04 in gap below the
    # ~0.18 in tall suptitle text.
    fig.legend(
        handles=legendItems,
        loc='upper left',
        bbox_to_anchor=(0.07, 1.0 - 0.26 / figHeight),
        ncol=max(2, (len(legendItems) + 1) // 2),
        fontsize=8, framealpha=0.85, borderaxespad=0.0,
    )

    if len(dataList) == 1:
        d = dataList[0]
        xlabel = (
            f'absolute read index  '
            f'(cube spans firstRead={d.firstRead}..lastRead={d.lastRead})'
        )
    else:
        lo = min(d.firstRead for d in dataList)
        hi = max(d.lastRead for d in dataList)
        xlabel = (
            f'absolute read index  (sources span reads {lo + 1}..{hi})'
        )
    axes[-1].set_xlabel(xlabel)
    if residAxes is not None:
        residAxes[-1].set_xlabel('read offset from start of half')
        residAxes[0].set_title(
            'cubeLin − common slope × offset  (per-pixel joint fit)',
            fontsize=10,
        )
        axes[0].set_title('ramp curves', fontsize=10)

    visit = dataList[0].visit
    cam = dataList[0].cam
    suffix = '' if len(dataList) == 1 else f'  ({len(dataList)} sources)'
    # Suptitle anchored near the top edge of the figure; the legend below
    # is offset further to leave room for the ~0.18 in tall title text.
    fig.suptitle(f'visit={visit} {cam}  per-pixel ramps{suffix}',
                 fontsize=11, y=1.0 - 0.04 / figHeight)
    # Manual layout — tight_layout can't model the bbox-anchored external
    # legend (warns and overshoots the top reserve). subplots_adjust uses
    # the header band absolutely and keeps the per-row data area as wide
    # and as tall as possible.
    topMargin = 1.0 - headerHeightIn / figHeight
    fig.subplots_adjust(
        top=topMargin,
        bottom=0.05 if nPx == 1 else 0.04 + 0.04 / nPx,
        left=0.07,
        right=0.97,
        wspace=0.06,
        hspace=0.15,
    )
    return fig


def plotResidualHist(
    data: Union[str, HalvesIsrData],
    *,
    rateMin: float = 5.0,
    fig=None,
):
    """Two-panel histogram: ``relDiff`` and ``addResidRel`` over unmasked bright pixels.

    Useful for at-a-glance shape — symmetric? heavy-tailed? bias?
    """
    import matplotlib.pyplot as plt

    d = _asData(data)
    unmasked = (d.maskUnion == 0) & np.isfinite(d.relDiff) & (d.avgRate > rateMin)
    if not unmasked.any():
        raise RuntimeError(f"no unmasked pixels with avgRate > {rateMin}")

    rd = d.relDiff[unmasked]
    rr = d.addResidRel[unmasked]

    if fig is None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    else:
        axes = fig.subplots(1, 2)

    for ax, vals, lbl in [(axes[0], rd, 'relDiff'),
                          (axes[1], rr, 'addResid / |imgF|')]:
        clipped = vals[np.abs(vals) <= 0.5]
        med = float(np.median(vals))
        mad = float(1.4826 * np.median(np.abs(vals - med)))
        ax.hist(clipped, bins=120, color='0.3')
        ax.axvline(med, color='C3', lw=1.5,
                   label=f'median={med:+.4f}  MAD={mad:.4f}  N={vals.size:,}')
        ax.axvline(0.0, color='k', lw=0.5, ls=':')
        ax.set_xlabel(lbl)
        ax.set_ylabel(f'pixels with avgRate > {rateMin}')
        ax.legend(loc='upper right', fontsize=9)

    fig.suptitle(f"visit={d.visit} {d.cam}  half-vs-half residual histograms",
                 fontsize=11)
    fig.tight_layout()
    return fig


def randomCRPixels(result, n: int = 10, *, rng=None):
    """Return up to ``n`` random ``(x, y)`` pixels flagged as CR by ``result``.

    Parameters
    ----------
    result : `cr.IterativeRepairResult`
    n : int
    rng : np.random.Generator or int, optional
    """
    if isinstance(rng, (int, np.integer)):
        rng = np.random.default_rng(int(rng))
    elif rng is None:
        rng = np.random.default_rng()
    pix2D = result.crFlagMask.any(axis=0)
    ys, xs = np.where(pix2D)
    if ys.size == 0:
        raise RuntimeError("no CR-flagged pixels in result.")
    n = min(int(n), ys.size)
    pick = rng.choice(ys.size, size=n, replace=False)
    return list(zip(xs[pick].tolist(), ys[pick].tolist()))


def randomGlitchPixels(result, n: int = 10, *, rng=None):
    """Return up to ``n`` random ``(x, y)`` pixels flagged as ASIC glitch.

    Parameters
    ----------
    result : `cr.IterativeRepairResult`
    n : int
    rng : np.random.Generator or int, optional
    """
    if isinstance(rng, (int, np.integer)):
        rng = np.random.default_rng(int(rng))
    elif rng is None:
        rng = np.random.default_rng()
    pix2D = result.glitchFlagMask.any(axis=0)
    ys, xs = np.where(pix2D)
    if ys.size == 0:
        raise RuntimeError("no ASIC-glitch-flagged pixels in result.")
    n = min(int(n), ys.size)
    pick = rng.choice(ys.size, size=n, replace=False)
    return list(zip(xs[pick].tolist(), ys[pick].tolist()))


def plotCRRamps(
    cubeOriginal: np.ndarray,
    cubeRepaired: np.ndarray,
    result,
    coords,
    *,
    space: str = 'flux',
    rowHeight: float = 2.0,
    width: float = 6.0,
    fig=None,
):
    """Plot N pixels comparing pre- and post-repair ramps.

    One row per pixel, two columns:

    - Left: ``cubeOriginal`` (pre-repair linearized) with the per-pixel UTR
      rate from ``result.rate``. Flux space overlays the line
      ``cubeRepaired[0] + rate * k``; delta space overlays a horizontal
      line at ``rate``. The line is anchored on the repaired cube's first
      sample so a CR/glitch in an early delta does not offset it.
    - Right: ``cubeRepaired`` (post-CR/glitch repair) drawn the same way.

    Flagged samples are highlighted in both columns: CR as red ``x``,
    ASIC glitch as open orange ``o``. In flux space the marker sits on
    the cumulative *after* the flagged delta (index ``k+1``); in delta
    space it sits on the flagged delta itself.

    Parameters
    ----------
    cubeOriginal, cubeRepaired : np.ndarray
        ``(N, H, W)`` cumulative cubes (pre- and post-repair).
    result : `cr.IterativeRepairResult`
        Source of ``rate``, ``crFlagMask``, ``glitchFlagMask``.
    coords : (x, y) tuple, or sequence of tuples.
    space : {'flux', 'delta'}
    rowHeight, width : float
    fig : matplotlib.figure.Figure, optional
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    if space not in ('flux', 'delta'):
        raise ValueError(f"space must be 'flux' or 'delta'; got {space!r}.")
    if cubeOriginal.shape != cubeRepaired.shape:
        raise ValueError(
            f"cubeOriginal {cubeOriginal.shape} and cubeRepaired "
            f"{cubeRepaired.shape} must share shape."
        )
    if (
        isinstance(coords, tuple)
        and len(coords) == 2
        and not isinstance(coords[0], (tuple, list))
    ):
        coords = [coords]
    coords = list(coords)
    nPx = len(coords)
    if nPx == 0:
        raise ValueError("coords is empty.")

    if fig is None:
        fig, axesGrid = plt.subplots(
            nPx, 2, figsize=(width * 2, rowHeight * nPx),
            squeeze=False, sharex='col',
        )
    else:
        axesGrid = fig.subplots(nPx, 2, squeeze=False, sharex='col')

    nReads = cubeOriginal.shape[0]
    reads = np.arange(nReads)
    deltaIdx = np.arange(nReads - 1)

    for i, (x, y) in enumerate(coords):
        rate = float(result.rate[y, x])
        crFlag = np.asarray(result.crFlagMask[:, y, x], dtype=bool)
        glFlag = np.asarray(result.glitchFlagMask[:, y, x], dtype=bool)

        if space == 'flux':
            seriesOrig = cubeOriginal[:, y, x]
            seriesRep = cubeRepaired[:, y, x]
            # A flag at delta index k -> highlight cumulative sample k+1.
            crSamp = np.concatenate(([False], crFlag))
            glSamp = np.concatenate(([False], glFlag))
            xVals = reads
            xLabel = 'read index'
            yLabel = 'ADU'
        else:
            seriesOrig = np.diff(cubeOriginal[:, y, x])
            seriesRep = np.diff(cubeRepaired[:, y, x])
            crSamp = crFlag
            glSamp = glFlag
            xVals = deltaIdx
            xLabel = 'delta index'
            yLabel = 'ΔADU'

        for col, vals, title in ((0, seriesOrig, 'linearized'),
                                 (1, seriesRep, 'repaired')):
            ax = axesGrid[i, col]
            ax.plot(xVals, vals, '.-', color='C3', lw=1.0, markersize=4)
            if space == 'flux':
                # Anchor on the repaired cube's first sample: cubeOriginal[0]
                # is contaminated when an early delta is the flagged CR/glitch.
                ax.plot(reads, float(seriesRep[0]) + rate * reads,
                        color='k', lw=0.8, ls='--', alpha=0.7)
            else:
                ax.axhline(rate, color='k', lw=0.8, ls='--', alpha=0.7)
            if crSamp.any():
                ax.plot(xVals[crSamp], vals[crSamp], 'x',
                        color='red', markersize=8, mew=1.5, linestyle='none')
            if glSamp.any():
                ax.plot(xVals[glSamp], vals[glSamp], 'o',
                        mfc='none', color='orange', markersize=8, mew=1.5,
                        linestyle='none')
            if i == 0:
                ax.set_title(title)
            ax.grid(True, alpha=0.3)
        axesGrid[i, 0].set_ylabel(f'(x, y) = ({x}, {y})\n{yLabel}', fontsize=9)

    axesGrid[-1, 0].set_xlabel(xLabel)
    axesGrid[-1, 1].set_xlabel(xLabel)

    legendItems = [
        Line2D([0], [0], marker='.', color='C3', lw=1.0, label='cube'),
        Line2D([0], [0], color='k', lw=0.8, ls='--', label='UTR rate'),
        Line2D([0], [0], marker='x', color='red', lw=0, markersize=8,
               mew=1.5, label='CR'),
        Line2D([0], [0], marker='o', color='orange', mfc='none', lw=0,
               markersize=8, mew=1.5, label='ASIC glitch'),
    ]
    axesGrid[0, 0].legend(handles=legendItems, loc='upper left', fontsize=8,
                          framealpha=0.85)
    fig.suptitle(f'CR / ASIC-glitch correction  ({space} space)',
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    return fig


_DEFAULT_DS9_PLANE_COLORS = {
    'BAD': 'red',
    'SAT': 'yellow',
    'CR': 'magenta',
    'ASIC_GLITCH': 'orange',
    'INTRP': 'green',
    'NO_DATA': 'cyan',
    'EDGE': 'cyan',
    'BORDER': 'cyan',
    'SUSPECT': 'purple',
}


def displayDS9(
    image,
    mask=None,
    *,
    maskPlaneDict=None,
    planes=('BAD', 'SAT', 'CR', 'ASIC_GLITCH'),
    planeColors=None,
    ds9=None,
    transparency: int = 70,
    scale: str = 'zscale',
):
    """Display a 2D image plus colored mask planes on a pyds9 DS9 frame.

    Opens a new frame on the given (or newly created) DS9 connection,
    sends the image, then overlays each requested mask plane as its own
    colored mask layer.

    Parameters
    ----------
    image : np.ndarray
        2D image to display (e.g. UTR rate, CDS, single read).
    mask : np.ndarray, optional
        2D integer mask with the per-plane bits set. If None, only the
        image is displayed.
    maskPlaneDict : dict, optional
        ``{plane_name: bit_index}`` mapping. Required when ``mask`` is
        given. PixelRampData carries this as ``data.maskPlaneDict``.
    planes : sequence of str
        Mask plane names to overlay. Default ``('BAD', 'SAT', 'CR',
        'ASIC_GLITCH')``. Names missing from ``maskPlaneDict`` are
        silently skipped.
    planeColors : dict, optional
        Override the default color map. Built-in defaults: BAD=red,
        SAT=yellow, CR=magenta, ASIC_GLITCH=orange, INTRP=green,
        NO_DATA/EDGE/BORDER=cyan, SUSPECT=purple.
    ds9 : pyds9.DS9, optional
        Existing DS9 connection. If None, a new one is created (which
        launches DS9 if it isn't running).
    transparency : int
        Mask transparency, 0 (opaque) to 100 (invisible). DS9 default 70.
    scale : str
        DS9 scale spec, e.g. ``'zscale'``, ``'linear'``, ``'log'``.

    Returns
    -------
    ds9 : pyds9.DS9
        The (possibly newly created) DS9 connection. Pass back in to
        layer more frames on the same DS9 window.
    """
    try:
        import pyds9
    except ImportError as e:
        raise ImportError(
            "pyds9 is required for displayDS9; install it (pip install pyds9) "
            "or run via the LSST stack environment that provides it."
        ) from e

    if ds9 is None:
        ds9 = pyds9.DS9()

    ds9.set('frame new')

    img = np.asarray(image, dtype=np.float32)
    ds9.set_np2arr(img)
    ds9.set(f'scale {scale}')
    ds9.set('zoom to fit')

    if mask is None or maskPlaneDict is None or not planes:
        return ds9

    colors = {**_DEFAULT_DS9_PLANE_COLORS, **(planeColors or {})}
    ds9.set(f'mask transparency {int(transparency)}')

    maskArr = np.asarray(mask)
    H, W = maskArr.shape
    for plane in planes:
        if plane not in maskPlaneDict:
            continue
        bit = int(maskPlaneDict[plane])
        # Each plane is sent as its own DS9 mask layer; the most recent
        # `mask color` applies to the next array-mask push.
        planeArr = (maskArr & maskArr.dtype.type(1 << bit)).astype('u2')
        if planeArr.sum() == 0:
            continue
        color = colors.get(plane, 'gray')
        ds9.set(f'mask color {color}')
        ds9.set(
            f'array mask [xdim={W},ydim={H},bitpix=16]', planeArr,
        )

    return ds9


def _selMask(data, excludePlanes):
    """``(H, W)`` bool mask of pixels with NONE of ``excludePlanes`` set.

    Falls back to all-True when the data has no mask. Used by the
    halves-comparison helpers below.
    """
    rateShape = None
    if getattr(data, 'avgRate', None) is not None:
        rateShape = data.avgRate.shape
    if getattr(data, 'mask', None) is None or getattr(data, 'maskPlaneDict', None) is None:
        return np.ones(rateShape, dtype=bool) if rateShape is not None else None
    excludeBit = data.mask.dtype.type(0)
    for name in excludePlanes:
        if name in data.maskPlaneDict:
            excludeBit = excludeBit | (
                data.mask.dtype.type(1) << data.mask.dtype.type(int(data.maskPlaneDict[name]))
            )
    if excludeBit == 0:
        return np.ones(data.mask.shape, dtype=bool)
    return (data.mask & excludeBit) == 0


def _halvesRateMetric(first, second, metric: str):
    """Return ``(metricArr, xlabel)`` for the requested comparison metric."""
    r1 = np.asarray(first.avgRate, dtype=np.float64)
    r2 = np.asarray(second.avgRate, dtype=np.float64)
    if metric == 'diff':
        return r1 - r2, 'rate diff  (first − second)  [ADU/read]'
    if metric == 'reldiff':
        with np.errstate(invalid='ignore', divide='ignore'):
            m = 2.0 * (r1 - r2) / (r1 + r2)
        return m, 'rel rate diff  2(first − second)/(first + second)'
    if metric == 'ratio':
        with np.errstate(invalid='ignore', divide='ignore'):
            m = r1 / r2
        return m, 'rate ratio  first / second'
    raise ValueError(
        f"unknown metric {metric!r}; expected 'diff', 'reldiff', or 'ratio'."
    )


def plotHalvesRateHistogram(
    first,
    second,
    *,
    metric: str = 'diff',
    bins: int = 100,
    excludePlanes=_DEFAULT_EXCLUDE_PLANES,
    logy: bool = True,
    range_=None,
    ax=None,
):
    """Histogram of per-pixel UTR-rate comparison between two ramp halves.

    Parameters
    ----------
    first, second : PixelRampData
        Half-ramp results from ``validate.collectPixelRampData``.
    metric : {'diff', 'reldiff', 'ratio'}
        - ``'diff'``: ``first.avgRate − second.avgRate`` (ADU/read).
        - ``'reldiff'``: ``2*(first - second)/(first + second)``.
        - ``'ratio'``: ``first / second``.
    bins : int
    excludePlanes : sequence of str
        Drop pixels with any of these planes set in EITHER half.
        Defaults to ``('BAD', 'SAT', 'BORDER')``.
    logy : bool
    range_ : tuple, optional
        ``(lo, hi)`` x-axis range for the histogram. Default is set
        from the 1st / 99th percentile span of the metric.
    ax : matplotlib axis, optional

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    metricArr, xlabel = _halvesRateMetric(first, second, metric)
    sel1 = _selMask(first, excludePlanes)
    sel2 = _selMask(second, excludePlanes)
    valid = np.isfinite(metricArr) & sel1 & sel2
    vals = metricArr[valid]

    if ax is None:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.figure

    if range_ is None:
        if vals.size:
            p1, p99 = np.percentile(vals, [1, 99])
            if metric == 'ratio':
                pad = max(p99 - p1, 0.1)
                lo, hi = max(p1 - 0.5 * pad, 0.0), p99 + 0.5 * pad
            else:
                span = max(abs(p1), abs(p99))
                lo, hi = -2.5 * span, 2.5 * span
        else:
            lo, hi = -1.0, 1.0
        range_ = (float(lo), float(hi))

    ax.hist(
        vals, bins=bins, range=range_,
        color='steelblue', edgecolor='black', lw=0.4,
    )
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('pixels')
    ax.grid(True, alpha=0.3, axis='y')

    # Annotate percentiles of |deviation from center|.
    center = 1.0 if metric == 'ratio' else 0.0
    dev = np.abs(vals - center)
    if dev.size:
        pcts = (50, 84, 95, 99, 99.9)
        pVals = np.percentile(dev, pcts)
        pStr = '   '.join(f'{p}%={v:.3f}' for p, v in zip(pcts, pVals))
    else:
        pStr = 'no valid pixels'

    visit = getattr(first, 'visit', '?')
    cam = getattr(first, 'cam', '?')
    nValid = int(valid.sum())
    nExcl = int(valid.size - nValid)
    ax.set_title(
        f'v{visit} {cam}  half-vs-half UTR rate ({metric})  '
        f'[N={nValid:,}, excluded={nExcl:,}]\n'
        f'|Δ| pct  {pStr}',
        fontsize=10,
    )
    fig.tight_layout()
    return fig


def summarizeHalves(
    first,
    second,
    *,
    metric: str = 'diff',
    threshold: Optional[float] = None,
    excludePlanes=_DEFAULT_EXCLUDE_PLANES,
    maxList: int = 20,
    printOutput: bool = True,
):
    """Mask-plane counts + outlier coords from a halves comparison.

    Parameters
    ----------
    first, second : PixelRampData
    metric : {'diff', 'reldiff', 'ratio'}
    threshold : float, optional
        Outlier threshold for the chosen metric. Defaults: diff=5
        ADU/read, reldiff=0.2, ratio=2 (an outlier is also <1/threshold
        for ratio).
    excludePlanes : sequence of str
        Pixels with any of these mask planes set in EITHER half are
        excluded from outlier counting.
    maxList : int
        Cap on the number of outliers listed in the printed summary.
        The returned ``outliers`` list is full-length.
    printOutput : bool
        Print a text summary; defaults True.

    Returns
    -------
    dict with keys:
      - ``outliers``: list of ``(x, y)`` in descending |deviation| order
      - ``metricValues``: np.ndarray of the metric at each outlier,
        same order as ``outliers``
      - ``maskCounts``: ``{plane: (first, second, union)}`` pixel counts
      - ``percentiles``: ``{pct: |Δ|}`` over valid pixels
      - ``nValid``: count of pixels that passed the excludePlanes filter
    """
    metricArr, label = _halvesRateMetric(first, second, metric)
    sel1 = _selMask(first, excludePlanes)
    sel2 = _selMask(second, excludePlanes)
    valid = np.isfinite(metricArr) & sel1 & sel2

    if threshold is None:
        threshold = {'diff': 5.0, 'reldiff': 0.2, 'ratio': 2.0}.get(metric, 0.0)

    center = 1.0 if metric == 'ratio' else 0.0
    if metric == 'ratio':
        outlierMask = valid & (
            (metricArr > threshold) | (metricArr < 1.0 / max(threshold, 1e-9))
        )
    else:
        outlierMask = valid & (np.abs(metricArr - center) > threshold)

    ys, xs = np.where(outlierMask)
    diffs = metricArr[outlierMask]
    order = np.argsort(np.abs(diffs - center))[::-1]
    outliers = [(int(xs[i]), int(ys[i])) for i in order]
    metricValues = diffs[order]

    # Mask-plane counts per half + union.
    planeNames = set()
    for d in (first, second):
        if getattr(d, 'maskPlaneDict', None):
            planeNames |= set(d.maskPlaneDict)
    maskCounts = {}
    for name in sorted(planeNames):
        cf = cs = cu = 0
        bf = bs = 0
        if first.mask is not None and name in (first.maskPlaneDict or {}):
            bf = 1 << first.maskPlaneDict[name]
            cf = int(((first.mask & bf) != 0).sum())
        if second.mask is not None and name in (second.maskPlaneDict or {}):
            bs = 1 << second.maskPlaneDict[name]
            cs = int(((second.mask & bs) != 0).sum())
        if first.mask is not None and second.mask is not None:
            cu = int((((first.mask & bf) != 0) | ((second.mask & bs) != 0)).sum())
        maskCounts[name] = (cf, cs, cu)

    vals = metricArr[valid]
    dev = np.abs(vals - center)
    pcts = (50, 84, 95, 99, 99.9)
    percentiles = (
        {p: float(np.percentile(dev, p)) for p in pcts} if dev.size
        else {p: float('nan') for p in pcts}
    )

    if printOutput:
        nValid = int(valid.sum())
        nExcl = int(valid.size - nValid)
        print(f'half-vs-half {label}')
        print(f'  valid pixels (after excludePlanes={tuple(excludePlanes)}): '
              f'{nValid:,}   excluded: {nExcl:,}')
        print(f'  outlier threshold: {threshold}')
        print()
        print('mask-plane counts  (first / second / union):')
        for name, (cf, cs, cu) in maskCounts.items():
            print(f'  {name:14s}  {cf:>10,} / {cs:>10,} / {cu:>10,}')
        print()
        print('|Δ| percentiles over valid pixels:')
        for p in pcts:
            print(f'  {p:5}%: {percentiles[p]:.4f}')
        print()
        print(f'outliers (|metric − {center}| > {threshold}): {len(outliers):,}')
        if outliers:
            shown = min(maxList, len(outliers))
            print(f'top {shown} by |Δ|:')
            print(f'  {"(x, y)":>16}     metric')
            for (x, y), v in zip(outliers[:shown], metricValues[:shown]):
                print(f'  {f"({x}, {y})":>16}   {v:>8.3f}')

    return {
        'outliers': outliers,
        'metricValues': metricValues,
        'maskCounts': maskCounts,
        'percentiles': percentiles,
        'nValid': int(valid.sum()),
    }


def displayPixelRampDS9(data, image=None, *, ds9=None, **kwargs):
    """Send a `PixelRampData`'s mask + rate to DS9 in one call.

    Convenience wrapper around `displayDS9` that pulls ``mask`` and
    ``maskPlaneDict`` straight off ``data``.

    Parameters
    ----------
    data : PixelRampData
    image : np.ndarray or str, optional
        2D image to display. If None (default), uses ``data.avgRate``.
        If a string, looks up that attribute on ``data`` (handy for
        ``'avgRate'``; cube fields are 3D so pass a slice like
        ``data.cubeLin[10]`` directly instead).
    ds9 : pyds9.DS9, optional
        Existing DS9 connection; passes through to `displayDS9`.
    **kwargs : forwarded to `displayDS9` (planes, planeColors,
        transparency, scale, ...).

    Returns
    -------
    ds9 : pyds9.DS9
    """
    if image is None:
        image = data.avgRate
    elif isinstance(image, str):
        image = getattr(data, image)
    return displayDS9(
        image, data.mask,
        maskPlaneDict=data.maskPlaneDict,
        ds9=ds9,
        **kwargs,
    )


def filterCoordsByChannel(
    coords,
    channels,
    *,
    mode: str = 'select',
    channelHeight: int = H4_AMP_WIDTH,
):
    """Filter ``(x, y)`` pixel coordinates by ASIC channel.

    H4 detectors have 32 ASIC channels stacked along y, each
    ``channelHeight`` rows tall. The channel of pixel ``(x, y)`` is
    ``y // channelHeight``.

    Parameters
    ----------
    coords : sequence of ``(x, y)`` tuples
        Pixel coordinates in the project's user-facing convention.
    channels : int or sequence of ints
        Channel indices (0..31 for H4) to act on.
    mode : {'select', 'cull'}
        ``'select'`` (default) keeps coords whose channel is in
        ``channels``; ``'cull'`` drops them.
    channelHeight : int
        Rows per ASIC channel; 128 for H4.

    Returns
    -------
    list of ``(x, y)`` tuples, same order as the input.
    """
    if mode not in ('select', 'cull'):
        raise ValueError(f"mode must be 'select' or 'cull'; got {mode!r}.")
    if isinstance(channels, (int, np.integer)):
        channels = (int(channels),)
    channelSet = {int(c) for c in channels}
    keepInside = mode == 'select'
    out = []
    for x, y in coords:
        inChannels = (int(y) // channelHeight) in channelSet
        if inChannels == keepInside:
            out.append((int(x), int(y)))
    return out


def plotGlitchHistogramPerChannel(
    source,
    *,
    channelHeight: int = H4_AMP_WIDTH,
    countPixels: bool = False,
    ax=None,
    color: str = '#ff7f0e',
):
    """Bar histogram of ASIC glitches per H4 ASIC channel.

    H4 detectors have 32 ASIC channels along the y-axis, each
    ``channelHeight`` rows tall. This plot shows how many glitches the
    iterative detector flagged in each channel.

    Parameters
    ----------
    source : PixelRampData, cr.IterativeRepairResult, or np.ndarray
        Anything carrying a ``glitchFlagMask`` attribute of shape
        ``(N-1, H, W)``, or the array itself.
    channelHeight : int
        Rows per ASIC channel. 128 for H4.
    countPixels : bool
        Default (False) counts glitch PAIRS — distinct ASIC-glitch
        events. True counts unique flagged pixels per channel; useful
        for spatial-density inspection but obscures pixels that took
        multiple hits.
    ax : matplotlib axis, optional
    color : str

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    if isinstance(source, np.ndarray):
        flags = source
    else:
        flags = getattr(source, 'glitchFlagMask', None)
    if flags is None:
        raise ValueError(
            "source has no `glitchFlagMask` (pass an IterativeRepairResult, "
            "a PixelRampData populated by validate.collectPixelRampData, or "
            "the (N-1, H, W) bool array directly)."
        )
    flags = np.asarray(flags, dtype=bool)
    if flags.ndim != 3:
        raise ValueError(
            f"glitchFlagMask must be (N-1, H, W); got shape {flags.shape}."
        )

    if countPixels:
        # Per-pixel "had at least one glitch hit" mask, summed per channel.
        perPixel = flags.any(axis=0)
        ylabel = 'unique flagged pixels'
    else:
        # Glitch pairs come as two adjacent True deltas; count each pair
        # once at its starting delta index.
        pairStart = flags[:-1] & flags[1:]
        perPixel = pairStart.sum(axis=0).astype(np.int64)
        ylabel = 'ASIC-glitch pair count'

    H, W = perPixel.shape
    nChannels = H // channelHeight
    if H % channelHeight != 0:
        raise ValueError(
            f"H={H} not divisible by channelHeight={channelHeight}; "
            f"can't bin by channel cleanly."
        )
    perChannel = perPixel.reshape(nChannels, channelHeight, W).sum(axis=(1, 2))

    if ax is None:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.figure

    ax.bar(np.arange(nChannels), perChannel,
           color=color, edgecolor='black', lw=0.4)
    ax.set_xlabel(f'ASIC channel  (y // {channelHeight})')
    ax.set_ylabel(ylabel)
    titleSuffix = ''
    cam = getattr(source, 'cam', None)
    visit = getattr(source, 'visit', None)
    if cam and visit:
        titleSuffix = f'  (visit={visit} {cam})'
    ax.set_title(f'ASIC glitches per channel{titleSuffix}')
    ax.set_xticks(np.arange(nChannels))
    ax.tick_params(axis='x', labelsize=8)
    ax.grid(True, axis='y', alpha=0.3)
    ax.set_xlim(-0.5, nChannels - 0.5)

    # Annotate total in the top-right.
    total = int(perChannel.sum())
    ax.text(
        0.995, 0.95, f'total: {total:,}',
        transform=ax.transAxes, ha='right', va='top',
        fontsize=9, color='0.3',
        bbox=dict(boxstyle='round,pad=0.2',
                  facecolor='white', edgecolor='none', alpha=0.85),
    )

    fig.tight_layout()
    return fig
