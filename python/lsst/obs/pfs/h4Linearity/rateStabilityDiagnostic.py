"""Tuning diagnostic for the rate-stability rejection.

Runs the linearized ISR ramp path on a real exposure, applies
:func:`rateStability.detectRateInstability`, and plots the resulting
fractional-disagreement distribution and per-pixel deltas so the
``rateStabilityThreshold`` can be checked on real ramps.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np

from . import rateStability
from .validate import processRamp, _resolveCam

__all__ = [
    "RateStabilityDiagnosticData",
    "runRateStabilityDiagnostic",
    "plotRateStabilityFraction",
    "plotUnstableLocations",
    "plotPixelRampGrid",
    "makeUnstablePixelReport",
]


@dataclass
class RateStabilityDiagnosticData:
    """Bundle returned by :func:`runRateStabilityDiagnostic`.

    Carries everything the rate-stability plots need: the
    rate-stability result, the linearized CR-corrected per-read deltas
    and the matching cumulative cube (for per-pixel delta and flux
    panels), the per-delta CR/glitch flag mask (for marker styling),
    the 2-D ISR image (location-map background), and the run metadata.

    Attributes
    ----------
    result : rateStability.RateStabilityResult
        The rate-stability test output.
    deltas : np.ndarray
        ``(H, W, N-1)`` float32 linearized, CR-corrected per-read deltas.
    cube : np.ndarray
        ``(H, W, N)`` float32 linearized, CR-corrected cumulative
        ramp. Same data as ``deltas`` but in cumulative form for the
        flux panels in :func:`plotPixelRampGrid`.
    flagMask : np.ndarray
        ``(H, W, N-1)`` bool; True at deltas flagged as CR or glitch
        (for visually distinct markers in the delta panels).
    image : np.ndarray
        ``(H, W)`` float32 2-D ISR science image (location-map
        background).
    cam : str
        Detector name (e.g. ``"n3"``).
    visit : int
        Exposure visit number.
    threshold : float
        Rejection threshold for the fractional disagreement.
    rateFloorADU : float
        Lower bound on the denominator of the fractional metric.
    minDeltasPerSegment : int
        Minimum un-flagged deltas for a half to be testable.
    """

    result: "rateStability.RateStabilityResult"
    deltas: np.ndarray
    cube: np.ndarray
    flagMask: np.ndarray
    image: np.ndarray
    cam: str
    visit: int
    threshold: float
    rateFloorADU: float
    minDeltasPerSegment: int


def _paginatePixels(coords, *, maxPixels, panelsPerPage):
    """Cap a list of pixel coordinates and split it into per-page chunks.

    Parameters
    ----------
    coords : sequence
        Pixel coordinates (e.g. ``(y, x)`` tuples), in display order.
    maxPixels : int
        Keep at most this many coordinates; the rest are dropped.
    panelsPerPage : int
        Chunk size for pagination (>= 1).

    Returns
    -------
    pages : list of list
        ``coords`` capped to ``maxPixels`` then split into chunks of
        ``panelsPerPage`` (the last chunk may be shorter).
    nShown : int
        Number of coordinates kept (``min(len(coords), maxPixels)``).
    nTotal : int
        ``len(coords)`` -- the full count before capping.
    """
    if maxPixels < 0 or panelsPerPage < 1:
        raise ValueError("maxPixels must be >= 0 and panelsPerPage >= 1")
    coords = list(coords)
    nTotal = len(coords)
    shown = coords[:maxPixels]
    nShown = len(shown)
    pages = [shown[i:i + panelsPerPage]
             for i in range(0, nShown, panelsPerPage)]
    return pages, nShown, nTotal


def _rejectedCoordsByFraction(rejectMask, fraction):
    """Rejected-pixel ``(y, x)`` coordinates, ordered by fraction descending.

    The fraction is the rate-stability disagreement metric, so this
    orders the pixels most-unstable first.

    Parameters
    ----------
    rejectMask : np.ndarray
        ``(H, W)`` bool; True at rejected pixels.
    fraction : np.ndarray
        ``(H, W)`` float; the per-pixel fractional disagreement.

    Returns
    -------
    list of (int, int)
        ``(y, x)`` tuples, largest fraction first.
    """
    ys, xs = np.where(rejectMask)
    coords = list(zip(ys.tolist(), xs.tolist()))
    coords.sort(key=lambda yx: float(fraction[yx[0], yx[1]]), reverse=True)
    return coords


def runRateStabilityDiagnostic(
    butler,
    dataId,
    *,
    cam: Optional[str] = None,
    threshold: float = 0.20,
    rateFloorADU: float = 5.0,
    minDeltasPerSegment: int = 3,
    firstRead: Optional[int] = None,
    lastRead: Optional[int] = None,
    nirDark=None,
    defects=None,
    linearity=None,
    log=None,
) -> RateStabilityDiagnosticData:
    """Run the ISR ramp path and the rate-stability test on one exposure.

    The pipeline runs at its configured rate-stability defaults and
    captures the result; this function then *re-runs* the metric at
    the threshold / rate floor / segment-size knobs requested in this
    call, against the same pixel set the pipeline used (echoed back as
    ``RateStabilityResult.goodPixelMask``). This lets a single ramp
    sweep multiple knob settings without re-doing linearization or CR
    repair.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
    dataId : mapping
        Must resolve ``visit`` (and ``arm``/``spectrograph`` as needed).
    cam : str, optional
        Detector name (e.g. ``"n3"``); inferred from the raw if omitted.
    threshold, rateFloorADU, minDeltasPerSegment
        Forwarded to :func:`rateStability.detectRateInstability` for
        the diagnostic re-run; the pipeline's production defaults
        determined which pixels were considered good in the first
        place.
    firstRead, lastRead : int, optional
        Read range; forwarded to ``processRamp``.
    nirDark, defects, linearity, log : optional
        Forwarded to ``processRamp``.

    Returns
    -------
    RateStabilityDiagnosticData
        The result, the linearized CR-corrected deltas, the per-delta
        flag mask, the 2-D ISR image, and the run metadata -- ready
        for the plotting functions.
    """
    cam, _ = _resolveCam(cam, None, butler, dataId)
    intermediates = dict.fromkeys(
        ['crCorrected', 'crResult', 'rateStabilityResult'])
    exposure, cube, crResult = processRamp(
        butler, dataId, cam=cam,
        doLinearize=True, doCR=True, repairCR=True,
        firstRead=firstRead, lastRead=lastRead,
        nirDark=nirDark, defects=defects, linearity=linearity, log=log,
        intermediates=intermediates,
    )
    rsProd = intermediates.get('rateStabilityResult')
    if rsProd is None:
        raise RuntimeError(
            "makeNirExposure did not populate intermediates['rateStabilityResult']; "
            "config.h4.doRateStability must be True and the linearized UTR path "
            "must run for the rate-stability diagnostic."
        )

    # cube is (H, W, N) post-CR cumulative. Build deltas in the same
    # (H, W, N-1) layout the production rate-stability test consumed.
    # Keep the cube around for the flux panels.
    cube = np.asarray(cube, dtype=np.float32)
    deltas = np.diff(cube, axis=-1)
    if crResult is not None:
        flagMask = crResult.crFlagMask | crResult.glitchFlagMask
    else:
        flagMask = np.zeros(deltas.shape, dtype=bool)

    result = rateStability.detectRateInstability(
        deltas, flagMask,
        goodPixelMask=rsProd.goodPixelMask,
        threshold=threshold,
        rateFloorADU=rateFloorADU,
        minDeltasPerSegment=minDeltasPerSegment,
    )
    return RateStabilityDiagnosticData(
        result=result,
        deltas=deltas,
        cube=cube,
        flagMask=flagMask,
        image=np.asarray(exposure.image.array, dtype=np.float32),
        cam=str(cam),
        visit=int(dataId["visit"]),
        threshold=threshold,
        rateFloorADU=rateFloorADU,
        minDeltasPerSegment=minDeltasPerSegment,
    )


def plotRateStabilityFraction(data, *, ax=None):
    """Histogram the rate-stability fraction with the rejection threshold.

    Tested pixels (those with finite ``fraction``) contribute to the
    histogram; the red dashed line marks ``threshold`` -- pixels above
    it are rejected.

    Parameters
    ----------
    data : RateStabilityDiagnosticData
        Output of :func:`runRateStabilityDiagnostic`.
    ax : matplotlib axis, optional
        Drawn into if supplied; otherwise a new figure is made.

    Returns
    -------
    matplotlib axis
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _fig, ax = plt.subplots(figsize=(7, 4))

    result = data.result
    thr = float(data.threshold)
    vals = result.fraction[np.isfinite(result.fraction)]

    upper = float(np.percentile(vals, 99.9)) if vals.size else thr
    hi = max(thr * 1.5, upper, 1.0)
    bins = np.linspace(0.0, hi, 100)
    if vals.size:
        ax.hist(vals, bins=bins, density=True, histtype="step",
                label=f"fraction (n={vals.size})")
        ax.set_yscale("log")
    else:
        ax.text(0.5, 0.5, "no testable pixels", ha="center", va="center",
                transform=ax.transAxes, color="gray")
    ax.axvline(thr, color="red", linestyle="--",
               label=f"reject thr = {thr:.0%}")
    ax.set_xlabel("rate-stability fraction")
    ax.set_ylabel("density")
    ax.set_title(f"rate-stability fraction  nRejected={result.nRejected}  "
                 f"nUntestable={result.nUntestable}")
    ax.legend()
    return ax


def plotUnstableLocations(data, *, ax=None):
    """Plot the detector image with the rate-stability-rejected pixels marked.

    Parameters
    ----------
    data : RateStabilityDiagnosticData
        Output of :func:`runRateStabilityDiagnostic`.
    ax : matplotlib axis, optional
        Drawn into if supplied; otherwise a new figure is made.

    Returns
    -------
    matplotlib axis
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _fig, ax = plt.subplots(figsize=(7, 7))

    img = data.image
    lo, hi = np.nanpercentile(img, [1.0, 99.0])
    ax.imshow(img, origin="lower", cmap="gray", vmin=lo, vmax=hi,
              interpolation="nearest")
    ys, xs = np.where(data.result.rejectMask)
    ax.scatter(xs, ys, s=3, c="red", alpha=0.35,
               edgecolors="none", linewidths=0,
               label=f"RATE_UNSTABLE ({len(ys)})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(loc="upper right")
    ax.set_title(
        f"RATE_UNSTABLE locations  {data.cam} visit={data.visit}  "
        f"N={data.result.nRejected}  "
        f"(threshold={data.threshold:.0%}, "
        f"rateFloor={data.rateFloorADU:.1f})"
    )
    return ax


def plotPixelRampGrid(data, coords, *, ncols=4, suptitle=None):
    """Plot a grid of per-pixel delta + flux panels.

    Layout: ``ncols`` pixel-cells per row, each cell composed of two
    sub-axes side-by-side — a delta panel on the left and a flux panel
    on the right. The standard report call passes 32 pixel coords at
    ``ncols=4`` to fill an 8×4 grid.

    Every rate plotted on either panel is the production UTR-weighted
    rate from :meth:`PfsIsrTask.calcUTRrateFromDeltas` applied to the
    relevant delta span: the full ramp (dashed grey reference), the
    first half (colour C0), and the second half (colour C1). CR/glitch-
    flagged deltas are drawn as grey crosses so the eye can tell which
    samples the test excluded.

    Parameters
    ----------
    data : RateStabilityDiagnosticData
        Output of :func:`runRateStabilityDiagnostic`.
    coords : sequence of (y, x)
        Pixel coordinates to plot, one cell each.
    ncols : int
        Pixel-cells per row. Each pixel cell consumes 2 sub-axis
        columns, so the matplotlib grid is ``(nRows, 2 * ncols)``.
    suptitle : str, optional
        Figure title.

    Returns
    -------
    matplotlib figure
    """
    import matplotlib.pyplot as plt
    from lsst.obs.pfs.isrTask import PfsIsrTask

    coords = list(coords)
    if not coords:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.axis("off")
        ax.text(0.5, 0.5, "no pixels in bin", ha="center", va="center",
                fontsize=14, color="red", transform=ax.transAxes)
        if suptitle:
            fig.suptitle(suptitle, fontsize=10)
        return fig

    nRows = (len(coords) + ncols - 1) // ncols
    # 2 sub-axes per pixel cell; figure size scales with the grid.
    fig, axes = plt.subplots(
        nRows, 2 * ncols,
        figsize=(2.0 * 2 * ncols, 1.7 * nRows),
        squeeze=False,
    )

    deltas = data.deltas
    cube = data.cube
    flagMask = data.flagMask
    result = data.result
    nDeltas = deltas.shape[-1]
    nReads = cube.shape[-1]
    deltaIdx = np.arange(nDeltas)
    readIdx = np.arange(nReads)
    bounds = np.linspace(0, nDeltas, 3).astype(int)
    midDelta = int(bounds[1])
    # First half spans cumulative reads [0, midDelta]; second half
    # spans [midDelta, nReads-1]. Cumulative read k corresponds to the
    # state after delta k-1, so the half boundary on the flux axis is
    # at read midDelta.

    for cellIdx in range(nRows * ncols):
        rowIdx = cellIdx // ncols
        colIdx = cellIdx % ncols
        axD = axes[rowIdx, 2 * colIdx]      # delta panel
        axF = axes[rowIdx, 2 * colIdx + 1]  # flux panel

        if cellIdx >= len(coords):
            axD.axis("off")
            axF.axis("off")
            continue

        y, x = coords[cellIdx]
        d = deltas[y, x, :]
        c = cube[y, x, :]
        f = flagMask[y, x, :]

        # Production rate (UTR-weighted) on each span.
        fullRate = float(PfsIsrTask.calcUTRrateFromDeltas(d))
        h1Rate = float(PfsIsrTask.calcUTRrateFromDeltas(d[:midDelta]))
        h2Rate = float(PfsIsrTask.calcUTRrateFromDeltas(d[midDelta:]))

        # --- Delta panel -------------------------------------------
        axD.plot(deltaIdx[~f], d[~f], "k.", ms=3)
        if f.any():
            axD.plot(deltaIdx[f], d[f], "x", color="0.6", ms=4,
                     markeredgewidth=0.8)
        axD.hlines(h1Rate, bounds[0] - 0.5, bounds[1] - 0.5,
                   colors="C0", lw=1.5)
        axD.hlines(h2Rate, bounds[1] - 0.5, bounds[2] - 0.5,
                   colors="C1", lw=1.5)
        axD.axhline(fullRate, color="0.4", ls="--", lw=0.7, alpha=0.7)
        axD.axvline(bounds[1] - 0.5, color="0.85", lw=0.5)
        axD.set_title(f"(x={x}, y={y})", fontsize=7)
        axD.tick_params(labelsize=6)

        # --- Flux panel --------------------------------------------
        axF.plot(readIdx, c, "k.", ms=3)
        # Full-ramp UTR rate as a line anchored on the first cumulative
        # value. Slope is the production rate; this is what the
        # pipeline's nirImage / nReads tracks visually.
        axF.plot(readIdx, c[0] + fullRate * (readIdx - readIdx[0]),
                 color="0.4", ls="--", lw=0.7, alpha=0.7)
        # Per-half UTR rates as line segments anchored at each half's
        # own first cumulative value, so a vertical step at the
        # boundary on a broken pixel shows up immediately.
        axF.plot(readIdx[: midDelta + 1],
                 c[0] + h1Rate * (readIdx[: midDelta + 1] - readIdx[0]),
                 color="C0", lw=1.5)
        axF.plot(readIdx[midDelta:],
                 c[midDelta] + h2Rate * (readIdx[midDelta:] - readIdx[midDelta]),
                 color="C1", lw=1.5)
        axF.axvline(midDelta - 0.5, color="0.85", lw=0.5)
        axF.tick_params(labelsize=6)

        # Compact stats box on the delta panel.
        frac = float(result.fraction[y, x])
        statlines = [
            f"r1={h1Rate:7.2f}",
            f"r2={h2Rate:7.2f}",
            f"r={fullRate:7.2f}",
            f"frac={frac:5.1%}",
        ]
        axD.text(0.03, 0.97, "\n".join(statlines), transform=axD.transAxes,
                 va="top", ha="left", fontsize=5.5, family="monospace",
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.85))

    fig.supxlabel("delta idx  (left in each pair)   |   read idx  (right)",
                  fontsize=8)
    fig.supylabel("ADU / read   (delta)   |   cumulative ADU   (flux)",
                  fontsize=8)
    if suptitle:
        fig.suptitle(suptitle, fontsize=10)
    fig.tight_layout()
    return fig


def _rejectsByRateBand(data, *, lo, hi):
    """Reject coords with production rate in ``[lo, hi)``, fraction-descending.

    "Production rate" here means
    :meth:`PfsIsrTask.calcUTRrateFromDeltas` on the pixel's full delta
    ramp — the same number the pipeline writes to ``exposure.image``
    divided by ``nDeltas``.
    """
    from lsst.obs.pfs.isrTask import PfsIsrTask

    rejectMask = data.result.rejectMask
    fraction = data.result.fraction
    ys, xs = np.where(rejectMask)
    if ys.size == 0:
        return []
    # Production rate per rejected pixel.
    rate = PfsIsrTask.calcUTRrateFromDeltas(data.deltas[ys, xs, :])
    rateAbs = np.abs(rate)
    inBand = (rateAbs >= lo) & (rateAbs < hi)
    if not inBand.any():
        return []
    ys = ys[inBand]
    xs = xs[inBand]
    frac = fraction[ys, xs]
    order = np.argsort(-frac)            # fraction descending
    return [(int(ys[i]), int(xs[i])) for i in order]


def makeUnstablePixelReport(
    data, outputPath, *,
    faintMaxRate=20.0,
    brightMinRate=20.0,
    nFaintPages=4,
    nBrightPages=4,
    panelsPerPage=32,
    ncols=4,
):
    """Write a multi-page PDF report of the rate-stability rejects.

    Page 1 carries the location map (rejects as small translucent
    markers) and the rate-stability fraction histogram. The remaining
    pages are split into two banded sections, each paginated through
    :func:`plotPixelRampGrid` (delta + flux panels per pixel):

      - Fainter band: ``rateFloorADU <= |rate| < faintMaxRate`` —
        the population that's above the noise floor but well short
        of "bright."
      - Brighter band: ``|rate| >= brightMinRate``.

    The rate used for binning is the production
    :meth:`PfsIsrTask.calcUTRrateFromDeltas` applied per pixel — the
    same number the pipeline writes to ``exposure.image / nDeltas``.
    Both bands sort pixels by rate-stability fraction descending so
    the most-unstable members come first.

    Parameters
    ----------
    data : RateStabilityDiagnosticData
        Output of :func:`runRateStabilityDiagnostic`.
    outputPath : str
        Destination PDF path.
    faintMaxRate : float
        Upper edge of the faint band (ADU/read).
    brightMinRate : float
        Lower edge of the bright band (ADU/read).
    nFaintPages, nBrightPages : int
        Cap on pages drawn from each band.
    panelsPerPage : int
        Pixel cells per page. Default 32 → 8 rows × 4 ``ncols`` cells.
    ncols : int
        Pixel cells per row.

    Returns
    -------
    str
        ``outputPath``.
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    rateFloor = float(data.rateFloorADU)
    faintCoords = _rejectsByRateBand(data, lo=rateFloor, hi=float(faintMaxRate))
    brightCoords = _rejectsByRateBand(data, lo=float(brightMinRate), hi=np.inf)

    def _pagesFor(coords, nPages):
        # Cap so we draw at most ``nPages`` pages of ``panelsPerPage``.
        cap = panelsPerPage * nPages
        pages, nShown, nTotal = _paginatePixels(
            coords, maxPixels=cap, panelsPerPage=panelsPerPage)
        return pages, nShown, nTotal

    faintPages, faintShown, faintTotal = _pagesFor(faintCoords, nFaintPages)
    brightPages, brightShown, brightTotal = _pagesFor(brightCoords, nBrightPages)

    with PdfPages(outputPath) as pdf:
        # Cover page: location map + fraction histogram.
        fig, (axTop, axBot) = plt.subplots(2, 1, figsize=(8, 11))
        plotUnstableLocations(data, ax=axTop)
        plotRateStabilityFraction(data, ax=axBot)
        if data.result.nRejected == 0:
            axTop.text(0.5, 0.5, "no RATE_UNSTABLE pixels",
                       transform=axTop.transAxes, ha="center", va="center",
                       fontsize=14, color="red")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        def _emitBand(label, pages, nShown, nTotal, loRate, hiRate):
            for pageno, page in enumerate(pages, start=1):
                suptitle = (
                    f"{label} rejects  {data.cam} visit={data.visit}  "
                    f"rate ∈ [{loRate:g}, {hiRate:g}) ADU/rd  "
                    f"page {pageno}/{len(pages)}"
                )
                if nShown < nTotal:
                    suptitle += f"  [showing {nShown} of {nTotal}]"
                fig = plotPixelRampGrid(data, page, ncols=ncols,
                                        suptitle=suptitle)
                pdf.savefig(fig)
                plt.close(fig)

        _emitBand("FAINTER", faintPages, faintShown, faintTotal,
                  rateFloor, float(faintMaxRate))
        _emitBand("BRIGHTER", brightPages, brightShown, brightTotal,
                  float(brightMinRate), float("inf"))

    return outputPath
