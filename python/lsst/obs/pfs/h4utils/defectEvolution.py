"""Panel of small cutouts showing a defect (and its IPC halo) evolving in time.

Notebook-first. Point it at a pixel on one detector and a set of epochs
(reductions in a collection); it lays out one small cutout per epoch, sorted
by time, on a shared diverging scale so over-subtraction (blue) and
under-subtraction (red) read at a glance. The BAD mask is outlined per panel,
so you can see when a drifting defect finally trips a masking gate.

    from lsst.obs.pfs.h4utils.defectEvolution import plotDefectEvolution
    # every postISRCCD for n4 in the collection, sorted by epoch:
    fig, recs = plotDefectEvolution(
        butler, x=1999, y=428, spectrograph=4,
        collections="u/cloomis/PIPE2D-1855/ndark-drift/all4",
        half=5, datasetType="postISRCCD",
        savePath="/work/cloomis/outputs/PIPE2D-1855/ndark-drift/defect_1999_428.pdf")

``recs`` is the per-epoch table (date, dt, center value, cutout min, IPC-ring
sum) so the same call feeds both the picture and a quantitative trend.
"""
import numpy as np

__all__ = ["defectCutouts", "plotDefectEvolution"]


def _camDataId(arm, spectrograph):
    return dict(arm=arm, spectrograph=int(spectrograph))


def defectCutouts(butler, x, y, spectrograph, collections, *, arm="n",
                  datasetType="postISRCCD", half=5, visits=None,
                  maskPlanes=("BAD", "INTRP"), reference=None, instrument="PFS"):
    """Gather epoch-sorted cutouts around ``(x, y)`` for one detector.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
    x, y : `int`
        Defect pixel in image (column, row) coordinates; the cutout is
        centred here and the array is indexed ``[y, x]``.
    spectrograph : `int`
        NIR spectrograph 1-4.
    collections : collection(s) to read from.
    arm : `str`, optional
        NIR arm (default ``"n"``).
    datasetType : `str`, optional
        ``"postISRCCD"`` (default), ``"calexp"``, ...
    half : `int`, optional
        Cutout half-width; the panel is ``(2*half+1)`` square (default 5 →
        11×11, enough to hold the defect plus its IPC nearest-neighbour halo).
    visits : iterable of `int`, optional
        Restrict to these visits; otherwise every matching dataset in
        ``collections`` is used.
    maskPlanes : iterable of `str`, optional
        Planes to record as per-pixel boolean cutouts (for overlays / stats).
    reference : `str`, `numpy.ndarray`, or `None`, optional
        Subtract a common reference cutout from every panel, so the panels
        show *deviations* rather than absolute level — the way to inspect
        ``doDark=False`` products (uncorrected level, but within-set changes
        are meaningful). ``"median"`` subtracts the per-pixel median across
        the gathered stack (≈ the combined master); an ``(2*half+1)`` square
        array subtracts that directly; ``None`` (default) leaves cutouts raw.
        ``center``/``cmin``/``ipc4`` are measured on the referenced cutout.
    instrument : `str`, optional

    Returns
    -------
    records : `list` of `dict`
        One per epoch, sorted by MJD, with keys: ``visit``, ``mjd``,
        ``date``, ``dt`` (days since the earliest epoch), ``nread``,
        ``cut`` (image cutout, reference-subtracted if requested),
        ``cutRaw`` (the un-referenced cutout), ``planes`` (dict plane→bool
        cutout), ``center`` (value at the defect pixel), ``cmin`` (cutout
        minimum), ``ipc4`` (sum of the 4 edge-adjacent pixels — the IPC
        cross; diagonal coupling is negligible so the 4-neighbour sum is the
        cleanest IPC-printthrough measure).
    """
    where = (f"instrument='{instrument}' AND arm='{arm}' "
             f"AND spectrograph={int(spectrograph)}")
    if visits is not None:
        where += " AND visit IN (" + ",".join(str(int(v)) for v in visits) + ")"
    # findFirst so a chained collection (or several RUNs) yields one dataset
    # per dataId, not a duplicate per RUN.
    refs = list(butler.registry.queryDatasets(datasetType, collections=collections,
                                              where=where, findFirst=True))
    records = []
    for ref in refs:
        exp = butler.get(ref)
        img = exp.image.array
        msk = exp.mask.array
        y0, y1 = y - half, y + half + 1
        x0, x1 = x - half, x + half + 1
        cut = img[y0:y1, x0:x1].copy()
        planes = {}
        for p in maskPlanes:
            try:
                bit = exp.mask.getPlaneBitMask(p)
            except Exception:
                continue
            planes[p] = ((msk[y0:y1, x0:x1] & bit) != 0)
        md = exp.getMetadata()
        v = int(ref.dataId["visit"])
        rec = list(butler.registry.queryDimensionRecords(
            "visit", where=f"instrument='{instrument}' AND visit={v}"))[0]
        records.append(dict(
            visit=v, mjd=rec.timespan.begin.mjd,
            date=rec.timespan.begin.utc.isot[:10],
            nread=md.get("H4NREAD") or md.get("H4NTOT"),
            cutRaw=cut, planes=planes))
    records.sort(key=lambda r: r["mjd"])

    # Resolve the reference and derive the referenced cutout + stats.
    if isinstance(reference, str) and reference == "median":
        ref_cut = np.nanmedian(np.stack([r["cutRaw"] for r in records]), axis=0)
    elif reference is None:
        ref_cut = None
    else:
        ref_cut = np.asarray(reference)
    c = half
    for r in records:
        cut = r["cutRaw"] if ref_cut is None else (r["cutRaw"] - ref_cut)
        r["cut"] = cut
        r["center"] = float(cut[c, c])
        r["cmin"] = float(np.nanmin(cut))
        # IPC cross: the 4 edge-adjacent neighbours (up/down/left/right).
        r["ipc4"] = float(cut[c-1, c] + cut[c+1, c] + cut[c, c-1] + cut[c, c+1])
    if records:
        mjd0 = records[0]["mjd"]
        for r in records:
            r["dt"] = r["mjd"] - mjd0
    return records


def plotDefectEvolution(butler, x, y, spectrograph, collections, *, arm="n",
                        datasetType="postISRCCD", half=5, visits=None,
                        reference=None, ncol=5, vmax=None, cmap="RdBu_r",
                        perPanelScale=False, maskOverlay="BAD", markCenter=True,
                        title=None, savePath=None, instrument="PFS"):
    """Lay out epoch-sorted cutouts of a defect + IPC halo as a panel grid.

    Parameters
    ----------
    butler, x, y, spectrograph, collections, arm, datasetType, half, visits,
    instrument
        See `defectCutouts`.
    ncol : `int`, optional
        Panels per row (default 5).
    vmax : `float`, optional
        Symmetric colour limit (± ``vmax``). Default: the 99th percentile of
        ``|cutout|`` across all epochs (shared scale). Ignored when
        ``perPanelScale`` is set.
    cmap : `str`, optional
        Diverging colormap (default ``"RdBu_r"``: blue = negative =
        over-subtracted).
    perPanelScale : `bool`, optional
        Autoscale each panel independently instead of a shared scale. Useful
        when one epoch's amplitude dwarfs the rest.
    maskOverlay : `str` or `None`, optional
        Outline pixels carrying this mask plane in each panel (default
        ``"BAD"``); ``None`` disables the overlay.
    markCenter : `bool`, optional
        Draw a thin box around the defect pixel (default True).
    title : `str`, optional
        Figure suptitle; a sensible default is built from the coordinates.
    savePath : `str`, optional
        If given, write the figure there (PDF inferred from the extension).

    Returns
    -------
    fig : `matplotlib.figure.Figure`
    records : `list` of `dict`
        The `defectCutouts` table (already epoch-sorted).
    """
    import matplotlib.pyplot as plt

    recs = defectCutouts(butler, x, y, spectrograph, collections, arm=arm,
                         datasetType=datasetType, half=half, visits=visits,
                         maskPlanes=(maskOverlay,) if maskOverlay else (),
                         reference=reference, instrument=instrument)
    if not recs:
        raise ValueError(f"no {datasetType} found for n{spectrograph} in "
                         f"{collections}")

    if vmax is None and not perPanelScale:
        stack = np.stack([r["cut"] for r in recs])
        vmax = float(np.nanpercentile(np.abs(stack), 99)) or 1.0

    # Label the elapsed time in hours for a short (intra-acquisition) span,
    # otherwise in days.
    span = recs[-1]["dt"] if len(recs) > 1 else 0.0
    useHours = 0 < span < 2.0

    n = len(recs)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.15 * ncol, 2.35 * nrow),
                             squeeze=False)
    im = None
    for r, ax in zip(recs, axes.flat):
        if perPanelScale:
            v = float(np.nanpercentile(np.abs(r["cut"]), 99)) or 1.0
        else:
            v = vmax
        im = ax.imshow(r["cut"], origin="lower", cmap=cmap, vmin=-v, vmax=v)
        if maskOverlay and maskOverlay in r["planes"]:
            bad = r["planes"][maskOverlay]
            if bad.any():
                ax.contour(bad.astype(float), levels=[0.5], colors="lime",
                           linewidths=0.8)
        if markCenter:
            ax.add_patch(plt.Rectangle((half - 0.5, half - 0.5), 1, 1,
                                       fill=False, ec="k", lw=0.6))
        dtxt = f"{r['dt']*24:.1f}h" if useHours else f"{r['dt']:.0f}d"
        ax.set_title(f"{r['date']}  {dtxt}\ncen={r['center']:.0f} "
                     f"min={r['cmin']:.0f}", fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes.flat[n:]:
        ax.axis("off")

    if title is None:
        ov = f"; green = {maskOverlay}" if maskOverlay else ""
        refn = {"median": " (minus stack median)"}.get(
            reference, " (minus reference)" if reference is not None else "")
        title = (f"n{spectrograph} defect ({x},{y}) + IPC halo vs epoch — "
                 f"{datasetType}{refn}\nblue = over-subtracted (negative){ov}")
    fig.suptitle(title, fontsize=10)
    if im is not None and not perPanelScale:
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.5, label="e$^-$")
    else:
        fig.tight_layout(rect=(0, 0, 1, 0.96))

    if savePath:
        fig.savefig(savePath)
        print(f"wrote {savePath}")
    return fig, recs
