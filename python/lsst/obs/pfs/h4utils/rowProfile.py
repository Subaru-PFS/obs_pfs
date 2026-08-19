"""Row-background profiles: median a column band along axis=1 and plot vs row.

Notebook-first. The "middle region" (a wide inter-bundle dark strip, default
cols 2000:2100) medianed across columns gives a per-row background profile —
handy for comparing how cleanly different darks / reductions subtract.

    from lsst.obs.pfs.h4utils.rowProfile import plotMedianRowProfile
    plotMedianRowProfile(butler, [
        dict(label="new dark", visit=144959, spectrograph=1,
             collections="u/cloomis/1843-irp4_all"),
        dict(label="old dark", visit=144959, spectrograph=1,
             collections="u/cloomis/144959_olddark"),
    ])
"""
import numpy as np

__all__ = ["medianRowProfile", "plotMedianRowProfile"]


def medianRowProfile(butler, visit, spectrograph, collections, *, arm="n",
                     datasetType="postISRCCD", cols=(2000, 2100),
                     maskPlanes=("BAD",), instrument="PFS"):
    """Median of ``datasetType`` over ``cols`` (axis=1), per row.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
    visit, spectrograph : `int`
    collections : collection(s) to read from.
    arm : `str`, optional
        NIR arm (default ``"n"``).
    datasetType : `str`, optional
        ``"postISRCCD"`` (default) or ``"calexp"`` etc.
    cols : `(int, int)`, optional
        Column band ``[x0, x1)`` to median (default the (2000, 2100) dark strip).
    maskPlanes : iterable of `str`, optional
        Planes to exclude (NaN) before the median. Default ``("BAD",)``;
        pass ``()`` to median all pixels.
    instrument : `str`, optional

    Returns
    -------
    rows : `numpy.ndarray`
        Row indices ``0..H-1``.
    profile : `numpy.ndarray`
        Per-row median over the column band (NaN where a row is all-masked).
    """
    exp = butler.get(datasetType, instrument=instrument, visit=visit, arm=arm,
                     spectrograph=spectrograph, collections=collections)
    x0, x1 = cols
    img = exp.image.array[:, x0:x1].astype(float)
    if maskPlanes:
        d = exp.mask.getMaskPlaneDict()
        band = exp.mask.array[:, x0:x1]
        bad = np.zeros(band.shape, bool)
        for p in maskPlanes:
            if p in d:
                bad |= (band & (1 << d[p])) != 0
        img = np.where(bad, np.nan, img)
    with np.errstate(all="ignore"):
        profile = np.nanmedian(img, axis=1)
    return np.arange(len(profile)), profile


def plotMedianRowProfile(butler, series, *, cols=(2000, 2100),
                         datasetType="postISRCCD", arm="n", maskPlanes=("BAD",),
                         rows=None, ylim=None, ax=None, instrument="PFS"):
    """Overlay median-along-axis=1 row profiles for several series.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
    series : iterable of `dict`
        Each: ``visit``, ``spectrograph``, ``collections`` (required); optional
        ``label``, ``arm``, ``datasetType``.
    cols : `(int, int)`, optional
        Column band to median (default (2000, 2100)).
    datasetType, arm, maskPlanes : defaults applied when a series omits them.
    rows : `(int, int)`, optional
        Restrict the plotted x (row) range ``[y0, y1]`` (x-zoom).
    ylim : `(float, float)`, optional
        Fix the y (intensity) range, e.g. ``(-25, 60)`` — needed to keep a
        dead row (e.g. the LINEARITY_DEFECT row 752) from blowing up the
        autoscale when interactive zoom isn't available.
    ax : `matplotlib.axes.Axes`, optional
        Draw into this axes; a new figure is made if omitted.

    Returns
    -------
    ax : `matplotlib.axes.Axes`
    """
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))
    for s in series:
        y, prof = medianRowProfile(
            butler, s["visit"], s["spectrograph"], s["collections"],
            arm=s.get("arm", arm), datasetType=s.get("datasetType", datasetType),
            cols=cols, maskPlanes=maskPlanes, instrument=instrument)
        label = s.get("label", f"{s['visit']} n{s['spectrograph']}")
        ax.plot(y, prof, lw=0.6, label=label)
    ax.axhline(0, color="k", lw=0.3, ls=":")
    ax.set_xlabel("row (y)")
    ax.set_ylabel(f"median over cols {cols[0]}:{cols[1]} ({datasetType}, e-)")
    if rows is not None:
        ax.set_xlim(*rows)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend(fontsize=8)
    return ax
