"""Display a grid of reduced H4 exposures in ds9.

Notebook-first: open ds9, then

    from lsst.obs.pfs.h4utils.displayGrid import displayGrid
    displayGrid(butler, [138211, 144959, 144955], collections=COLL)

lays out one ds9 frame per (visit, camera) cell — visits down the rows,
cameras across the columns — with the requested mask planes overlaid
faintly and all other planes hidden. Choose ``postISRCCD`` or ``calexp``
with ``datasetType``.

Two flavours:

- `displayGrid` uses ``lsst.afw.display`` with the ``ds9`` backend. The afw
  display backends are not set up by default — first
  ``setup display_ds9`` (and ``display_matplotlib`` / ``display_astrowidgets``
  for those backends).
- `displayGridDs9` drives a `pyds9.DS9` handle you pass in, via XPA directly.
  It needs only ``pyds9`` — no ``display_*`` products.
- `displayGridCollectionsDs9` is the same as `displayGridDs9` but lays out a
  single visit with one **collection** per row (e.g. the same visit reduced
  with different calibs), cameras across the columns.
"""
import io

import lsst.afw.display as afwDisplay

__all__ = ["displayGrid", "displayGridDs9", "displayGridCollectionsDs9"]

DEFAULT_MASK_PLANES = ("BAD", "CR", "INTRP")   # INTRP is the LSST interpolated plane
DEFAULT_CAMERAS = ("n1", "n2", "n3", "n4")
DEFAULT_PLANE_COLORS = {"BAD": "red", "CR": "green", "INTRP": "blue"}


def _camDataId(cam):
    """``"n2"`` -> ``dict(arm="n", spectrograph=2)``."""
    return dict(arm=cam[0], spectrograph=int(cam[1:]))


def _tileDs9(nrows, ncols):
    """Best-effort: put ds9 into an ``ncols`` x ``nrows`` tile grid."""
    try:
        from lsst.display.ds9 import ds9Cmd
        ds9Cmd(f"tile grid layout {ncols} {nrows}", flush=True)
        ds9Cmd("tile grid mode manual", flush=True)
        ds9Cmd("tile yes", flush=True)
    except Exception as exc:  # ds9 command API varies; tiling is cosmetic
        print(f"(could not auto-tile ds9: {exc!r}; set View > Tile manually)")


def displayGrid(butler, visits, cameras=DEFAULT_CAMERAS, *,
                datasetType="calexp", collections=None,
                maskPlanes=DEFAULT_MASK_PLANES, transparency=85,
                scaleAlgorithm="asinh", scaleType="zscale",
                frame0=1, backend="ds9", instrument="PFS"):
    """Display an ``nvisit`` x ``ncam`` grid of reduced exposures in ds9.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler to read from.
    visits : iterable of `int`
        Visits; one grid row each.
    cameras : iterable of `str`, optional
        Camera names like ``"n2"`` (arm + spectrograph); one column each.
        Defaults to the four NIR detectors ``n1``..``n4``.
    datasetType : `str`, optional
        ``"calexp"`` (default) or ``"postISRCCD"``.
    collections : optional
        Collections to read from; defaults to the butler's.
    maskPlanes : iterable of `str`, optional
        Mask planes to overlay faintly. Default ``("BAD", "CR", "INTERP")``.
        All other planes are hidden.
    transparency : `float`, optional
        Transparency (0..100) for the shown planes; higher is fainter.
        Default 85.
    scaleAlgorithm, scaleType : `str`, optional
        Passed to ``Display.scale`` (default asinh / zscale).
    frame0 : `int`, optional
        ds9 frame number of the first (top-left) cell.
    backend : `str`, optional
        afw display backend (default ``"ds9"``).
    instrument : `str`, optional
        Instrument name for the dataId (default ``"PFS"``).

    Returns
    -------
    displays : `dict`
        ``{(visit, camera): lsst.afw.display.Display}`` for the cells that
        were shown (missing datasets are skipped with a message).
    """
    afwDisplay.setDefaultBackend(backend)
    visits = list(visits)
    cameras = list(cameras)
    ncam = len(cameras)

    displays = {}
    for iv, visit in enumerate(visits):
        for ic, cam in enumerate(cameras):
            frame = frame0 + iv * ncam + ic
            try:
                exp = butler.get(datasetType, instrument=instrument, visit=visit,
                                 collections=collections, **_camDataId(cam))
            except Exception as exc:
                print(f"  {visit} {cam}: no {datasetType} "
                      f"({type(exc).__name__}: {exc})")
                continue
            disp = afwDisplay.Display(frame=frame, backend=backend)
            disp.setMaskTransparency(100)                    # hide all planes
            # Only touch planes this exposure actually has: e.g. postISRCCD
            # predates interpolation, so it carries no INTRP plane.
            existing = exp.mask.getMaskPlaneDict()
            for plane in maskPlanes:
                if plane in existing:
                    disp.setMaskTransparency(transparency, plane)  # faint
            disp.mtv(exp, title=f"{visit} {cam} {datasetType}")
            disp.scale(scaleAlgorithm, scaleType)
            displays[(visit, cam)] = disp

    _tileDs9(nrows=len(visits), ncols=ncam)
    return displays


def _maskFitsBytes(boolArr):
    """FITS bytes for a boolean mask array (non-zero pixels are the mask)."""
    from astropy.io import fits as afits
    buf = io.BytesIO()
    afits.HDUList([afits.PrimaryHDU(boolArr.astype("uint8"))]).writeto(buf)
    return buf.getvalue()


def _drawCellDs9(ds9, exp, frame, *, maskPlanes, colors, transparency,
                 scaleAlgorithm, scaleMode, label=None, labelColor="green"):
    """Load one exposure into ``frame`` of ``ds9`` and overlay its mask planes.

    Shared by `displayGridDs9` and `displayGridCollectionsDs9`. When ``label``
    is given it is burned into the frame as a text region near the top-left, so
    tiled rows are unambiguous without reading a separate legend.
    """
    ds9.set(f"frame {frame}")
    ds9.set_np2arr(exp.image.array)
    ds9.set(f"scale {scaleAlgorithm}")
    ds9.set(f"scale mode {scaleMode}")
    ds9.set("mask clear")
    planes = exp.mask.getMaskPlaneDict()
    for plane in maskPlanes:
        if plane not in planes:
            continue
        bit = exp.mask.getPlaneBitMask(plane)
        setPix = (exp.mask.array & bit) != 0
        if not setPix.any():          # e.g. INTRP on postISRCCD
            continue
        ds9.set(f"mask color {colors.get(plane, 'red')}")
        ds9.set(f"mask transparency {transparency:g}")
        ds9.set("fits mask", _maskFitsBytes(setPix))
    if label:
        # Text-region point is the label centre; place near the top-left.
        # ds9's image origin is lower-left, so a large y is the top.
        h, w = exp.image.array.shape
        x, y = int(0.15 * w), int(0.95 * h)
        ds9.set("regions", f'image; text {x} {y} # text={{{label}}} '
                           f'color={labelColor} font="helvetica 14 bold"')


def displayGridDs9(ds9, butler, visits, cameras=DEFAULT_CAMERAS, *,
                   datasetType="calexp", collections=None,
                   maskPlanes=DEFAULT_MASK_PLANES, planeColors=None,
                   transparency=85, scaleAlgorithm="asinh", scaleMode="zscale",
                   frame0=1, instrument="PFS", showLabels=True,
                   labelColor="green"):
    """Draw an ``nvisit`` x ``ncam`` grid into a passed-in `pyds9.DS9`.

    Like `displayGrid`, but instead of going through ``lsst.afw.display`` it
    drives the caller's ``pyds9`` handle directly via XPA — so the caller keeps
    the handle for further control. The image is loaded with ``set_np2arr`` and
    each requested mask plane is layered on as a ds9 native mask
    (``fits mask`` + ``mask color`` / ``mask transparency``).

    Parameters
    ----------
    ds9 : `pyds9.DS9`
        The already-connected ds9 handle to drive.
    butler : `lsst.daf.butler.Butler`
        Butler to read from.
    visits : iterable of `int`
        Visits; one grid row each.
    cameras : iterable of `str`, optional
        Camera names like ``"n2"``; one column each. Default ``n1``..``n4``.
    datasetType : `str`, optional
        ``"calexp"`` (default) or ``"postISRCCD"``.
    collections : optional
        Collections to read from; defaults to the butler's.
    maskPlanes : iterable of `str`, optional
        Mask planes to overlay. Default ``("BAD", "CR", "INTRP")``. A plane the
        exposure lacks, or one with no set pixels (e.g. INTRP on ``postISRCCD``),
        is skipped.
    planeColors : `dict`, optional
        ``{plane: ds9-color}`` overrides for the overlay colors.
    transparency : `float`, optional
        ds9 mask transparency (0..100); higher is fainter. Default 85.
    scaleAlgorithm, scaleMode : `str`, optional
        ds9 ``scale`` algorithm / mode (default asinh / zscale).
    frame0 : `int`, optional
        ds9 frame number of the first (top-left) cell.
    instrument : `str`, optional
        Instrument name for the dataId (default ``"PFS"``).
    showLabels : `bool`, optional
        Burn a ``"<visit> <camera>"`` text label into each frame (default
        ``True``) so tiled cells are unambiguous.
    labelColor : `str`, optional
        ds9 color for the label text (default ``"green"``).

    Returns
    -------
    frames : `dict`
        ``{(visit, camera): ds9-frame-number}`` for the cells that were shown.
    """
    colors = dict(DEFAULT_PLANE_COLORS)
    if planeColors:
        colors.update(planeColors)
    visits = list(visits)
    cameras = list(cameras)
    ncam = len(cameras)

    ds9.set(f"tile grid layout {ncam} {len(visits)}")   # <cols> <rows>
    ds9.set("tile grid mode manual")            # honor our layout, don't auto-arrange
    ds9.set("tile yes")

    frames = {}
    for iv, visit in enumerate(visits):
        for ic, cam in enumerate(cameras):
            frame = frame0 + iv * ncam + ic
            try:
                exp = butler.get(datasetType, instrument=instrument, visit=visit,
                                 collections=collections, **_camDataId(cam))
            except Exception as exc:
                print(f"  {visit} {cam}: no {datasetType} "
                      f"({type(exc).__name__}: {exc})")
                continue
            _drawCellDs9(ds9, exp, frame, maskPlanes=maskPlanes, colors=colors,
                         transparency=transparency, scaleAlgorithm=scaleAlgorithm,
                         scaleMode=scaleMode,
                         label=f"{visit} {cam}" if showLabels else None,
                         labelColor=labelColor)
            frames[(visit, cam)] = frame

    return frames


def displayGridCollectionsDs9(ds9, butler, visit, collections,
                              cameras=DEFAULT_CAMERAS, *,
                              datasetType="calexp", maskPlanes=DEFAULT_MASK_PLANES,
                              planeColors=None, transparency=85,
                              scaleAlgorithm="asinh", scaleMode="zscale",
                              frame0=1, instrument="PFS", showLabels=True,
                              labelColor="green"):
    """Draw an ``ncollection`` x ``ncam`` grid for one visit into a `pyds9.DS9`.

    Like `displayGridDs9`, but the rows are **collections** (the same visit
    reduced different ways) rather than visits. Cameras run across the columns.

    Parameters
    ----------
    ds9 : `pyds9.DS9`
        The already-connected ds9 handle to drive.
    butler : `lsst.daf.butler.Butler`
        Butler to read from.
    visit : `int`
        The single visit shown in every cell.
    collections : iterable
        One row per entry. Each entry is either a collection (a `str` or a list
        of collections) or a ``(label, collection)`` pair; the label is used
        only for the returned dict and the printed legend.
    cameras : iterable of `str`, optional
        Camera names like ``"n2"``; one column each. Default ``n1``..``n4``.
    datasetType : `str`, optional
        ``"calexp"`` (default) or ``"postISRCCD"``.
    maskPlanes : iterable of `str`, optional
        Mask planes to overlay. Default ``("BAD", "CR", "INTRP")``. A plane the
        exposure lacks, or one with no set pixels, is skipped.
    planeColors : `dict`, optional
        ``{plane: ds9-color}`` overrides for the overlay colors.
    transparency : `float`, optional
        ds9 mask transparency (0..100); higher is fainter. Default 85.
    scaleAlgorithm, scaleMode : `str`, optional
        ds9 ``scale`` algorithm / mode (default asinh / zscale).
    frame0 : `int`, optional
        ds9 frame number of the first (top-left) cell.
    instrument : `str`, optional
        Instrument name for the dataId (default ``"PFS"``).
    showLabels : `bool`, optional
        Burn a ``"<row-label> <camera>"`` text label into each frame (default
        ``True``) so the collection rows are unambiguous.
    labelColor : `str`, optional
        ds9 color for the label text (default ``"green"``).

    Returns
    -------
    frames : `dict`
        ``{(label, camera): ds9-frame-number}`` for the cells that were shown.
    """
    colors = dict(DEFAULT_PLANE_COLORS)
    if planeColors:
        colors.update(planeColors)
    cameras = list(cameras)
    ncam = len(cameras)

    # Normalize each row to (label, collection).
    rows = []
    for entry in collections:
        if isinstance(entry, tuple) and len(entry) == 2 and not isinstance(entry[1], int):
            rows.append((str(entry[0]), entry[1]))
        else:
            rows.append((str(entry), entry))

    ds9.set(f"tile grid layout {ncam} {len(rows)}")     # <cols> <rows>
    ds9.set("tile grid mode manual")            # honor our layout, don't auto-arrange
    ds9.set("tile yes")
    print(f"visit {visit}  rows (top->bottom): "
          + ", ".join(f"{i}:{lab}" for i, (lab, _) in enumerate(rows)))

    frames = {}
    for ir, (label, coll) in enumerate(rows):
        for ic, cam in enumerate(cameras):
            frame = frame0 + ir * ncam + ic
            try:
                exp = butler.get(datasetType, instrument=instrument, visit=visit,
                                 collections=coll, **_camDataId(cam))
            except Exception as exc:
                print(f"  {label} {cam}: no {datasetType} "
                      f"({type(exc).__name__}: {exc})")
                continue
            _drawCellDs9(ds9, exp, frame, maskPlanes=maskPlanes, colors=colors,
                         transparency=transparency, scaleAlgorithm=scaleAlgorithm,
                         scaleMode=scaleMode,
                         label=f"{label} {cam}" if showLabels else None,
                         labelColor=labelColor)
            frames[(label, cam)] = frame

    return frames
