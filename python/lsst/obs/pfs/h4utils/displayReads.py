"""Display the individual reads of one H4 ramp in ds9.

Notebook-first, like `displayGrid`:

    import pyds9
    from lsst.obs.pfs.h4utils.displayReads import displayReadsDs9
    ds9 = pyds9.DS9()
    displayReadsDs9(ds9, butler, 145076, "n1", collections=COLL)

lays out one ds9 frame per read, tiled and locked, with the measured
illumination state burned into each frame's label.

Reads, not the reduced image, are what persistence questions are actually
about: which reads were lit, what the ramp did in the read after the lamp went
off, whether a feature grows or decays along the ramp. Two views, and the
difference matters:

``delta``
    what arrived *during* each read, i.e. the difference between consecutive
    reads. Use this for illumination and decay.

    The first frame has no preceding read to difference against, so it is
    shown with its own median removed. That puts it on the same zero-based
    footing as the others -- it is charge accumulated since read 0, and
    carries a pedestal the differences do not -- while preserving any
    structure in it, which for the first interval is usually persistence
    from the preceding exposure rather than the lamp.
``cumulative``
    the stored planes, each the charge since read 0. Monotonic, so it hides
    both the lamp switching off and anything decaying.

Labels come from `measureIllumination`, so they report the lamp state actually
measured from the flux rather than from ``EXPTIME``, which does not mark it --
the lamp starts about two reads into the ramp.
"""

import numpy as np

from .rampTiming import measureIllumination, lampsOn

__all__ = ["displayReadsDs9", "parseReads"]

DEFAULT_SCALE_ALGORITHM = "asinh"
DEFAULT_SCALE_MODE = "zscale"


def _camDataId(cam):
    """``"n2"`` -> ``dict(arm="n", spectrograph=2)``."""
    return dict(arm=cam[0], spectrograph=int(cam[1:]))


def parseReads(spec, nRead):
    """Turn a read selection into a list of indices.

    Parameters
    ----------
    spec : `str` or `None`
        ``"all"`` (or `None`), a comma-separated list ``"1,3,5"``, or ranges
        ``"0-4"`` / ``"2:"``. Negative indices count back from the end, so
        ``"-3:"`` is the last three reads -- the ones after the lamp.
    nRead : `int`
        Number of reads in the ramp.

    Returns
    -------
    reads : `list` of `int`
        Selected indices, in order, without duplicates.
    """
    if spec in (None, "", "all"):
        return list(range(nRead))
    out = []
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        sep = ":" if ":" in part else ("-" if "-" in part.lstrip("-") else None)
        if sep is None:
            index = int(part)
            out.append(index + nRead if index < 0 else index)
            continue
        lo, _, hi = part.partition(sep)
        lo = int(lo) if lo.strip() else 0
        hi = int(hi) if hi.strip() else nRead - 1
        if lo < 0:
            lo += nRead
        if hi < 0:
            hi += nRead
        out.extend(range(lo, min(hi, nRead - 1) + 1))
    return [k for k in dict.fromkeys(out) if 0 <= k < nRead]


def displayReadsDs9(ds9, butler, visit, camera="n1", *, reads="all",
                    mode="delta", collections=None, box=None,
                    scaleAlgorithm=DEFAULT_SCALE_ALGORITHM,
                    scaleMode=DEFAULT_SCALE_MODE, scaleLimits=None,
                    frame0=1, instrument="PFS", showLabels=True,
                    labelColor="green", tile=True, lock=True):
    """Draw one ds9 frame per read of a ramp.

    Parameters
    ----------
    ds9 : `pyds9.DS9`
        The already-connected ds9 handle to drive.
    butler : `lsst.daf.butler.Butler`
        Butler to read from.
    visit : `int`
        Visit to display.
    camera : `str`, optional
        Camera name like ``"n1"``. Default ``"n1"``.
    reads : `str`, optional
        Read selection; see `parseReads`. Default all of them.
    mode : `str`, optional
        ``"delta"`` (default) or ``"cumulative"``.
    collections : optional
        Collections to read from; defaults to the butler's.
    box : `tuple` of `slice`, optional
        Sub-region to display, and to measure the illumination over. A whole
        4096x4096 frame per read is a lot of ds9 frames; a box is usually
        what you want.
    scaleAlgorithm, scaleMode : `str`, optional
        ds9 ``scale`` algorithm / mode (default asinh / zscale).
    scaleLimits : `tuple`, optional
        ``(lo, hi)`` manual limits, overriding ``scaleMode``. Worth setting
        when comparing reads, since per-frame zscale hides the very changes
        along the ramp that are being looked for.
    frame0 : `int`, optional
        ds9 frame number of the first read shown.
    instrument : `str`, optional
        Instrument name for the dataId (default ``"PFS"``).
    showLabels : `bool`, optional
        Burn ``"<read> <time> <lit fraction>"`` into each frame.
    labelColor : `str`, optional
        Colour of that label.
    tile, lock : `bool`, optional
        Tile the frames, and lock frame/scale/colorbar together.

    Returns
    -------
    frames : `dict`
        ``{read index: ds9 frame number}``.
    """
    if mode not in ("delta", "cumulative"):
        raise ValueError(f"mode must be 'delta' or 'cumulative', got {mode!r}")

    dataId = dict(instrument=instrument, visit=visit, **_camDataId(camera))
    kwargs = dict(collections=collections) if collections is not None else {}
    cube = butler.get("rawISRCube", dataId, **kwargs)
    nRead = cube.getNumReads()
    window = measureIllumination(cube, box=box)
    lamps = "+".join(sorted(lampsOn(cube.metadata))) or "none"
    wanted = parseReads(reads, nRead)
    if not wanted:
        raise ValueError(f"no reads selected by {reads!r}; ramp has {nRead}")

    frames = {}
    previous = None
    frame = frame0
    for index in range(nRead):
        array = np.asarray(cube.getReadArray(index), dtype=np.float32)
        if mode == "delta":
            if previous is None:
                # Nothing to difference against: show it relative to its own
                # median so its pedestal does not swamp the display, and it
                # shares a zero point with the differences that follow.
                shown = array[box] if box is not None else array
                current = array - float(np.median(shown))
            else:
                current = array - previous
            previous = array
        else:
            current, previous = array, array
        if index not in wanted:
            continue

        ds9.set(f"frame {frame}")
        ds9.set_np2arr(current[box] if box is not None else current)
        ds9.set(f"scale {scaleAlgorithm}")
        if scaleLimits is not None:
            lo, hi = scaleLimits
            ds9.set(f"scale limits {lo:g} {hi:g}")
        else:
            ds9.set(f"scale mode {scaleMode}")
        if showLabels:
            # Name the interval, not just its lower index. Plane i is
            # read[i+1] - read[0], so what is displayed is the charge that
            # arrived BETWEEN reads i and i+1 -- labelling it "read i" invites
            # the reading that it is a single read, and the first interval in
            # particular is routinely mistaken for the lamp when what it holds
            # is persistence from the previous exposure.
            if window is None:
                label = f"{visit} {camera} reads {index}-{index + 1}  ({lamps})"
            else:
                end = window.ends[index]
                lit = 100*window.litFraction(index)
                # The first interval is excluded from lamp detection, so any
                # flux in it is persistence released by the preceding
                # exposures. It follows the same traces as the illumination,
                # which is what makes it so easy to misread as the lamp.
                note = ("  PERSISTENCE, not lamp"
                        if index == 0 and lit > 2 else "")
                label = (f"{visit} {camera} reads {index}-{index + 1}  "
                         f"{end - window.frameTime:.0f}-{end:.0f}s  "
                         f"{lit:.0f}% lit  ({lamps}){note}")
            shape = (current[box] if box is not None else current).shape
            x, y = int(0.30*shape[1]), int(0.95*shape[0])
            ds9.set("regions", f'image; text {x} {y} # text={{{label}}} '
                               f'color={labelColor} font="helvetica 14 bold"')
        frames[index] = frame
        frame += 1
    del cube

    if tile:
        ds9.set("tile yes")
    if lock:
        ds9.set("lock frame image")
        ds9.set("lock scale yes")
        ds9.set("lock colorbar yes")
    return frames
