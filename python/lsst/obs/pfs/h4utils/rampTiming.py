"""Measure when a ramp was actually illuminated, from its flux.

The headers do not say. ``DARKTIME`` is ``nreads * frameTime`` by definition
and ``EXPTIME`` is the *commanded* lamp time, while ``MJD-END`` is merely
``MJD-STR + EXPTIME``. None of them records when the lamp switched, and the
commanding turns it on roughly two reads into the ramp, so on a 30 s quartz
flat ``MJD-END`` sits about 24 s before the lamp actually went off. Anything
that needs an illumination epoch -- persistence decay above all -- has to
measure it.

The flux does record it. A bright flat's ramp is piecewise linear: dark before
the lamp, rising at a constant rate while it is on, dark after. The reads that
straddle each transition are partially illuminated, and the fraction of the
plateau rate they carry places that transition to a fraction of a read.

**A read is not an instantaneous sample.** A row pointer sweeps continuously
down the detector, from the top row to the bottom, then wraps to the top; one
pass is one read. The images as stored are rotated relative to that, so a
detector row maps to an image *column* and time runs from the left of the
image to the right. A lamp switching mid-read therefore leaves a *spatial*
boundary rather than a uniform partial level: the columns scanned before the
switch show nothing and those after show light, with the ratio to a fully
illuminated read rising linearly with column. `measureEdgeByScan` fits that
ramp and places the transition to a fraction of a read, which is far better
than anything a whole-read median can do -- a median averages the lit and
unlit halves together and can only place an edge to the nearest read.

Do not mistake that spatial gradient for a dose-dependence. Comparing bright
against faint pixels across the whole frame samples different scan positions,
not different doses, and produces a difference that looks like persistence and
is not: within a fixed column band, bright and faint pixels agree to a few
percent.

The cameras run independent frame clocks, so their ramps begin up to about one
frame apart. Do **not** try to remove that with ``MJD-STR``: measured against
the lamp, which is common to all four, applying those header offsets makes the
cameras agree worse, roughly doubling their spread. Work in each camera's own
ramp coordinates.
"""

from dataclasses import dataclass, field

import numpy as np

from ..utils import getLamps

__all__ = ["IlluminationWindow", "measureIllumination", "lampsOn",
           "TRUSS_LAMP_KEYS"]

#: Telescope truss lamps. Unlike the AIT/IIS lamps, which `getLamps` reports
#: from boolean cards, these carry commanded and actual **voltages**, and an
#: exposure using them is timed by the shutter rather than by switching the
#: lamp. A ramp lit by these has every AIT/IIS boolean False, so `getLamps`
#: alone returns an empty set for it.
TRUSS_LAMP_KEYS = tuple(f"W_TFF{i}{suffix}" for i in (1, 2, 3, 4)
                        for suffix in ("VV", "VC"))

#: A read carrying more than this fraction of the plateau constrains an edge.
EDGE_FRACTION = 0.02
#: A partially illuminated read must show the same fraction across bright and
#: faint pixels, within this tolerance, to count as illumination rather than
#: persistence. See `_isIllumination`.
BRIGHTNESS_TOLERANCE = 0.15
#: A read above this fraction of the plateau counts as fully illuminated.
LIT_FRACTION = 0.5


def lampsOn(md):
    """Lamps this exposure used, including the telescope truss lamps.

    `lsst.obs.pfs.utils.getLamps` covers the AIT and IIS lamps across several
    header generations and falls back on ``DATA-TYP``; it does not know about
    the truss lamps, which are on separate voltage cards.

    Parameters
    ----------
    md : `dict`
        Header metadata.

    Returns
    -------
    lamps : `set` of `str`
        Lit lamps, truss entries named ``"truss<N>"``.
    """
    lamps = set(getLamps(md))
    for i in (1, 2, 3, 4):
        volts = md.get(f"W_TFF{i}VV", None) or md.get(f"W_TFF{i}VC", None)
        if volts:
            lamps.add(f"truss{i}")
    return lamps


@dataclass
class IlluminationWindow:
    """When a ramp was illuminated, measured from its flux.

    All times are seconds from the ramp's read 0. Read ``k`` completes at
    ``(k + 1) * frameTime``, and cube plane ``i`` is ``read[i+1] - read[0]``,
    so plane ``i`` ends at ``(i + H4READ0 + 2) * frameTime``.

    Attributes
    ----------
    on : `float`
        Lamp-on, placed within a read by that read's illuminated fraction.
    off : `float`
        Lamp-off, as ``on + EXPTIME``. EXPTIME is the commanded lamp time and
        the leading edge is measured, so this is more reliable than a
        flux-derived trailing edge.
    offFlux, offUpper : `float`
        The flux-derived trailing edge and its upper bound, for comparison.
        Both are biased: persistence in the straddling read pulls them late,
        saturation of the brightest pixels pulls them early.
    plateau : `float`
        Illuminated rate, e-/s.
    rate : `np.ndarray`
        Per-read rate, e-/s.
    ends : `np.ndarray`
        Per-read end time, seconds from read 0.
    leadingDark, trailingDark : `int`
        Whole reads entirely before lamp-on / after lamp-off. The commanding
        intends at least one of each, preferably two trailing, so that the
        concurrent b/r CCD readouts get their full commanded illumination.
    lampStillOn : `bool`
        True if the ramp ended before the lamp did, so there is no trailing
        unilluminated read at all.
    """

    on: float
    off: float
    offFlux: float
    offUpper: float
    plateau: float
    frameTime: float
    nRead: int
    rate: np.ndarray = field(repr=False)
    ends: np.ndarray = field(repr=False)
    leadingDark: int
    trailingDark: int
    lampStillOn: bool
    exptime: float
    darktime: float
    mjdStart: float

    @property
    def duration(self):
        """Measured illumination duration, seconds."""
        return self.off - self.on

    @property
    def offMjd(self):
        """Absolute MJD of lamp-off, on this camera's own ramp clock."""
        return self.mjdStart + self.off/86400.0

    def litFraction(self, read):
        """Fraction of ``read`` that was illuminated, from its rate."""
        return float(np.clip(self.rate[read]/self.plateau, 0.0, 1.0))


def _isIllumination(delta, reference, plateauBright, plateauFaint,
                    tolerance=BRIGHTNESS_TOLERANCE):
    """Is a partly-filled read illumination, or persistence?

    **This test is unreliable and is kept only as a diagnostic.** Its premise
    -- that a lamp leaves every pixel with the same fraction of its
    illuminated rate -- is false, because the readout scans across the frame
    time: a pixel's illuminated fraction depends on where it sits in the scan.
    Comparing bright against faint pixels drawn from the whole frame therefore
    compares scan positions as much as brightnesses, and on 145076 that
    produced an apparent 45%-against-1% split which looks like a
    dose-dependence and is not. Within a fixed column band the same comparison
    gives agreement to a few percent.

    Use `measureEdgeByScan` to place an edge. This function can only say
    whether a read's illuminated fraction varies with brightness *at fixed
    scan position*, which the caller must arrange.

    Parameters
    ----------
    delta : `np.ndarray`
        The read in question, as a rate.
    reference : `np.ndarray`
        A read known to be fully illuminated, as a rate.
    plateauBright, plateauFaint : `float`
        The illuminated rate over the bright and faint pixel samples.

    Returns
    -------
    isIllumination : `bool`
    """
    # Both samples must be genuinely illuminated, or the comparison is
    # against noise: on a dithered flat the median pixel sees nothing, so a
    # "faint" band taken from the whole distribution has no plateau to take a
    # fraction of. Two bands from within the lit pixels, differing severalfold
    # in brightness, is what the test needs.
    bright = reference >= np.percentile(reference, 99.0)
    faint = ((reference >= np.percentile(reference, 90.0))
             & (reference < np.percentile(reference, 99.0)))
    if not bright.any() or not faint.any() or plateauFaint <= 0:
        return True
    fracBright = float(np.median(delta[bright]))/plateauBright
    fracFaint = float(np.median(delta[faint]))/plateauFaint
    # Signed, because the two contaminants push opposite ways. Persistence
    # scales with prior dose, so it lifts the BRIGHT pixels above the faint
    # ones. Saturation on the brightest pixels does the reverse, suppressing
    # them late in a bright ramp -- that is still illumination, so only an
    # excess on the bright side disqualifies a read.
    return (fracBright - fracFaint) <= tolerance


def measureIllumination(cube, box=None, minPlateau=5.0,
                        minPlateauReads=2, litPercentile=99.0,
                        firstUsableRead=0):
    """Measure a ramp's illumination window from its flux.

    Parameters
    ----------
    cube : `lsst.obs.pfs.ImageCube`
        The ramp, as read from a ``rawISRCube``.
    box : `tuple` of `slice`, optional
        Region to measure over. Defaults to the whole detector, which is
        usually wasteful; a region on illuminated fibers is enough and much
        faster.
    litPercentile : `float`, optional
        Percentile of pixels used to define the illuminated level. A fiber
        flat lights a few percent of the detector, so a **median** over any
        region is dominated by the gaps between traces and reports no
        illumination at all while the traces are at hundreds of e-/s. The
        default follows the brightest few percent instead.
    minPlateau : `float`, optional
        Reject ramps whose peak rate is below this (e-/s) as too faint to
        locate an edge. Arcs frequently are, and so is any flat measured over
        a region that happens to lie between fibers -- measure over a region
        that is actually illuminated.
    minPlateauReads : `int`, optional
        Require the plateau to be held across at least this many reads.
        Guards against a single anomalous read -- typically read 0, carrying
        reset behaviour -- being mistaken for the illuminated level.
    firstUsableRead : `int`, optional
        Ignore reads before this when locating the plateau. Defaults to 0:
        the lamp does arrive in the first interval on these ramps -- measured
        by the scan, at 14.4-18.2 s -- so skipping it forces the answer onto
        a read boundary. Raise it only if a ramp's first read is known to be
        corrupted.

    Returns
    -------
    window : `IlluminationWindow` or `None`
        `None` if no usable plateau was found.
    """
    meta = cube.metadata
    nRead = cube.getNumReads()
    frameTime = float(meta["W_H4FRMT"])
    read0 = int(meta.get("H4READ0", 0) or 0)
    if box is None:
        box = (slice(None), slice(None))
    # Track the illuminated pixels, not the typical ones. Which pixels those
    # are is decided once, from the read carrying the most signal, so every
    # read is then measured over the same set and the rates are comparable.
    planes = [np.asarray(cube.getReadArray(i), np.float32)[box]
              for i in range(nRead)]
    deltas = [planes[0]] + [planes[i] - planes[i - 1] for i in range(1, nRead)]
    brightest = max(range(nRead), key=lambda i: np.percentile(deltas[i],
                                                              litPercentile))
    threshold = np.percentile(deltas[brightest], litPercentile)
    lit = deltas[brightest] >= threshold
    if not lit.any():
        return None
    rate = np.array([float(np.median(d[lit]))/frameTime for d in deltas])
    ends = np.array([(i + read0 + 2)*frameTime for i in range(nRead)])

    usable = np.zeros(nRead, dtype=bool)
    usable[firstUsableRead:] = True
    peak = rate[usable].max() if usable.any() else 0.0
    if peak <= minPlateau:
        return None
    plateau = float(np.median(rate[usable & (rate > LIT_FRACTION*peak)]))
    full = np.where(usable & (rate > LIT_FRACTION*plateau))[0]
    # A genuine plateau is held across several reads. One read alone above the
    # threshold means the "plateau" is an artefact -- most often the first
    # read, which carries reset and settling behaviour and can exceed every
    # illuminated read on a faint flat. Taking it as the plateau reports the
    # lamp coming on in the first read, which it did not.
    if len(full) < minPlateauReads:
        return None
    firstFull, lastFull = int(full[0]), int(full[-1])

    lit = np.where(usable & (rate > EDGE_FRACTION*plateau))[0]
    first = int(lit[0]) if len(lit) and lit[0] <= firstFull else firstFull
    on = float(ends[first] - frameTime*float(np.clip(rate[first]/plateau, 0, 1)))
    # The lamp-off edge is taken as on + EXPTIME, not from the flux. EXPTIME
    # is the commanded lamp time and the leading edge is now solid, whereas a
    # flux-derived trailing edge is caught between two contaminants pulling
    # opposite ways: prompt persistence in the straddling read biases it late,
    # while the brightest pixels -- the ones the lit mask selects -- saturate
    # late in a bright ramp and bias it early, by as much as a whole read.
    # `offFlux` keeps the flux estimate so the two can be compared.
    exptime = float(meta["EXPTIME"])
    off = float(on + exptime)
    offFlux = float(ends[lastFull])
    offUpper = (offFlux if lastFull == nRead - 1 else
                offFlux + frameTime*float(np.clip(rate[lastFull + 1]/plateau,
                                                  0, 1)))

    return IlluminationWindow(
        on=on, off=off, offFlux=offFlux, offUpper=float(offUpper),
        plateau=plateau,
        frameTime=frameTime, nRead=nRead, rate=rate, ends=ends,
        leadingDark=int(np.floor(on/frameTime - 1.0 + 1e-6)),
        trailingDark=int(np.floor((float(meta["DARKTIME"]) - off)/frameTime
                                  + 1e-6)),
        lampStillOn=(lastFull == nRead - 1),
        exptime=exptime, darktime=float(meta["DARKTIME"]),
        mjdStart=float(meta["MJD-STR"]))


def measureEdgeByScan(cube, box=None, litPercentile=99.0, nBins=16):
    """Locate the lamp transitions within a read, using the readout scan.

    A read is not an instantaneous sample. The detector is scanned over the
    frame time, so a lamp switching mid-read leaves a *spatial* boundary: the
    part scanned before the switch shows nothing, the part after shows light,
    and in between the ratio to a fully illuminated read rises linearly with
    scan position. On 146285/n1 that ratio runs 0.006 below image column 2000
    and then climbs steadily to 0.49 at the last column.

    That structure carries far more timing information than a whole-read
    median, which averages the lit and unlit halves together and can only ever
    place an edge to the nearest read. Fitting the ramp gives the transition to
    a fraction of a read.

    For a pixel scanned at time ``t`` within the read, with the lamp switching
    on at ``T``, the fraction of its integration that was illuminated is
    ``(t - T)/frameTime``, clipped to [0, 1]. Scan position maps to image
    column, so fitting ratio against column and solving for the zero crossing
    gives ``T``.

    Returns
    -------
    result : `dict` or `None`
        ``onFraction`` is where in the read the lamp came on, as a fraction of
        the frame time; ``onColumn`` the corresponding scan position;
        ``slope`` the fitted rise per column, which should be about 1/ncols if
        the scan spans one frame time; ``readIndex`` the read it happened in.
    """
    meta = cube.metadata
    nRead = cube.getNumReads()
    frameTime = float(meta["W_H4FRMT"])
    if box is None:
        box = (slice(None), slice(None))
    planes = [np.asarray(cube.getReadArray(i), np.float32)[box]
              for i in range(nRead)]
    deltas = [planes[0]] + [planes[i] - planes[i - 1] for i in range(1, nRead)]
    order = sorted(range(nRead),
                   key=lambda i: np.percentile(deltas[i], litPercentile))
    full = order[-1]
    lit = deltas[full] >= np.percentile(deltas[full], litPercentile)
    if not lit.any():
        return None
    ncols = deltas[full].shape[1]
    edges = np.linspace(0, ncols, nBins + 1).astype(int)

    def profile(delta):
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = np.zeros_like(lit)
            m[:, lo:hi] = True
            m &= lit
            if m.sum() < 50:
                out.append(np.nan)
                continue
            ref = float(np.median(deltas[full][m]))
            out.append(float(np.median(delta[m]))/ref if ref else np.nan)
        return np.array(out)

    centres = 0.5*(edges[:-1] + edges[1:])
    for index in range(nRead):
        if index == full:
            continue
        ratio = profile(deltas[index])
        rising = np.isfinite(ratio) & (ratio > 0.05) & (ratio < 0.95)
        if rising.sum() < 3:
            continue
        slope, intercept = np.polyfit(centres[rising], ratio[rising], 1)
        if slope <= 0:
            continue
        onColumn = -intercept/slope
        if not 0 <= onColumn <= ncols:
            continue
        return dict(readIndex=index, onColumn=float(onColumn),
                    onFraction=float(onColumn/ncols), slope=float(slope),
                    expectedSlope=1.0/ncols, frameTime=frameTime,
                    ratio=ratio, centres=centres)
    return None
