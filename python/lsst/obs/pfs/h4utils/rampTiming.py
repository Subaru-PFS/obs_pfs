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

Two cautions, both learned the hard way:

- Use a *low* threshold to decide which read carries an edge. Reads only
  9-29% illuminated are common; treating them as dark truncates the answer to
  a read boundary and can make cameras look a whole read apart.
- The read straddling lamp-off also contains prompt persistence -- on a bright
  flat that is ~16% of the plateau -- so the trailing edge cannot be placed
  exactly from flux alone. `IlluminationWindow.off` is the start of the first
  unilluminated read, a lower bound, and ``offUpper`` is the naive flux
  estimate; the truth lies between.

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

    A lamp switching partway through a read leaves every pixel with the same
    *fraction* of its illuminated rate, whatever its brightness. Persistence
    does not: it scales with what the pixel received from earlier exposures,
    so the brightest pixels -- which took the most dose then -- show the
    largest fraction, and the faint ones show none.

    On 145076 read 0 runs at 47% of plateau on the brightest 0.1% of pixels,
    45% on the brightest 1%, 19% on the brightest 10% and 0% on the median
    pixel. That is a persistence signature, and reading it as illumination
    moves the inferred lamp-on a whole read early.

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
                        firstUsableRead=1):
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
        Ignore reads before this when locating the lamp. Read 0 is not usable:
        it carries persistence released by the preceding exposures, which
        follows the same fiber traces as the illumination and so is easily
        mistaken for it. On 145076 read 0 sits at 45% of plateau on the
        brightest pixels while the median pixel sees nothing -- a
        dose-proportional signature, not a lamp. The commanding never turns
        the lamp on that early, so nothing is lost by skipping it, and
        `leadingDark` still counts it as unilluminated, which for the lamp it
        is.

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
