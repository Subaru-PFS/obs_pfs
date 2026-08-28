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
        Lamp-off: the end of the last fully illuminated read. A lower bound.
    offUpper : `float`
        Lamp-off ignoring prompt persistence in the straddling read. An upper
        bound; the truth lies between this and ``off``.
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


def measureIllumination(cube, box=None, minPlateau=5.0):
    """Measure a ramp's illumination window from its flux.

    Parameters
    ----------
    cube : `lsst.obs.pfs.ImageCube`
        The ramp, as read from a ``rawISRCube``.
    box : `tuple` of `slice`, optional
        Region to measure over. Defaults to the whole detector, which is
        usually wasteful; a region on illuminated fibers is enough and much
        faster.
    minPlateau : `float`, optional
        Reject ramps whose peak rate is below this (e-/s) as too faint to
        locate an edge. Arcs frequently are.

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
    cumulative = [float(np.median(np.asarray(cube.getReadArray(i), np.float32)[box]))
                  for i in range(nRead)]
    rate = np.diff([0.0] + cumulative)/frameTime
    ends = np.array([(i + read0 + 2)*frameTime for i in range(nRead)])

    peak = rate.max()
    if peak <= minPlateau:
        return None
    plateau = float(np.median(rate[rate > LIT_FRACTION*peak]))
    full = np.where(rate > LIT_FRACTION*plateau)[0]
    if len(full) == 0:
        return None
    firstFull, lastFull = int(full[0]), int(full[-1])

    lit = np.where(rate > EDGE_FRACTION*plateau)[0]
    first = int(lit[0]) if len(lit) and lit[0] <= firstFull else firstFull
    on = float(ends[first] - frameTime*float(np.clip(rate[first]/plateau, 0, 1)))
    off = float(ends[lastFull])
    offUpper = (off if lastFull == nRead - 1 else
                off + frameTime*float(np.clip(rate[lastFull + 1]/plateau, 0, 1)))

    return IlluminationWindow(
        on=on, off=off, offUpper=float(offUpper), plateau=plateau,
        frameTime=frameTime, nRead=nRead, rate=rate, ends=ends,
        leadingDark=int(np.floor(on/frameTime - 1.0 + 1e-6)),
        trailingDark=int(np.floor((float(meta["DARKTIME"]) - off)/frameTime
                                  + 1e-6)),
        lampStillOn=(lastFull == nRead - 1),
        exptime=float(meta["EXPTIME"]), darktime=float(meta["DARKTIME"]),
        mjdStart=float(meta["MJD-STR"]))
