"""Simulation framework for H4 ramp defects (CRs, ASIC glitches).

Generates synthetic cumulative-ADU ramps with controllable noise and
injects known defects so detection algorithms can be exercised against
ground truth, in either the linearized or the raw-domain view.

Two ramp builders are provided:

- `makeRamp` — linear response, no nonlinearity. Both CRs and glitches
  inject directly into the cumulative cube. Use this when you only need
  to exercise post-linearization detection (e.g.
  `cr.iterativeUtrDetectAndRepair`).
- `makeRawRamp` — applies a forward `Nonlinearity` to the true-linear
  signal so the result is a "raw" (detector-domain) cube. CRs are
  injected in true space (before the ADC nonlinearity); ASIC glitches
  are injected in measured space (after, as digital electronics
  artifacts).

`makeRawAndLinearRamps` is a convenience that returns both a raw cube
and its linearized counterpart (the analytic inverse of the same
nonlinearity), so detection algorithms can be run on either side and
compared against the same injected truth.

Pixel response (constant rate, simple model)::

    trueLin[k] = bias + rate * k    (+ Poisson noise if enabled)
    measured   = nl.forward(trueLin) + N(0, readNoise)    (raw path)

Defect injection
----------------

CR (cosmic-ray, step up in cumulative)::

    trueLin[k:, y, x] += amount      for k >= cr.read

Persists through the end of the ramp. Multiple CRs per pixel allowed.

ASIC glitch (single-read offset, digital)::

    measured[k, y, x] += amount      at exactly one read

Returns to the unglitched cumulative at the next read. In delta space
(``np.diff(cum, axis=0)``) it shows as a matched +A / −A pair at
adjacent read indices. Amplitudes are conventionally digital — ±2**N
for some N (commonly N ∈ {10..13} → 1024..8192 ADU), reflecting bit-
flips in the ASIC.

Pure numpy + dataclasses; no LSST stack imports.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Union

import numpy as np


@dataclass
class RampParams:
    """Synthetic ramp generator parameters.

    Attributes
    ----------
    nReads : int
        Number of reads in the cube.
    H, W : int
        Spatial dimensions.
    rate : float
        Per-read flux (ADU/read), constant across the ramp.
    readNoise : float
        Per-read Gaussian noise (ADU).
    bias : float
        Constant pedestal added to every read (ADU).
    poisson : bool
        Include shot noise (Gaussian approximation with variance =
        cumulative signal).
    """

    nReads: int = 45
    H: int = 16
    W: int = 16
    rate: float = 50.0
    readNoise: float = 5.0
    bias: float = 0.0
    poisson: bool = True


@dataclass(frozen=True)
class CR:
    """Cosmic-ray hit (step up in cumulative starting at read ``read``)."""

    y: int
    x: int
    read: int
    amount: float


@dataclass(frozen=True)
class AsicGlitch:
    """Digital ASIC glitch (single-read offset at one read).

    Amplitude is signed: positive for up-glitches, negative for down-glitches.
    Typically a power of 2 or sum of powers of 2.
    """

    y: int
    x: int
    read: int
    amount: float


@dataclass(frozen=True)
class Nonlinearity:
    """Quadratic-compressive detector nonlinearity for sim.

    Models the H4 detector's response as a quadratic deviation from
    linear:

        forward (true_linear → measured):
            measured = t − α * t² / qMax

        inverse (measured → true_linear, analytic root):
            t = (1 − sqrt(max(0, 1 − 4 * α * m / qMax))) / (2 * α / qMax)

    The inverse is exact in float64 arithmetic for t < qMax / (2α) (the
    monotonic regime). With ``alpha=0.05`` and ``qMax=60000`` the
    response is 5% compressed at full well, the monotonic limit is
    600,000 ADU (well above any real cube value), and a complete
    round-trip ``inverse(forward(t))`` returns ``t`` to ~1 ADU
    precision in float32.

    ``alpha=0`` is the linear case (forward and inverse both identity);
    use this to disable nonlinearity without removing the parameter.
    """

    alpha: float = 0.05
    qMax: float = 60000.0

    def forward(self, t: np.ndarray) -> np.ndarray:
        """True linear cumulative → measured (raw-domain) cumulative."""
        if self.alpha == 0.0:
            return np.asarray(t, dtype=np.float32).copy()
        t64 = np.asarray(t, dtype=np.float64)
        m = t64 - self.alpha * t64 * t64 / self.qMax
        return m.astype(np.float32, copy=False)

    def inverse(self, m: np.ndarray) -> np.ndarray:
        """Measured (raw-domain) cumulative → true linear cumulative."""
        if self.alpha == 0.0:
            return np.asarray(m, dtype=np.float32).copy()
        m64 = np.asarray(m, dtype=np.float64)
        a = self.alpha / self.qMax
        disc = np.maximum(0.0, 1.0 - 4.0 * a * m64)
        t = (1.0 - np.sqrt(disc)) / (2.0 * a)
        return t.astype(np.float32, copy=False)


def _asGenerator(rng):
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(rng)


def makeRamp(
    params: RampParams,
    crs: Optional[Sequence[CR]] = None,
    glitches: Optional[Sequence[AsicGlitch]] = None,
    rng: Optional[Union[int, np.random.Generator]] = None,
) -> np.ndarray:
    """Build a synthetic cumulative-ADU cube with optional injected defects.

    Parameters
    ----------
    params : RampParams
    crs : sequence of CR, optional
        Cosmic-ray injections.
    glitches : sequence of AsicGlitch, optional
        ASIC-glitch injections.
    rng : np.random.Generator | int | None
        Reproducible random state. ``int`` seeds a fresh generator,
        ``None`` uses an unseeded one.

    Returns
    -------
    cube : np.ndarray
        Shape ``(nReads, H, W)``, float32. Cumulative ADU per read.
    """
    rng = _asGenerator(rng)

    nReads, H, W = params.nReads, params.H, params.W
    # Per-read flux increments. Each read accumulates 'rate' ADU of new
    # signal (plus Poisson sqrt(rate) noise if enabled). CRs are charge
    # deposits at a specific read; they inject into the per-read flux.
    if params.poisson and params.rate > 0:
        perRead = rng.normal(params.rate,
                             np.sqrt(params.rate),
                             (nReads, H, W))
    else:
        perRead = np.full((nReads, H, W), params.rate, dtype=np.float64)

    if crs:
        for c in crs:
            perRead[c.read, c.y, c.x] += c.amount

    # Cumulate to get the true-charge cumulative, then add per-read read
    # noise (independent at each readout) and ASIC glitches (digital
    # post-ADC offsets).
    cum = np.cumsum(perRead, axis=0) + params.bias
    cum = cum + rng.normal(0.0, params.readNoise, (nReads, H, W))
    cum = cum.astype(np.float32)

    if glitches:
        for g in glitches:
            cum[g.read, g.y, g.x] += np.float32(g.amount)

    return cum


def utrRate(
    cube: np.ndarray,
    readMask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Per-pixel UTR (Up The Ramp) rate from a cumulative cube.

    Computes the least-squares slope of cumulative-vs-read for each
    pixel. With ``readMask`` you can exclude specific reads from the fit
    (e.g. CR-affected reads during an iterative detector); without it
    the fit uses all reads with equal weight.

    Parameters
    ----------
    cube : np.ndarray
        Shape ``(N, H, W)``. Cumulative ADU per read.
    readMask : np.ndarray, optional
        Boolean shape ``(N, H, W)`` (per-pixel per-read), or ``(N,)``
        (whole-read mask), or ``(H, W)`` (whole-pixel mask). True means
        the read is used in the fit. If omitted, all reads are used.

    Returns
    -------
    rate : np.ndarray
        Shape ``(H, W)``, float32. ADU per read.

    Notes
    -----
    For PFS H4 the production rate estimator (``PfsIsrTask.calcUTRrates``)
    applies optimal weights derived from a noise model. The simple
    unweighted LSQ slope used here is adequate for sim closure tests
    and for the iterative-detection rate reference.
    """
    nReads, H, W = cube.shape
    k = np.arange(nReads, dtype=np.float64)
    cube64 = cube.astype(np.float64, copy=False)

    if readMask is None:
        # Unweighted closed form: slope = sum((k-mean)*y) / sum((k-mean)**2)
        kMean = k.mean()
        kCentered = k - kMean
        kVar = float((kCentered * kCentered).sum())
        slope = np.einsum('k,khw->hw', kCentered, cube64) / kVar
        return slope.astype(np.float32, copy=False)

    # With a mask, broadcast it to (N, H, W) bool and do a per-pixel
    # weighted fit (weights are 0/1). Same closed form, just with weights.
    mask = np.asarray(readMask, dtype=bool)
    if mask.shape == (nReads,):
        mask = np.broadcast_to(mask[:, None, None], (nReads, H, W))
    elif mask.shape == (H, W):
        mask = np.broadcast_to(mask[None], (nReads, H, W))
    elif mask.shape != (nReads, H, W):
        raise ValueError(f"readMask shape {readMask.shape} not compatible with "
                         f"cube {cube.shape}.")
    w = mask.astype(np.float64)
    kCube = np.broadcast_to(k[:, None, None], (nReads, H, W))
    sumW = w.sum(axis=0)
    sumK = (w * kCube).sum(axis=0)
    sumY = (w * cube64).sum(axis=0)
    sumKK = (w * kCube * kCube).sum(axis=0)
    sumKY = (w * kCube * cube64).sum(axis=0)
    denom = sumW * sumKK - sumK * sumK
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = np.where(denom > 0, (sumW * sumKY - sumK * sumY) / denom, 0.0)
    return slope.astype(np.float32, copy=False)


def makeRawRamp(
    params: RampParams,
    nonlinearity: Optional[Nonlinearity] = None,
    crs: Optional[Sequence[CR]] = None,
    asicGlitches: Optional[Sequence[AsicGlitch]] = None,
    rng: Optional[Union[int, np.random.Generator]] = None,
) -> np.ndarray:
    """Build a synthetic *raw-domain* (post-ADC) cumulative-ADU cube.

    Physically motivated injection order:

    1. ``trueLin = rate * k + bias`` (noiseless linear ramp).
    2. (optional) add Poisson shot noise to ``trueLin`` (in photon-charge space).
    3. Inject ``crs`` into ``trueLin`` (charge deposits, before the ADC).
    4. Apply ``nonlinearity.forward`` (true → measured).
    5. Add per-read Gaussian read noise (in measured space).
    6. Inject ``asicGlitches`` (digital, single-read; post-ADC artifacts).

    Setting ``nonlinearity=None`` is equivalent to ``alpha=0``: the cube
    is linear (same as `makeRamp`, but glitches still inject post-noise).

    Parameters
    ----------
    params : RampParams
    nonlinearity : Nonlinearity, optional
    crs : sequence of CR, optional
        Injected before the nonlinearity, i.e. on the true-linear cube.
    asicGlitches : sequence of AsicGlitch, optional
        Injected after the nonlinearity, i.e. on the measured cube.
    rng : int | np.random.Generator | None

    Returns
    -------
    measured : np.ndarray
        ``(nReads, H, W)`` float32 cube in raw (detector-domain) ADU.
    """
    rng = _asGenerator(rng)

    nReads, H, W = params.nReads, params.H, params.W
    # Per-read flux increments with Poisson noise (if enabled); CRs add
    # charge to a specific read.
    if params.poisson and params.rate > 0:
        perRead = rng.normal(params.rate,
                             np.sqrt(params.rate),
                             (nReads, H, W))
    else:
        perRead = np.full((nReads, H, W), params.rate, dtype=np.float64)

    if crs:
        for c in crs:
            perRead[c.read, c.y, c.x] += c.amount

    trueLin = np.cumsum(perRead, axis=0) + params.bias

    if nonlinearity is None or nonlinearity.alpha == 0.0:
        measured = trueLin.astype(np.float32, copy=False)
    else:
        measured = nonlinearity.forward(trueLin)

    measured = measured + rng.normal(
        0.0, params.readNoise, (nReads, H, W)
    ).astype(np.float32)

    if asicGlitches:
        for g in asicGlitches:
            measured[g.read, g.y, g.x] += np.float32(g.amount)

    return measured.astype(np.float32, copy=False)


def makeRawAndLinearRamps(
    params: RampParams,
    nonlinearity: Nonlinearity,
    crs: Optional[Sequence[CR]] = None,
    asicGlitches: Optional[Sequence[AsicGlitch]] = None,
    rng: Optional[Union[int, np.random.Generator]] = None,
) -> tuple:
    """Return ``(raw, linearized)`` cubes for one consistent ground truth.

    The raw cube is the output of `makeRawRamp` with the supplied
    nonlinearity; the linearized cube is ``nonlinearity.inverse(raw)``.
    Both share the same injected CRs (in true space) and ASIC glitches
    (in measured space), so detection on either side can be checked
    against the same truth catalog. In the linearized cube:

    - CRs are clean step-ups (the round-trip nonlinearity recovers them).
    - ASIC glitches are *slightly deformed* — their amplitude after the
      inverse depends on the local slope of the inverse curve at the
      glitched read. The user can compare this deformation to glitch
      detection performance.

    Returns
    -------
    raw, linearized : np.ndarray
        Both ``(nReads, H, W)`` float32.
    """
    raw = makeRawRamp(params, nonlinearity=nonlinearity, crs=crs,
                      asicGlitches=asicGlitches, rng=rng)
    linearized = nonlinearity.inverse(raw)
    return raw, linearized


# ----- Convenience generators -----

def digitalGlitchAmounts(
    n: int,
    bits: Sequence[int] = (10, 11, 12, 13),
    signed: bool = True,
    rng: Optional[Union[int, np.random.Generator]] = None,
) -> np.ndarray:
    """Generate ``n`` digital glitch amplitudes drawn from ``±2**bit``.

    Default bits 10..13 give magnitudes 1024, 2048, 4096, 8192 ADU,
    which roughly matches the observed range on PFS H4 channel-24
    glitches.

    Returns
    -------
    amounts : np.ndarray
        Shape ``(n,)``, float32. Signed if ``signed=True``.
    """
    rng = _asGenerator(rng)
    bits_arr = np.asarray(bits)
    chosen = rng.choice(bits_arr, size=n)
    amounts = (2.0 ** chosen).astype(np.float32)
    if signed:
        signs = rng.choice([-1.0, 1.0], size=n).astype(np.float32)
        amounts = amounts * signs
    return amounts


def randomCRs(
    n: int,
    H: int,
    W: int,
    nReads: int,
    amountRange: tuple = (50.0, 5000.0),
    readRange: Optional[tuple] = None,
    rng: Optional[Union[int, np.random.Generator]] = None,
) -> list:
    """Generate ``n`` random CR hits with log-uniform amplitudes.

    ``readRange=None`` means ``(1, nReads - 1)`` (skip first and last
    reads so the step is fully visible in cumulative space).
    """
    rng = _asGenerator(rng)
    if readRange is None:
        readRange = (1, nReads - 1)
    ys = rng.integers(0, H, size=n)
    xs = rng.integers(0, W, size=n)
    reads = rng.integers(readRange[0], readRange[1], size=n)
    logLo, logHi = np.log(amountRange[0]), np.log(amountRange[1])
    amounts = np.exp(rng.uniform(logLo, logHi, size=n))
    return [
        CR(int(y), int(x), int(r), float(a))
        for y, x, r, a in zip(ys, xs, reads, amounts)
    ]


def randomAsicGlitches(
    n: int,
    H: int,
    W: int,
    nReads: int,
    bits: Sequence[int] = (10, 11, 12, 13),
    readRange: Optional[tuple] = None,
    rng: Optional[Union[int, np.random.Generator]] = None,
) -> list:
    """Generate ``n`` random ASIC glitches with digital amplitudes.

    ``readRange=None`` means ``(1, nReads - 2)`` — skip both endpoints so
    the glitch and its return-to-baseline are both inside the ramp.
    """
    rng = _asGenerator(rng)
    if readRange is None:
        readRange = (1, nReads - 2)
    ys = rng.integers(0, H, size=n)
    xs = rng.integers(0, W, size=n)
    reads = rng.integers(readRange[0], readRange[1], size=n)
    amounts = digitalGlitchAmounts(n, bits=bits, signed=True, rng=rng)
    return [
        AsicGlitch(int(y), int(x), int(r), float(a))
        for y, x, r, a in zip(ys, xs, reads, amounts)
    ]
