"""Row-banded traversal of the H4 ramp cubes.

Several whole-cube passes in the H4 ISR reduce along the time axis one pixel
at a time: the CR candidate scan, the BAD-pixel outlier count, the final
UTR-weight accumulation, the glitch-height maximum and the rate-stability
segment sums. Each is independent per pixel, so cutting the cube into bands of
image rows cannot change any pixel's result -- what it changes is the working
set. A ``[..., k]`` slice of a 9 GB cube reads one element per cache line and
holds the whole cube live; a band small enough to stay resident is read once
and hit repeatedly, and any temporary built from it is bounded by the band
rather than by the cube.

Two things a caller must respect for the result to stay bit-identical:

- Reduce only along the time axis. Banding changes which pixels are handled
  together, never which values a pixel sees.
- Keep any float accumulation *over reads* in its original order. float32
  addition is not associative, so replacing an ascending-k loop with a single
  reduction changes the last bits.
"""

import numpy as np

__all__ = ("CHUNK_BYTES", "rowChunks")

#: Working-set target, in bytes, for a banded pass. Small enough that a band
#: of the delta cube stays in cache, large enough that the per-band numpy
#: overhead is negligible.
CHUNK_BYTES = 16 << 20


def rowChunks(shape, itemsize, budget=None):
    """Yield ``(lo, hi)`` row bands of a cube that fit within ``budget``.

    Parameters
    ----------
    shape : `tuple`
        Shape of the cube being traversed; axis 0 is banded.
    itemsize : `int`
        Bytes per element of that cube.
    budget : `int`, optional
        Target bytes per band. Defaults to `CHUNK_BYTES`.

    Yields
    ------
    lo, hi : `int`
        Half-open row range. The bands tile ``range(shape[0])`` in order.
    """
    if budget is None:
        budget = CHUNK_BYTES
    nRows = shape[0]
    perRow = int(np.prod(shape[1:])) * itemsize
    step = max(1, int(budget) // max(1, perRow))
    for lo in range(0, nRows, step):
        yield lo, min(lo + step, nRows)
