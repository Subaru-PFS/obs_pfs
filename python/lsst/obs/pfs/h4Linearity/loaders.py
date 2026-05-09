"""Development convenience loaders. Production callers supply their own loaders."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .types import Ramp


def loadNpz(path: str | Path) -> tuple[Ramp, np.ndarray]:
    """Load a ``.npz`` with ``deltas`` and ``photodiode`` arrays.

    The on-disk format stores per-read deltas with an implicit read0 = 0.
    This loader prepends the zero read and accumulates, yielding ``N+1``
    cumulative reads from ``N`` deltas. The caller is expected to apply
    the photodiode correction before passing the ramp into
    :func:`nirLinearity.fit.fit`.
    """
    path = Path(path)
    with np.load(path) as data:
        deltas = np.asarray(data["deltas"], dtype=np.float32)
        photodiode = np.asarray(data["photodiode"])
    nDeltas, h, w = deltas.shape
    reads = np.empty((nDeltas + 1, h, w), dtype=np.float32)
    reads[0] = 0.0
    np.cumsum(deltas, axis=0, out=reads[1:])
    return Ramp(reads=reads), photodiode
