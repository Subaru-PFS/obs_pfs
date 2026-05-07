"""Development convenience loaders. Production callers supply their own loaders."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .types import Ramp


def loadNpz(path: str | Path) -> tuple[Ramp, np.ndarray]:
    """Load a ``.npz`` with ``deltas`` and ``photodiode`` arrays.

    The on-disk format stores per-read deltas; this loader converts to
    cumulative flux (``np.cumsum``) before constructing the :class:`Ramp`.
    The caller is expected to apply the photodiode correction before
    passing the ramp into :func:`nirLinearity.fit.fit`.
    """
    path = Path(path)
    with np.load(path) as data:
        deltas = np.asarray(data["deltas"], dtype=np.float32)
        photodiode = np.asarray(data["photodiode"])
    reads = np.cumsum(deltas, axis=0)
    return Ramp(reads=reads), photodiode
