"""Deprecated, disabled compatibility shim for the legacy NIR linearity calib.

The ``NirLinearity`` storage class name (pytype
``lsst.obs.pfs.nirLinearity.NirLinearity``) is retained so that butler
registries which already reference it remain valid across stack versions.
The class itself is **disabled**: the per-pixel H4 nonlinearity correction is
now provided by :mod:`lsst.obs.pfs.h4Linearity` (storage class ``H4Linearity``).
Any attempt to construct or read/write this legacy class raises.
"""

from __future__ import annotations

__all__ = ("NirLinearity",)

_DISABLED = (
    "lsst.obs.pfs.nirLinearity.NirLinearity is deprecated and disabled; "
    "use the h4Linearity calibration (lsst.obs.pfs.h4Linearity, storage class "
    "'H4Linearity') instead."
)


class NirLinearity:
    """Disabled legacy NIR linearity calibration.

    The name is retained only for butler schema compatibility; every entry
    point raises `RuntimeError`.
    """

    def __init__(self, *args, **kwargs):
        raise RuntimeError(_DISABLED)

    @classmethod
    def empty(cls, *args, **kwargs):
        raise RuntimeError(_DISABLED)

    @classmethod
    def fromFits(cls, *args, **kwargs):
        raise RuntimeError(_DISABLED)

    @classmethod
    def readFits(cls, *args, **kwargs):
        raise RuntimeError(_DISABLED)

    def toFits(self, *args, **kwargs):
        raise RuntimeError(_DISABLED)

    def writeFits(self, *args, **kwargs):
        raise RuntimeError(_DISABLED)
