from __future__ import annotations

import lsst.afw.image as afwImage

__all__ = ("CALIB_MASK_PLANES", "ISR_MASK_PLANES", "addObsPfsMaskPlanes")


CALIB_MASK_PLANES = (
    "PARTLY_VIGNETTED",
)
"""Planes a PFS exposure inherits from its calibs rather than from our own code.

Nothing in obs_pfs sets a bit in these; they are claimed only to fix where they
sit. ``lsst.cp.pipe.cpCombine.CalibCombineTask`` stamps ``PARTLY_VIGNETTED`` on
every combined bias and dark, so a visible-arm process acquires the plane the
moment those calibs are read, while a NIR process, which reads neither, never
does. Claiming it up front gives it the same bit on both arms, and claiming it
before `ISR_MASK_PLANES` makes that bit 21 -- where every file already written
carries it.
"""

ISR_MASK_PLANES = (
    "IPC",
    "DARK_DEFECT",
    "LINEARITY_DEFECT",
    "UNSTABLE",
    "RATE_UNSTABLE",
    "ASIC_GLITCH",
)
"""Planes obs_pfs ISR writes into, on either arm.

``IPC`` is claimed by ``pfs.drp.stella`` at bit 16, so naming it here changes
nothing; it is listed so this is the complete statement of what ISR publishes.
"""


def addObsPfsMaskPlanes(exposure: afwImage.Exposure | None = None) -> None:
    """Claim the mask planes every PFS exposure carries.

    Call once per process before any calib is read, and again on each output
    exposure.

    Order is significant. ``Mask.addMaskPlane`` is static -- there is no
    per-instance variant, and spelling it ``exposure.mask.addMaskPlane`` calls
    the same function -- and it hands out the lowest free bit. The name-to-bit
    map a product persists is therefore fixed by the order planes are first
    claimed in the process, which is why claiming them all up front is what
    makes two products comparable.

    Being static is also what makes the per-exposure call work: afw fits a newly
    claimed plane into every live ``Mask``. The check on ``exposure`` is a guard
    rather than a repair, since a mask reaching here with a dictionary of its own
    is a bug -- afw conforms every mask it reads to the process dictionary, so
    there is no legitimate way to acquire one.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`, optional
        Exposure whose mask must carry the full set.

    Raises
    ------
    RuntimeError
        If ``exposure``'s mask does not share the process mask plane dictionary.
    """
    mask = afwImage.Mask() if exposure is None else exposure.mask
    for plane in CALIB_MASK_PLANES + ISR_MASK_PLANES:
        if plane not in mask.getMaskPlaneDict():
            mask.addMaskPlane(plane)
    if exposure is not None:
        planes = exposure.mask.getMaskPlaneDict()
        expected = afwImage.Mask().getMaskPlaneDict()
        if planes != expected:
            raise RuntimeError(
                f"Exposure mask has its own plane dictionary: {planes} != {expected}")
