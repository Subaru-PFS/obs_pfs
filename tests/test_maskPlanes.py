import unittest

import numpy as np
import astropy.io.fits

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.utils.tests

from lsst.obs.pfs.maskPlanes import (
    CALIB_MASK_PLANES, ISR_MASK_PLANES, addObsPfsMaskPlanes,
)

from testUtils import HAS_DRP_STELLA, requireDrpStella

if HAS_DRP_STELLA:
    # drp_stella claims bits 11-20 as it imports, so the expected bit numbers
    # only hold once it has been imported. isrTask imports it in turn.
    import pfs.drp.stella  # noqa: F401
    from lsst.obs.pfs.isrTask import _projectInternalMask

ALL_PLANES = CALIB_MASK_PLANES + ISR_MASK_PLANES

PARTLY_VIGNETTED_BIT = 21
"""The bit every pfsArm and calib already on disk carries PARTLY_VIGNETTED in.

It follows drp_stella's ten planes at 11-20. A package claiming a plane ahead of
us would shift it, which is exactly what these tests exist to catch: the datamodel
`MaskHelper` compares bit numbers between files, so this value is a wire format.
"""


class AddObsPfsMaskPlanesTestCase(lsst.utils.tests.TestCase):
    """`addObsPfsMaskPlanes` claims the whole published set, on the process
    dictionary and on an individual exposure alike.
    """

    def setUp(self):
        addObsPfsMaskPlanes()

    def testEveryPlaneRegistered(self):
        planes = afwImage.Mask().getMaskPlaneDict()
        for name in ALL_PLANES:
            self.assertIn(name, planes)

    def testPlanesReachANewExposure(self):
        exposure = afwImage.ExposureF(geom.Extent2I(4, 4))
        addObsPfsMaskPlanes(exposure)
        planes = exposure.mask.getMaskPlaneDict()
        for name in ALL_PLANES:
            self.assertIn(name, planes)

    def testRepeatedCallsChangeNothing(self):
        """Every output exposure triggers a call, so it has to be idempotent."""
        before = afwImage.Mask().getMaskPlaneDict()
        addObsPfsMaskPlanes()
        exposure = afwImage.ExposureF(geom.Extent2I(4, 4))
        addObsPfsMaskPlanes(exposure)
        addObsPfsMaskPlanes(exposure)
        self.assertEqual(afwImage.Mask().getMaskPlaneDict(), before)

    def testNoTwoPlanesShareABit(self):
        """`MaskHelper` merges silently produce a helper with two names on one
        bit, so the dictionary we hand it must never contain one.
        """
        planes = afwImage.Mask().getMaskPlaneDict()
        self.assertEqual(len(set(planes.values())), len(planes))


@requireDrpStella
class MaskPlaneBitsTestCase(lsst.utils.tests.TestCase):
    """The bit numbers themselves, which are a persisted wire format."""

    def setUp(self):
        addObsPfsMaskPlanes()

    def testPartlyVignettedKeepsItsBit(self):
        self.assertEqual(afwImage.Mask().getMaskPlaneDict()["PARTLY_VIGNETTED"],
                         PARTLY_VIGNETTED_BIT)

    def testIsrPlanesFollowTheCalibPlanes(self):
        """The ISR planes are claimed second, so they cannot displace
        PARTLY_VIGNETTED however many of them there are.
        """
        planes = afwImage.Mask().getMaskPlaneDict()
        calibTop = max(planes[name] for name in CALIB_MASK_PLANES)
        for name in ("DARK_DEFECT", "LINEARITY_DEFECT", "UNSTABLE",
                     "RATE_UNSTABLE", "ASIC_GLITCH"):
            self.assertGreater(planes[name], calibTop)

    def testRoomToSpare(self):
        """afw allows 32 planes; leave the failure legible if we approach it."""
        planes = afwImage.Mask().getMaskPlaneDict()
        self.assertLess(max(planes.values()), afwImage.Mask.getNumPlanesMax())


class DivergentDictionaryTestCase(lsst.utils.tests.TestCase):
    """The per-exposure call refuses a mask carrying a dictionary of its own,
    rather than quietly returning an exposure whose bits mean something else.
    """

    def testGuardRaises(self):
        addObsPfsMaskPlanes()
        exposure = afwImage.ExposureF(geom.Extent2I(4, 4))
        # removeMaskPlane detaches the default and leaves live Masks alone, so
        # this leaves `exposure` holding a plane the process no longer knows.
        afwImage.Mask.removeMaskPlane("ASIC_GLITCH")
        try:
            with self.assertRaises(RuntimeError):
                addObsPfsMaskPlanes(exposure)
        finally:
            addObsPfsMaskPlanes()
        self.assertIn("ASIC_GLITCH", afwImage.Mask().getMaskPlaneDict())


class CalibConformanceTestCase(lsst.utils.tests.TestCase):
    """Claiming PARTLY_VIGNETTED before the calibs are read is lossless.

    A combined bias or dark carries its own numbering, and afw remaps by name as
    it reads. This is why the NIR arm can be handed the plane up front without
    disturbing what the visible arm reads out of its calibs.
    """

    def testForeignNumberingConformsOnRead(self):
        addObsPfsMaskPlanes()
        expected = afwImage.Mask().getMaskPlaneDict()["PARTLY_VIGNETTED"]
        foreign = 11
        self.assertNotEqual(foreign, expected)

        array = np.zeros((4, 4), dtype=np.int32)
        array[1, 2] = 1 << foreign
        hdu = astropy.io.fits.ImageHDU(array)
        hdu.header["MP_BAD"] = 0
        hdu.header["HIERARCH MP_PARTLY_VIGNETTED"] = foreign
        hdus = astropy.io.fits.HDUList([astropy.io.fits.PrimaryHDU(), hdu])
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            hdus.writeto(filename, overwrite=True)
            mask = afwImage.Mask.readFits(filename)

        self.assertEqual(mask.getMaskPlaneDict()["PARTLY_VIGNETTED"], expected)
        self.assertTrue(mask.array[1, 2] & mask.getPlaneBitMask("PARTLY_VIGNETTED"))


@requireDrpStella
class MaskPlaneParityTestCase(lsst.utils.tests.TestCase):
    """A visible-arm and an H4 output end up with the same dictionary, which is
    what lets their pfsArms be merged.
    """

    def testCcdAndH4ExposuresAgree(self):
        ccd = afwImage.ExposureF(geom.Extent2I(8, 8))
        addObsPfsMaskPlanes(ccd)

        h4 = afwImage.ExposureF(geom.Extent2I(8, 8))
        _projectInternalMask(h4, np.zeros((8, 8), dtype=np.uint16))
        addObsPfsMaskPlanes(h4)

        self.assertEqual(ccd.mask.getMaskPlaneDict(), h4.mask.getMaskPlaneDict())


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
