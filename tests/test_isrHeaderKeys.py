import unittest

import numpy as np

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.utils.tests
from lsst.obs.pfs import h4Linearity
from lsst.obs.pfs.isrTask import (
    _makeInternalMask, _projectInternalMask, _stampRampMetadata,
)


class StampRampMetadataTestCase(lsst.utils.tests.TestCase):
    """`_stampRampMetadata` writes the expected H4* keys on the
    exposure's metadata so downstream consumers can identify partial-ramp
    products from a full-ramp postISRCCD.
    """

    def _makeExposure(self):
        return afwImage.ExposureF(geom.Extent2I(4, 4))

    def testFullRampWithUTR(self):
        exp = self._makeExposure()
        _stampRampMetadata(exp, r0=0, r1=39, nTotal=40, appliedUTR=True)
        md = exp.getMetadata()
        self.assertEqual(int(md.getScalar('H4READ0')), 0)
        self.assertEqual(int(md.getScalar('H4READ1')), 39)
        self.assertEqual(int(md.getScalar('H4NREAD')), 40)
        self.assertEqual(int(md.getScalar('H4NTOT')), 40)
        self.assertTrue(bool(md.getScalar('H4UTRWT')))

    def testPartialRampMatchesFirstHalfRange(self):
        """The first-half slice carries the trimmed range, not the full ramp's."""
        exp = self._makeExposure()
        _stampRampMetadata(exp, r0=1, r1=22, nTotal=40, appliedUTR=True)
        md = exp.getMetadata()
        self.assertEqual(int(md.getScalar('H4READ0')), 1)
        self.assertEqual(int(md.getScalar('H4READ1')), 22)
        self.assertEqual(int(md.getScalar('H4NREAD')), 22)
        self.assertEqual(int(md.getScalar('H4NTOT')), 40)

    def testCdsLayout(self):
        """A CDS (non-UTR-weighted) product marks itself via H4UTRWT=False."""
        exp = self._makeExposure()
        _stampRampMetadata(exp, r0=0, r1=39, nTotal=40, appliedUTR=False)
        self.assertFalse(bool(exp.getMetadata().getScalar('H4UTRWT')))

    def testCommentsArePresent(self):
        """Each key carries a one-line comment for the FITS header."""
        exp = self._makeExposure()
        _stampRampMetadata(exp, r0=0, r1=39, nTotal=40, appliedUTR=True)
        md = exp.getMetadata()
        for key in ('H4READ0', 'H4READ1', 'H4NREAD', 'H4NTOT', 'H4UTRWT'):
            self.assertTrue(md.getComment(key), msg=f'{key} should carry a comment')


class MakeInternalMaskTestCase(lsst.utils.tests.TestCase):
    """`_makeInternalMask` seeds the H4 ISR internal mask with
    BORDER_PIX on the outer ring plus MASKED_BY_INPUT from the defects
    calib and whatever bits the linearity calib's badPixelMask carries.
    """

    def testBorderRingSet(self):
        m = _makeInternalMask((20, 24))
        # Outer 4 rows/cols carry BORDER_PIX; interior is zero.
        for region in (m[:4, :], m[-4:, :], m[:, :4], m[:, -4:]):
            self.assertTrue(((region & h4Linearity.BORDER_PIX)
                             == h4Linearity.BORDER_PIX).all())
        self.assertEqual(int(m[4:-4, 4:-4].sum()), 0)

    def testBorderWidthRespected(self):
        m = _makeInternalMask((16, 16), borderWidth=2)
        self.assertTrue(((m[:2, :] & h4Linearity.BORDER_PIX)
                         == h4Linearity.BORDER_PIX).all())
        self.assertEqual(int(m[2:-2, 2:-2].sum()), 0)

    def testInternalMaskDtypeUint16(self):
        m = _makeInternalMask((16, 16))
        self.assertEqual(m.dtype, np.uint16)


class ProjectInternalMaskTestCase(lsst.utils.tests.TestCase):
    """`_projectInternalMask` lifts internal-mask bits onto Exposure.mask
    planes following the canonical projection rule.
    """

    def _makeExposure(self, H=16, W=16):
        return afwImage.ExposureF(geom.Extent2I(W, H))

    def testBorderProjectsToBADOnly(self):
        """BORDER_PIX → BAD (no separate BORDER plane published)."""
        exp = self._makeExposure(H=12, W=12)
        internal = _makeInternalMask((12, 12))
        _projectInternalMask(exp, internal)
        badBit = exp.mask.getPlaneBitMask('BAD')
        # Outer ring all BAD; interior clean.
        for region in (exp.mask.array[:4, :], exp.mask.array[-4:, :],
                       exp.mask.array[:, :4], exp.mask.array[:, -4:]):
            self.assertTrue(((region & badBit) == badBit).all())
        self.assertEqual(int(exp.mask.array[4:-4, 4:-4].sum()), 0)
        # No BORDER plane registered.
        self.assertNotIn('BORDER', exp.mask.getMaskPlaneDict())

    def testDarkDefectProjection(self):
        """MASKED_BY_INPUT → DARK_DEFECT + BAD."""
        exp = self._makeExposure(H=8, W=8)
        internal = np.zeros((8, 8), dtype=np.uint16)
        internal[3, 4] |= h4Linearity.MASKED_BY_INPUT
        _projectInternalMask(exp, internal)
        darkDefectBit = exp.mask.getPlaneBitMask('DARK_DEFECT')
        badBit = exp.mask.getPlaneBitMask('BAD')
        self.assertEqual(exp.mask.array[3, 4] & darkDefectBit, darkDefectBit)
        self.assertEqual(exp.mask.array[3, 4] & badBit, badBit)

    def testDeadGroupProjection(self):
        """Any of INSUFFICIENT_POINTS / FIT_FAILED / NON_MONOTONIC →
        LINEARITY_DEFECT + BAD.
        """
        exp = self._makeExposure(H=8, W=8)
        internal = np.zeros((8, 8), dtype=np.uint16)
        internal[1, 1] |= h4Linearity.INSUFFICIENT_POINTS
        internal[2, 2] |= h4Linearity.FIT_FAILED
        internal[3, 3] |= h4Linearity.NON_MONOTONIC
        _projectInternalMask(exp, internal)
        ldBit = exp.mask.getPlaneBitMask('LINEARITY_DEFECT')
        for y, x in [(1, 1), (2, 2), (3, 3)]:
            self.assertEqual(exp.mask.array[y, x] & ldBit, ldBit)

    def testSatProjection(self):
        """ABOVE_VALID_RANGE → SAT + BAD."""
        exp = self._makeExposure(H=8, W=8)
        internal = np.zeros((8, 8), dtype=np.uint16)
        internal[4, 5] |= h4Linearity.ABOVE_VALID_RANGE
        _projectInternalMask(exp, internal)
        satBit = exp.mask.getPlaneBitMask('SAT')
        badBit = exp.mask.getPlaneBitMask('BAD')
        self.assertEqual(exp.mask.array[4, 5] & satBit, satBit)
        self.assertEqual(exp.mask.array[4, 5] & badBit, badBit)

    def testUnstableProjection(self):
        """UNSTABLE bit → UNSTABLE plane + BAD."""
        exp = self._makeExposure(H=8, W=8)
        internal = np.zeros((8, 8), dtype=np.uint16)
        internal[5, 2] |= h4Linearity.UNSTABLE
        _projectInternalMask(exp, internal)
        unstBit = exp.mask.getPlaneBitMask('UNSTABLE')
        badBit = exp.mask.getPlaneBitMask('BAD')
        self.assertEqual(exp.mask.array[5, 2] & unstBit, unstBit)
        self.assertEqual(exp.mask.array[5, 2] & badBit, badBit)

    def testUnclassifiedNotPublished(self):
        """UNCLASSIFIED bit folds into BAD but not into any standalone
        external plane.
        """
        exp = self._makeExposure(H=8, W=8)
        internal = np.zeros((8, 8), dtype=np.uint16)
        internal[6, 6] |= h4Linearity.UNCLASSIFIED
        _projectInternalMask(exp, internal)
        badBit = exp.mask.getPlaneBitMask('BAD')
        self.assertEqual(exp.mask.array[6, 6] & badBit, badBit)

    def testAsicGlitchPlanePublished(self):
        """A corrected glitch carries the published ASIC_GLITCH plane and is
        NOT BAD -- ASIC_GLITCH alone is a usable, informational flag.
        """
        exp = self._makeExposure(H=8, W=8)
        internal = np.zeros((8, 8), dtype=np.uint16)
        internal[2, 5] |= h4Linearity.ASIC_GLITCH
        _projectInternalMask(exp, internal)
        agBit = exp.mask.getPlaneBitMask('ASIC_GLITCH')
        badBit = exp.mask.getPlaneBitMask('BAD')
        self.assertEqual(exp.mask.array[2, 5] & agBit, agBit)
        self.assertEqual(exp.mask.array[2, 5] & badBit, 0)

    def testGlitchMaskedFoldsIntoBAD(self):
        """An uncorrectable glitch (ASIC_GLITCH | GLITCH_MASKED) is flagged in
        the ASIC_GLITCH plane and also BAD.
        """
        exp = self._makeExposure(H=8, W=8)
        internal = np.zeros((8, 8), dtype=np.uint16)
        internal[2, 5] |= h4Linearity.ASIC_GLITCH | h4Linearity.GLITCH_MASKED
        _projectInternalMask(exp, internal)
        agBit = exp.mask.getPlaneBitMask('ASIC_GLITCH')
        badBit = exp.mask.getPlaneBitMask('BAD')
        self.assertEqual(exp.mask.array[2, 5] & agBit, agBit)
        self.assertEqual(exp.mask.array[2, 5] & badBit, badBit)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
