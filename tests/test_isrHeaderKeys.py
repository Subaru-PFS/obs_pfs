import unittest

import numpy as np

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.utils.tests
from lsst.obs.pfs.isrTask import _stampBorderMask, _stampReadRangeMetadata


class StampReadRangeMetadataTestCase(lsst.utils.tests.TestCase):
    """`_stampReadRangeMetadata` writes the expected H4* keys on the
    exposure's metadata so downstream consumers can identify partial-ramp
    products from a full-ramp postISRCCD.
    """

    def _makeExposure(self):
        return afwImage.ExposureF(geom.Extent2I(4, 4))

    def testFullRampWithUTR(self):
        exp = self._makeExposure()
        _stampReadRangeMetadata(exp, r0=0, r1=39, nTotal=40, applyUTRWeights=True)
        md = exp.getMetadata()
        self.assertEqual(int(md.getScalar('H4READ0')), 0)
        self.assertEqual(int(md.getScalar('H4READ1')), 39)
        self.assertEqual(int(md.getScalar('H4NREAD')), 40)
        self.assertEqual(int(md.getScalar('H4NTOT')), 40)
        self.assertTrue(bool(md.getScalar('H4UTRWT')))

    def testPartialRampMatchesFirstHalfRange(self):
        """The first-half slice carries the trimmed range, not the full ramp's."""
        exp = self._makeExposure()
        _stampReadRangeMetadata(exp, r0=1, r1=22, nTotal=40, applyUTRWeights=True)
        md = exp.getMetadata()
        self.assertEqual(int(md.getScalar('H4READ0')), 1)
        self.assertEqual(int(md.getScalar('H4READ1')), 22)
        self.assertEqual(int(md.getScalar('H4NREAD')), 22)
        self.assertEqual(int(md.getScalar('H4NTOT')), 40)

    def testCdsLayout(self):
        """CDS / non-UTR-weighted products mark themselves via H4UTRWT=False."""
        exp = self._makeExposure()
        _stampReadRangeMetadata(exp, r0=0, r1=39, nTotal=40, applyUTRWeights=False)
        self.assertFalse(bool(exp.getMetadata().getScalar('H4UTRWT')))

    def testCommentsArePresent(self):
        """Each key carries a one-line comment for the FITS header."""
        exp = self._makeExposure()
        _stampReadRangeMetadata(exp, r0=0, r1=39, nTotal=40, applyUTRWeights=True)
        md = exp.getMetadata()
        for key in ('H4READ0', 'H4READ1', 'H4NREAD', 'H4NTOT', 'H4UTRWT'):
            self.assertTrue(md.getComment(key), msg=f'{key} should carry a comment')


class StampBorderMaskTestCase(lsst.utils.tests.TestCase):
    """`_stampBorderMask` registers a BORDER mask plane and stamps
    BORDER + BAD on the outer ``borderWidth``-pixel ring of an exposure.
    """

    def _makeExposure(self, H=16, W=16):
        return afwImage.ExposureF(geom.Extent2I(W, H))

    def testBorderPlaneRegistered(self):
        exp = self._makeExposure()
        # Plane may or may not pre-exist depending on previous tests; the
        # helper must guarantee it's present afterwards regardless.
        _stampBorderMask(exp)
        self.assertIn('BORDER', exp.mask.getMaskPlaneDict())

    def testOuterFourRingHasBorderAndBad(self):
        exp = self._makeExposure(H=20, W=24)
        _stampBorderMask(exp)
        borderBit = exp.mask.getPlaneBitMask('BORDER')
        badBit = exp.mask.getPlaneBitMask('BAD')
        m = exp.mask.array
        # All four outer-4 strips carry both bits.
        for region in (m[:4, :], m[-4:, :], m[:, :4], m[:, -4:]):
            self.assertTrue((region & borderBit == borderBit).all())
            self.assertTrue((region & badBit == badBit).all())

    def testInteriorUntouched(self):
        exp = self._makeExposure(H=20, W=24)
        _stampBorderMask(exp)
        borderBit = exp.mask.getPlaneBitMask('BORDER')
        badBit = exp.mask.getPlaneBitMask('BAD')
        interior = exp.mask.array[4:-4, 4:-4]
        self.assertEqual(int((interior & borderBit).sum()), 0,
                         'interior must not carry BORDER')
        self.assertEqual(int((interior & badBit).sum()), 0,
                         'interior must not carry BAD from this helper')

    def testCountMatchesExpectedRingSize(self):
        H, W = 20, 24
        exp = self._makeExposure(H=H, W=W)
        _stampBorderMask(exp)
        borderBit = exp.mask.getPlaneBitMask('BORDER')
        nBorder = int((exp.mask.array & borderBit > 0).sum())
        # Expected: total - interior = H*W - (H-8)*(W-8) for width=4 ring.
        expected = H * W - (H - 8) * (W - 8)
        self.assertEqual(nBorder, expected)

    def testRespectsBorderWidthArgument(self):
        exp = self._makeExposure(H=16, W=16)
        _stampBorderMask(exp, borderWidth=2)
        borderBit = exp.mask.getPlaneBitMask('BORDER')
        m = exp.mask.array
        # 2-pixel ring only.
        self.assertTrue((m[:2, :] & borderBit == borderBit).all())
        self.assertEqual(int((m[2:-2, 2:-2] & borderBit).sum()), 0)

    def testIdempotent(self):
        """Calling twice doesn't add extra bits or duplicate the plane."""
        exp = self._makeExposure()
        _stampBorderMask(exp)
        snapshot = exp.mask.array.copy()
        _stampBorderMask(exp)
        np.testing.assert_array_equal(exp.mask.array, snapshot)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
