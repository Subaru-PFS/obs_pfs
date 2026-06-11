import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs.h4Linearity.rateStabilityDiagnostic import (
    _paginatePixels,
    _rejectedCoordsByFraction,
)


class PaginatePixelsTestCase(lsst.utils.tests.TestCase):
    """Unit tests for the pure pixel-cap / pagination helper."""

    def testSplitsIntoPages(self):
        coords = [(0, i) for i in range(26)]
        pages, nShown, nTotal = _paginatePixels(
            coords, maxPixels=64, panelsPerPage=16)
        self.assertEqual(nTotal, 26)
        self.assertEqual(nShown, 26)
        self.assertEqual([len(p) for p in pages], [16, 10])

    def testCapsAtMaxPixels(self):
        coords = [(0, i) for i in range(100)]
        pages, nShown, nTotal = _paginatePixels(
            coords, maxPixels=64, panelsPerPage=16)
        self.assertEqual(nTotal, 100)
        self.assertEqual(nShown, 64)
        self.assertEqual([len(p) for p in pages], [16, 16, 16, 16])

    def testEmpty(self):
        pages, nShown, nTotal = _paginatePixels(
            [], maxPixels=64, panelsPerPage=16)
        self.assertEqual((pages, nShown, nTotal), ([], 0, 0))

    def testRejectsBadArgs(self):
        with self.assertRaises(ValueError):
            _paginatePixels([(0, 0)], maxPixels=64, panelsPerPage=0)
        with self.assertRaises(ValueError):
            _paginatePixels([(0, 0)], maxPixels=-1, panelsPerPage=16)


class RejectedCoordsByFractionTestCase(lsst.utils.tests.TestCase):
    """Unit tests for the fraction-descending pixel ordering."""

    def testOrdersByFractionDescending(self):
        rejectMask = np.zeros((4, 4), dtype=bool)
        fraction = np.full((4, 4), np.nan, dtype=np.float32)
        rejectMask[0, 1] = True
        fraction[0, 1] = 0.30
        rejectMask[2, 3] = True
        fraction[2, 3] = 0.90
        rejectMask[3, 0] = True
        fraction[3, 0] = 0.55
        coords = _rejectedCoordsByFraction(rejectMask, fraction)
        self.assertEqual(coords, [(2, 3), (3, 0), (0, 1)])

    def testEmpty(self):
        rejectMask = np.zeros((4, 4), dtype=bool)
        fraction = np.full((4, 4), np.nan, dtype=np.float32)
        self.assertEqual(_rejectedCoordsByFraction(rejectMask, fraction), [])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
