import sys
import unittest

import numpy as np
from numpy.typing import ArrayLike
import lsst.utils.tests

from lsst.obs.pfs.blackSpots import BlackSpotsConfig


class BlackSpotsTestCase(lsst.utils.tests.TestCase):
    """Tests the reading and using of black spots"""

    def setUp(self):
        self.spots = BlackSpotsConfig().read()
        self.rng = np.random.RandomState(12345)

    def checkNearest(self, x: ArrayLike, y: ArrayLike):
        """Checks that the result of findNearest is correct

        Parameters
        ----------
        x : array_like
            x-coordinates of point (mm).
        y : array_like
            y-coordinates of point (mm).
        """
        result = self.spots.findNearest(x, y)

        if np.isscalar(x):
            distance = np.hypot(self.spots.x - x, self.spots.y - y)
            best = np.argmin(distance)
            if not np.isfinite(x) or not np.isfinite(y):
                self.assertFalse(np.isfinite(result.distance))
                self.assertEqual(result.spotId, -1)
                self.assertTrue(np.isnan(result.x))
                self.assertTrue(np.isnan(result.y))
                self.assertTrue(np.isnan(result.r))
            else:
                self.assertFloatsAlmostEqual(result.distance, distance[best], atol=2.0e-16)
                self.assertEqual(result.spotId, self.spots.spotId[best])
                self.assertEqual(result.x, self.spots.x[best])
                self.assertEqual(result.y, self.spots.y[best])
                self.assertFloatsEqual(result.r, self.spots.r[best])
            return

        for ii in range(len(x)):
            distance = np.hypot(self.spots.x - x[ii], self.spots.y - y[ii])
            best = np.argmin(distance)
            if not np.isfinite(x[ii]) or not np.isfinite(y[ii]):
                self.assertFalse(np.isfinite(result.distance[ii]))
                self.assertEqual(result.spotId[ii], -1)
                self.assertTrue(np.isnan(result.x[ii]))
                self.assertTrue(np.isnan(result.y[ii]))
                self.assertTrue(np.isnan(result.r[ii]))
            else:
                self.assertFloatsAlmostEqual(result.distance[ii], distance[best], atol=2.0e-16)
                self.assertEqual(result.spotId[ii], self.spots.spotId[best])
                self.assertEqual(result.x[ii], self.spots.x[best])
                self.assertEqual(result.y[ii], self.spots.y[best])
                self.assertFloatsEqual(result.r[ii], self.spots.r[best])

    def testScalar(self):
        """Test scalar inputs"""
        self.checkNearest(0, 0)
        self.checkNearest(100, -100)
        self.checkNearest(0, np.nan)
        self.checkNearest(np.nan, 0)
        self.checkNearest(np.nan, np.nan)

    def testArray(self):
        """Test array inputs"""
        xy = self.rng.uniform(low=-200, high=200, size=(100, 2))
        xy[3, 0] = np.nan
        xy[4, 1] = np.nan
        xy[5, :] = np.nan
        self.checkNearest(xy[:, 0], xy[:, 1])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
