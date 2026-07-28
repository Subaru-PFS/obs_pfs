"""``_utrVariance`` seeds the per-pixel variance from the signal, but with the
photon term floored at zero: a negative pixel value is read/dark noise, not a
negative photon count, so it must not produce a negative variance (which
``maskNegativeVariance`` would then flag BAD). The variance is therefore always
positive; positive pixels are unchanged.
"""
import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs.isrTask import _utrVariance


def _utrCoeffs(nread):
    A = 6 * (nread * nread + 1) / (5 * nread * (nread + 1))
    floor = 12 * (nread - 1) / (nread * (nread + 1))
    return A, floor


class UtrVarianceTestCase(lsst.utils.tests.TestCase):
    def test_negativeSignalFlooredToReadNoiseFloor(self):
        # A negative pixel gets exactly the read-noise floor, never a negative
        # variance -- this is the fix that stops the spurious BAD masking.
        nread, rn = 90, 14.8
        img = np.array([[-1.0e6, -50.0, 0.0]], dtype=np.float32)
        var = _utrVariance(img, nread, rn, utrWeighted=True)
        _, floorCoeff = _utrCoeffs(nread)
        floor = floorCoeff * rn ** 2
        self.assertTrue((var > 0).all(), "variance must be strictly positive")
        np.testing.assert_allclose(var[0], floor, rtol=1e-6)

    def test_positiveSignalUnchanged(self):
        # Positive pixels keep A*signal + B*RN^2, i.e. identical to before the fix.
        nread, rn = 90, 14.8
        img = np.array([[100.0]], dtype=np.float32)
        var = _utrVariance(img, nread, rn, utrWeighted=True)
        A, floorCoeff = _utrCoeffs(nread)
        np.testing.assert_allclose(var[0, 0], A * 100.0 + floorCoeff * rn ** 2, rtol=1e-6)

    def test_cdsBranch(self):
        # CDS: var = max(signal,0) + 2*RN^2, always positive.
        rn = 14.8
        img = np.array([[-30.0, 50.0]], dtype=np.float32)
        var = _utrVariance(img, 0, rn, utrWeighted=False)
        np.testing.assert_allclose(var[0, 0], 2 * rn ** 2, rtol=1e-6)
        np.testing.assert_allclose(var[0, 1], 50.0 + 2 * rn ** 2, rtol=1e-6)
        self.assertTrue((var > 0).all())

    def test_inputNotMutated(self):
        # The input signal image must not be modified in place.
        img = np.array([[-5.0, 7.0]], dtype=np.float32)
        before = img.copy()
        _utrVariance(img, 90, 14.8, utrWeighted=True)
        np.testing.assert_array_equal(img, before)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(__import__("sys").modules["__main__"])
    unittest.main()
