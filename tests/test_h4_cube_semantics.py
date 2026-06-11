"""Sanity tests for the (N, H, W) ramp cube axis convention.

These tests are a canary for any change that accidentally swaps axes,
permutes a layout, or otherwise scrambles the convention used across
the H4 ISR pipeline: axis 0 is the read index, axes 1 and 2 are the
spatial (y, x) plane. They use a value-encoded cube where each sample
encodes its own ``(k, y, x)`` position, so any axis permutation yields
a measurably different (and diagnosable) value at any sample.

Tests in other files import ``valueEncodedRamp`` from here as a
fixture for integration tests that need to track per-sample identity
across pipeline stages.
"""
import unittest

import numpy as np

import lsst.utils.tests


K_FACTOR = 10_000
Y_FACTOR = 100
# x has factor 1 (no separate constant). With H <= 99 and W <= 99,
# the contribution from each axis is unambiguous: the floor-divide by
# 10_000 recovers k, by 100 (after subtracting k * K_FACTOR) recovers y,
# and the remainder is x. Any axis swap or permutation in pipeline code
# yields a value with the wrong axis's contribution in the wrong place.


def valueEncodedRamp(nReads, H, W, *, dtype=np.float32):
    """Build a cube where ``cube[k, y, x] = k*K_FACTOR + y*Y_FACTOR + x``.

    Returns a C-contiguous ``(nReads, H, W)`` array. Designed for tests
    that need to verify axis order is preserved across a transformation:
    if a function accidentally permutes axes, the value at any returned
    sample will not match the expected encoding.
    """
    k = np.arange(nReads, dtype=dtype)[:, None, None]
    y = np.arange(H, dtype=dtype)[None, :, None]
    x = np.arange(W, dtype=dtype)[None, None, :]
    return k * K_FACTOR + y * Y_FACTOR + x


class CubeSemanticsTestCase(lsst.utils.tests.TestCase):
    """Pin down what ``cube[k, y, x]`` means under the project convention."""

    def setUp(self):
        # Asymmetric N != H != W so axis swaps would change the shape.
        self.N = 11
        self.H = 5
        self.W = 7
        self.cube = valueEncodedRamp(self.N, self.H, self.W)

    def testShape(self):
        self.assertEqual(self.cube.shape, (self.N, self.H, self.W))

    def testEncoding(self):
        # cube[k, y, x] == k * K + y * Y + x. Spot-check asymmetric coords
        # so a (y, x) transposition would be loud.
        self.assertEqual(self.cube[5, 2, 3], 5 * K_FACTOR + 2 * Y_FACTOR + 3)
        self.assertEqual(self.cube[3, 4, 6], 3 * K_FACTOR + 4 * Y_FACTOR + 6)
        self.assertEqual(self.cube[10, 0, 0], 10 * K_FACTOR)
        self.assertEqual(self.cube[0, 4, 6], 4 * Y_FACTOR + 6)

    def testPerPixelRamp(self):
        # cube[:, y, x] is the per-pixel time series at (y, x).
        ramp = self.cube[:, 2, 3]
        self.assertEqual(ramp.shape, (self.N,))
        # Steps along the ramp: cube[k+1, y, x] - cube[k, y, x] = K_FACTOR.
        np.testing.assert_array_equal(
            np.diff(ramp), K_FACTOR * np.ones(self.N - 1, dtype=ramp.dtype)
        )

    def testPerReadFrame(self):
        # cube[k, :, :] is the 2-D frame for read k.
        frame = self.cube[5, :, :]
        self.assertEqual(frame.shape, (self.H, self.W))
        for y in range(self.H):
            for x in range(self.W):
                self.assertEqual(frame[y, x], 5 * K_FACTOR + y * Y_FACTOR + x)

    def testDiffAxis0(self):
        # np.diff(axis=0) computes per-read deltas; shape (N-1, H, W).
        deltas = np.diff(self.cube, axis=0)
        self.assertEqual(deltas.shape, (self.N - 1, self.H, self.W))
        # Every delta equals K_FACTOR (the per-read step).
        np.testing.assert_array_equal(
            deltas, K_FACTOR * np.ones_like(deltas)
        )

    def testCumsumAxis0Roundtrip(self):
        # diff -> cumsum reconstructs the cube modulo the initial frame.
        deltas = np.diff(self.cube, axis=0)
        reconstructed = np.empty_like(self.cube)
        reconstructed[0:1] = self.cube[0:1]
        np.cumsum(deltas, axis=0, out=reconstructed[1:])
        reconstructed[1:] += self.cube[0:1]
        np.testing.assert_array_equal(reconstructed, self.cube)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
