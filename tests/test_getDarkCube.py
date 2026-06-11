"""Tests for ``PfsIsrTask.getDarkCube`` axis order and read-indexing.

Uses a fake ``ImageCube`` whose ``getImageCube`` returns a
position-encoded ``(N, H, W)`` cube where each sample's value is
``k * 10_000 + y * 100 + x``. Any axis swap or off-by-one in the
read index produces an observably wrong sample value at every
returned position.
"""
import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs import isrTask as pfsIsrTask


N_READS = 12
H = 4
W = 7


K_FACTOR = 10_000
Y_FACTOR = 100


def _encodedCube(N, H, W, dtype=np.float32):
    k = np.arange(N, dtype=dtype)[:, None, None]
    y = np.arange(H, dtype=dtype)[None, :, None]
    x = np.arange(W, dtype=dtype)[None, None, :]
    return k * K_FACTOR + y * Y_FACTOR + x


class _FakeImageCube:
    """Minimum surface PfsIsrTask.getDarkCube needs:

    - ``nreads``: total reads in the cube.
    - ``getImageCube(nreads=k)``: returns ``cube[:k]`` (a fresh array).
    - ``metadata``: mapping with optional ``"GAIN"``.
    """

    def __init__(self, cube, *, gain=1.0):
        self._cube = cube
        self.nreads = cube.shape[0]
        self.metadata = {"GAIN": gain} if gain is not None else {}

    def getImageCube(self, nreads=None):
        n = self.nreads if nreads is None else nreads
        return self._cube[:n].copy()


def _makeIsrTask():
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doFlat = False
    config.doDark = True
    config.doDefect = True
    config.doSaturationInterpolation = False
    config.h4.quickCDS = False
    config.h4.doIPC = False
    config.h4.doWriteRawCube = True
    config.h4.doLinearize = False
    config.h4.doCR = False
    config.validate()
    return pfsIsrTask.PfsIsrTask(config=config)


class GetDarkCubeTestCase(lsst.utils.tests.TestCase):
    """The returned slice's last axis corresponds to the requested
    reads; axes 0 and 1 are the (H, W) spatial plane in that order
    (time axis last, matching the H4 ISR cube convention)."""

    def setUp(self):
        self.task = _makeIsrTask()
        self.encoded = _encodedCube(N_READS, H, W)
        # The fake ImageCube still serves (N, H, W) — getDarkCube
        # transposes at the boundary.
        self.imageCube = _FakeImageCube(self.encoded)
        # Reference slice in (H, W, N) form to compare against.
        self.encodedHWN = np.ascontiguousarray(self.encoded.transpose(1, 2, 0))

    def testFullCubeRoundTrip(self):
        out = self.task.getDarkCube(self.imageCube)
        self.assertEqual(out.shape, (H, W, N_READS))
        np.testing.assert_array_equal(out, self.encodedHWN)

    def testNreadsTruncatesFromLastAxis(self):
        # r0=0 (default), nreads=k → returned[:, :, :k] = encoded[:k] transposed.
        out = self.task.getDarkCube(self.imageCube, nreads=5)
        self.assertEqual(out.shape, (H, W, 5))
        np.testing.assert_array_equal(out, self.encodedHWN[:, :, :5])

    def testR0OffsetReturnsSliceFromR0(self):
        # r0=3, nreads=5 → encoded[3:8] transposed.
        out = self.task.getDarkCube(self.imageCube, nreads=5, r0=3)
        self.assertEqual(out.shape, (H, W, 5))
        np.testing.assert_array_equal(out, self.encodedHWN[:, :, 3:8])

    def testR0OffsetWithoutNreadsReturnsToEnd(self):
        # r0=4, nreads=None → encoded[4:nreads_total] transposed.
        out = self.task.getDarkCube(self.imageCube, r0=4)
        self.assertEqual(out.shape, (H, W, N_READS - 4))
        np.testing.assert_array_equal(out, self.encodedHWN[:, :, 4:])

    def testPositionEncodedSpotChecks(self):
        # Pin specific samples so an axis swap shows up as a value error
        # rather than a shape error.
        out = self.task.getDarkCube(self.imageCube, nreads=8, r0=2)
        # Returned out[..., 0] is encoded[2]; out[y, x, k] == (k+2)*K + y*Y + x.
        self.assertEqual(
            float(out[0, 0, 0]), 2 * K_FACTOR,
        )
        self.assertEqual(
            float(out[2, 5, 3]),
            (3 + 2) * K_FACTOR + 2 * Y_FACTOR + 5,
        )
        self.assertEqual(
            float(out[3, 6, 7]),
            (7 + 2) * K_FACTOR + 3 * Y_FACTOR + 6,
        )

    def testGainBackoutDivides(self):
        # If the dark cube was stored as electrons (gain applied), the
        # method divides by ``GAIN`` to back out the conversion.
        gain = 2.5
        imageCube = _FakeImageCube(self.encoded.copy(), gain=gain)
        out = self.task.getDarkCube(imageCube)
        np.testing.assert_allclose(out, self.encodedHWN / gain, atol=1e-3)

    def testGainOneNoChange(self):
        out = self.task.getDarkCube(self.imageCube)
        np.testing.assert_array_equal(out, self.encodedHWN)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
