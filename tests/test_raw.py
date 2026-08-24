"""``lsst.obs.pfs.raw.rotateImageBy90Striped`` must reproduce ``rotateImageBy90``.

Every H4 read is rotated on the way in, so the rotation runs 280 times per
140-read quantum on 67 MB planes. The afw whole-image rotation is
memory-latency-bound: each output column touches a separate source cache line,
so the working set is the entire image. The striped form is a cache-blocking
rewrite and must be exactly equal, not merely close.

The IRP4 reference plane is non-square -- ``(4096, 1024)`` on disk, ``(1024,
4096)`` after rotation -- so an H/W swap would survive a square-only test.
"""
import unittest

import numpy as np

import lsst.afw.math as afwMath
import lsst.utils.tests
from lsst.afw.image import ImageF

import lsst.obs.pfs.raw as pfsRawModule
from lsst.obs.pfs.raw import PfsRaw, rotateImageBy90Striped


def _afwReference(array, nQuarter):
    """The rotation as afw does it on the whole image at once."""
    return afwMath.rotateImageBy90(ImageF(array.copy()), nQuarter).getArray()


def _striped(array, nQuarter, **kwargs):
    return rotateImageBy90Striped(ImageF(array.copy()), nQuarter, **kwargs)


class RotateImageBy90StripedTestCase(lsst.utils.tests.TestCase):
    # (height, width) before rotation. The last is the IRP4 REF aspect ratio.
    SHAPES = ((1, 1), (1, 7), (7, 1), (4, 4), (5, 8), (8, 5), (64, 16), (32, 8))

    def testMatchesAfwForEveryShapeAndQuarter(self):
        rng = np.random.RandomState(1885)
        for shape in self.SHAPES:
            array = rng.uniform(-1e4, 1e4, size=shape).astype(np.float32)
            for nQuarter in range(4):
                with self.subTest(shape=shape, nQuarter=nQuarter):
                    expected = _afwReference(array, nQuarter)
                    got = _striped(array, nQuarter)
                    self.assertEqual(got.shape, expected.shape)
                    self.assertEqual(got.dtype, array.dtype)
                    np.testing.assert_array_equal(got, expected)

    def testProductionQuarterIsThree(self):
        # The H4 detectors carry yaw=270, i.e. nQuarter=3; that is the only
        # rotation the ramp path actually exercises.
        rng = np.random.RandomState(3)
        array = rng.uniform(0.0, 6e4, size=(96, 24)).astype(np.float32)
        np.testing.assert_array_equal(_striped(array, 3),
                                      _afwReference(array, 3))

    def testStripCountDoesNotChangeResult(self):
        # Strip count is a pure cache-blocking knob. Counts that do not
        # divide the image evenly must give the same answer as ones that do.
        rng = np.random.RandomState(17)
        array = rng.uniform(-1e3, 1e3, size=(48, 30)).astype(np.float32)
        for nQuarter in range(4):
            expected = _afwReference(array, nQuarter)
            for nStrips in (1, 2, 7, 16, 30, 48, 1000):
                with self.subTest(nQuarter=nQuarter, nStrips=nStrips):
                    np.testing.assert_array_equal(
                        _striped(array, nQuarter, nStrips=nStrips),
                        expected)

    def testResultIsContiguousAndIndependent(self):
        # The caller wraps the result in an ImageF and later mutates it
        # (replaceNansWith0, dark subtraction), so it must own its memory.
        array = np.arange(24, dtype=np.float32).reshape(4, 6)
        for nQuarter in range(4):
            got = _striped(array, nQuarter)
            self.assertTrue(got.flags["C_CONTIGUOUS"])
            got[...] = -1.0
            self.assertFalse((array == -1.0).any())

    def testNegativeAndLargeQuartersWrap(self):
        rng = np.random.RandomState(5)
        array = rng.uniform(-10.0, 10.0, size=(6, 10)).astype(np.float32)
        for nQuarter in range(4):
            expected = _striped(array, nQuarter)
            for equivalent in (nQuarter - 4, nQuarter + 4, nQuarter + 8):
                np.testing.assert_array_equal(
                    _striped(array, equivalent), expected)


class DetectorCachingTestCase(lsst.utils.tests.TestCase):
    """``PfsRaw.detector`` must be built once per instance.

    Every ramp read asks for the detector (to get its rotation), so an
    uncached property reloads the whole camera 280 times per 140-read quantum.
    ``__init__`` already allocates ``_detector`` for this.
    """

    def setUp(self):
        self.calls = []
        self.sentinel = object()

        def fakeLoadCamera(pfsCategory):
            self.calls.append(pfsCategory)
            return {7: self.sentinel}

        self.original = pfsRawModule.loadCamera
        pfsRawModule.loadCamera = fakeLoadCamera

    def tearDown(self):
        pfsRawModule.loadCamera = self.original

    def _makeRaw(self):
        # A CCD arm (W_ARM != 3) so the NIR gain rebuild is skipped; the
        # caching under test is common to both arms.
        raw = PfsRaw.__new__(PfsRaw)
        raw.path = "unused"
        raw.pfsCategory = None
        raw._metadata = {"DET-ID": 7, "W_ARM": 1}
        raw._detector = None
        raw._obsInfo = None
        raw._visitInfo = None
        return raw

    def testBuiltOnceAndReturnedByIdentity(self):
        raw = self._makeRaw()
        first = raw.detector
        self.assertIs(first, self.sentinel)
        for _ in range(9):
            self.assertIs(raw.detector, first)
        self.assertEqual(len(self.calls), 1,
                         f"camera loaded {len(self.calls)} times, expected 1")

    def testCacheIsPerInstance(self):
        a, b = self._makeRaw(), self._makeRaw()
        self.assertIs(a.detector, b.detector)   # same sentinel from the camera
        self.assertEqual(len(self.calls), 2,
                         "each PfsRaw must build its own detector")


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
