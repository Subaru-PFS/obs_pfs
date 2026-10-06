"""``PfsRaw.getNirRead`` dithers the integer H4 reads by U(-0.5, 0.5).

The H4 ASIC delivers integer ADU. Medians of integer-valued reads (the
read-by-read combine of many dark ramps, the IRP reference filter) snap to the
integer lattice, which shows up as read-to-read offsets in the combined dark.
Dithering each read as it is read in removes the lattice. The dither is seeded
from the visit, detector, read and HDU type, so that a read gives the same
values however many times it is read.
"""
import os
import shutil
import tempfile
import types
import unittest

import fitsio
import numpy as np

import lsst.utils.tests
from lsst.geom import Box2I, Point2I, Extent2I

from lsst.afw.image import ImageF
from lsst.obs.pfs.raw import PfsRaw, rotateImageBy90Striped

SHAPE = (16, 12)  # (height, width) on disk
NREADS = 3


def writeRaw(path, visit=12345, detId=3):
    """Write a minimal PFSB-like NIR ramp: IMAGE_n/REF_n uint16 HDU pairs."""
    rng = np.random.default_rng(1889)
    with fitsio.FITS(path, "rw", clobber=True) as fits:
        fits.write(None, header=dict(W_ARM=3, W_4FMTVR=2, W_H4NRED=NREADS,
                                     W_VISIT=visit, **{"DET-ID": detId}))
        planes = {}
        for r in range(1, NREADS + 1):
            for kind in ("IMAGE", "REF"):
                plane = rng.integers(1000, 30000, size=SHAPE).astype(np.uint16)
                fits.write(plane, extname=f"{kind}_{r}")
                planes[kind, r] = plane
    return planes


def makeRaw(path, visit=12345, detId=3, nQuarter=0):
    raw = PfsRaw.__new__(PfsRaw)
    raw.path = path
    raw.pfsCategory = None
    raw._metadata = {"W_ARM": 3, "W_4FMTVR": 2, "W_H4NRED": NREADS,
                     "W_VISIT": visit, "DET-ID": detId}
    orientation = types.SimpleNamespace(getNQuarter=lambda: nQuarter)
    raw._detector = types.SimpleNamespace(getOrientation=lambda: orientation)
    raw._obsInfo = None
    raw._visitInfo = None
    return raw


class RawDitherTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="rawDither-")
        self.path = os.path.join(self.root, "raw.fits")
        self.planes = writeRaw(self.path)

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def testUnditheredIsExact(self):
        raw = makeRaw(self.path)
        for kind in ("IMAGE", "REF"):
            got = raw.getNirRead(2, imageType=kind, doRotate=False, doDither=False).array
            np.testing.assert_array_equal(got, self.planes[kind, 2])

    def testDitherIsWithinHalfADU(self):
        raw = makeRaw(self.path)
        for kind in ("IMAGE", "REF"):
            with self.subTest(kind=kind):
                got = raw.getNirRead(2, imageType=kind, doRotate=False).array
                delta = got - self.planes[kind, 2]
                self.assertTrue(np.all(np.abs(delta) <= 0.5))
                # Not left on the integer lattice
                self.assertGreater(np.mean(np.abs(delta) > 1e-3), 0.9)
                self.assertLess(abs(np.mean(delta)), 0.1)

    def testDitherIsReproducible(self):
        a = makeRaw(self.path).getNirRead(1, doRotate=False).array
        b = makeRaw(self.path).getNirRead(1, doRotate=False).array
        np.testing.assert_array_equal(a, b)

    def testDitherDiffersByReadTypeVisitAndDetector(self):
        def delta(raw, readNum, kind):
            return raw.getNirRead(readNum, imageType=kind, doRotate=False).array - self.planes[kind, readNum]

        base = delta(makeRaw(self.path), 1, "IMAGE")
        others = dict(
            read=delta(makeRaw(self.path), 2, "IMAGE"),
            kind=delta(makeRaw(self.path), 1, "REF"),
            visit=delta(makeRaw(self.path, visit=12346), 1, "IMAGE"),
            detector=delta(makeRaw(self.path, detId=4), 1, "IMAGE"),
        )
        for name, other in others.items():
            with self.subTest(name=name):
                self.assertFalse(np.allclose(base, other))

    def testRotationAppliesToDitheredRead(self):
        unrotated = makeRaw(self.path).getNirRead(1, doRotate=False).array
        rotated = makeRaw(self.path, nQuarter=1).getNirRead(1).array
        np.testing.assert_array_equal(rotated, rotateImageBy90Striped(ImageF(unrotated.copy()), 1))

    def testBBoxMatchesFullRead(self):
        raw = makeRaw(self.path)
        full = raw.getNirRead(1, doRotate=False).array
        bbox = Box2I(Point2I(3, 5), Extent2I(4, 6))
        sub = raw.getNirRead(1, bbox=bbox, doRotate=False).array
        np.testing.assert_array_equal(sub, full[5:11, 3:7])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
