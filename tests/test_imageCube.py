import sys
import unittest
import fitsio
import numpy as np

from lsst.daf.base import PropertyList
from lsst.afw.image import ImageF
from lsst.obs.pfs.imageCube import ImageCube
import lsst.utils.tests

# Pixels are quantized on write, to within half a step.
STEP = abs(ImageCube.compression["qlevel"])
ATOL = 0.5*STEP + 1e-3


class ImageCubeTestCase(lsst.utils.tests.TestCase):
    """Tests the functionality of ImageCube"""
    def testBasic(self):
        """Test basic functionality"""
        numImages = 3

        metadata = PropertyList()
        metadata.set("FOO", "BAR")
        metadata.set("BAZ", 42)
        dimensions = (10, 10)

        cube = ImageCube.empty(metadata)
        for ii in range(numImages):
            cube[ii] = ImageF(np.full(dimensions, float(ii), dtype=np.float32))

        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            cube.writeFits(filename)

            with ImageCube.fromFile(filename) as new:
                for name in metadata.names():
                    self.assertEqual(new.metadata.get(name), metadata.get(name))
                for ii in range(numImages):
                    image = new[ii]
                    self.assertIsInstance(image, ImageF)
                    self.assertEqual(image.array.shape, dimensions)
                    self.assertFloatsAlmostEqual(image.array, ii, atol=ATOL, rtol=0)


class ImageCubeReadTestCase(lsst.utils.tests.TestCase):
    """The read path must not retain frames it was asked not to cache.

    ``getReadArray`` promises not to cache: a 139-read dark cube is read one
    frame at a time while the ramp cube is also live, and holding every frame
    would add 11.4 GB.
    """

    NUM = 6
    DIMS = (8, 5)

    def _write(self, filename):
        # fromCube is how the real nirDark cubes are built, and is what
        # stamps NREADS.
        metadata = PropertyList()
        metadata.set("GAIN", 1.0)
        data = np.stack([self._expected(ii) for ii in range(self.NUM)])
        ImageCube.fromCube(data, metadata).writeFits(filename)

    def _expected(self, ii):
        plane = np.full(self.DIMS, float(ii), dtype=np.float32)
        plane[0, 0] = 100.0 + ii
        return plane

    def testValuesRoundTrip(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                for ii in range(self.NUM):
                    np.testing.assert_allclose(cube.getReadArray(ii),
                                               self._expected(ii), atol=ATOL, rtol=0)
                    self.assertEqual(cube.getReadArray(ii).dtype, np.float32)

    def testGetReadArrayDoesNotCache(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                for ii in range(self.NUM):
                    cube.getReadArray(ii)
                self.assertEqual(len(cube._images), 0,
                                 "getReadArray must not populate the cache")

    def testGetReadArrayReturnsIndependentArrays(self):
        # The dark-subtract path scales the returned frame, so successive
        # calls must not hand back the same buffer.
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                first = cube.getReadArray(2)
                second = cube.getReadArray(2)
                self.assertFalse(np.shares_memory(first, second))
                first[...] = -999.0
                np.testing.assert_allclose(cube.getReadArray(2),
                                           self._expected(2), atol=ATOL, rtol=0)

    def testGetItemStillCaches(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                image = cube[3]
                self.assertIsInstance(image, ImageF)
                self.assertIs(cube[3], image)
                np.testing.assert_allclose(image.array, self._expected(3), atol=ATOL, rtol=0)
                # A cached read is what getReadArray hands back thereafter.
                self.assertTrue(np.shares_memory(cube.getReadArray(3),
                                                 image.array))

    def testGetImageCubeMatches(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                stack = cube.getImageCube()
                self.assertEqual(stack.shape, (self.NUM, *self.DIMS))
                for ii in range(self.NUM):
                    np.testing.assert_allclose(stack[ii], self._expected(ii), atol=ATOL, rtol=0)
                stack = cube.getImageCube(nreads=2)
                self.assertEqual(stack.shape, (2, *self.DIMS))

    def testReadAllCaches(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                cube.readAll()
                self.assertEqual(len(cube._images), self.NUM)
                for ii in range(self.NUM):
                    np.testing.assert_allclose(cube[ii].array,
                                               self._expected(ii), atol=ATOL, rtol=0)

    def testFlushClearsCacheAndKeepsReadCount(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                cube.readAll()
                nreads = cube.getNumReads()
                cube.flush()
                self.assertEqual(len(cube._images), 0)
                self.assertEqual(cube.getNumReads(), nreads,
                                 "flushing the cache must not lose the read count")


class ImageCubeCompressionTestCase(lsst.utils.tests.TestCase):
    """Cubes are written with RICE and a fixed 0.1 e- quantization step.

    cfitsio's default float quantization picks the step per tile from a noise
    estimate, which images with vertical spectra break: the step grows with
    the signal. A fixed step, subtractively dithered, bounds the error at half
    a step everywhere and leaves it unbiased.
    """

    NUM = 4
    DIMS = (37, 29)

    def setUp(self):
        rng = np.random.default_rng(1889)
        self.data = rng.normal(100.0, 7.3, size=(self.NUM, *self.DIMS)).astype(np.float32)
        self.metadata = PropertyList()
        self.metadata.set("GAIN", 1.7416)

    def assertQuantized(self, filename, data=None):
        data = self.data if data is None else data
        with ImageCube.fromFile(filename) as cube:
            self.assertEqual(cube.getNumReads(), len(data))
            err = np.array([cube.getReadArray(ii) for ii in range(len(data))], dtype=np.float64) - data
        self.assertLessEqual(np.abs(err).max(), ATOL)
        self.assertLess(abs(err.mean()), 0.003)
        self.assertFloatsAlmostEqual(err.std(), STEP/np.sqrt(12), rtol=0.1)
        with fitsio.FITS(filename) as fits:
            for ii in range(len(data)):
                header = fits[f"IMAGE_{ii + 1}"].read_header()
                # cfitsio writes the RICE_1 synonym RICE_ONE
                self.assertIn(header.get("ZCMPTYPE"), ("RICE_1", "RICE_ONE"))
                self.assertEqual(header.get("ZQUANTIZ"), "SUBTRACTIVE_DITHER_2")

    def testWriteFits(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            ImageCube.fromCube(self.data, self.metadata).writeFits(filename)
            self.assertQuantized(filename)

    def testWriteCubeData(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            ImageCube.writeCubeData(filename, self.data, self.metadata)
            self.assertQuantized(filename)

    def testWriteFitsOverwrites(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            ImageCube.fromCube(self.data[:2] + 5, self.metadata).writeFits(filename)
            ImageCube.fromCube(self.data, self.metadata).writeFits(filename)
            self.assertQuantized(filename)

    def testStepIsIndependentOfSignal(self):
        """Bright vertical traces do not coarsen the quantization."""
        x = np.arange(self.DIMS[1])
        profile = sum(np.exp(-0.5*((x - xc)/1.5)**2) for xc in (5.3, 14.3, 23.3))
        data = (self.data + 2e4*profile[None, None, :]).astype(np.float32)
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            ImageCube.fromCube(data, self.metadata).writeFits(filename)
            self.assertQuantized(filename, data)

    def testZerosPreserved(self):
        data = self.data.copy()
        data[:, :5, :] = 0.0
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            ImageCube.fromCube(data, self.metadata).writeFits(filename)
            with ImageCube.fromFile(filename) as cube:
                for ii in range(self.NUM):
                    np.testing.assert_array_equal(cube.getReadArray(ii)[:5, :], 0.0)

    def testReproducible(self):
        """The dither is seeded from the data, so a rewrite gives the same values."""
        arrays = []
        for _ in range(2):
            with lsst.utils.tests.getTempFilePath(".fits") as filename:
                ImageCube.fromCube(self.data, self.metadata).writeFits(filename)
                with ImageCube.fromFile(filename) as cube:
                    arrays.append(cube.getImageCube())
        np.testing.assert_array_equal(arrays[0], arrays[1])

    def testHeaderRoundTrip(self):
        metadata = PropertyList()
        metadata.set("GAIN", 1.7416)
        metadata.set("NAME", "n2")
        metadata.set("FLAG", True)
        metadata.set("COUNT", 42)
        metadata.set("W_H4IRPN", 4)
        metadata.set("CALIB_INPUT_0", 144635)  # longer than 8 characters
        metadata.set("LONGSTR", "x" * 100)
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            ImageCube.fromCube(self.data, metadata).writeFits(filename)
            with ImageCube.fromFile(filename) as cube:
                for name in metadata.names():
                    self.assertEqual(cube.metadata.get(name), metadata.get(name), name)
                self.assertEqual(cube.metadata.get("NREADS"), self.NUM)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
