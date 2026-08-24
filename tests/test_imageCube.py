import sys
import unittest
import numpy as np

from lsst.daf.base import PropertyList
from lsst.afw.image import ImageF
from lsst.obs.pfs.imageCube import ImageCube
import lsst.utils.tests


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
                    self.assertFloatsEqual(image.array, ii)


class ImageCubeReadTestCase(lsst.utils.tests.TestCase):
    """The read path must not retain frames it was asked not to cache.

    ``getReadArray`` promises not to cache. Reading through
    ``astropy.io.fits`` broke that promise -- astropy caches ``hdu.data`` on
    the HDU itself -- so a 139-read dark cube, read once and used once, stayed
    resident for the rest of the quantum: 11.4 GB arriving exactly when the
    ramp cube is also live.
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
                    np.testing.assert_array_equal(cube.getReadArray(ii),
                                                  self._expected(ii))
                    self.assertEqual(cube.getReadArray(ii).dtype, np.float32)

    def testGetReadArrayDoesNotCache(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                for ii in range(self.NUM):
                    cube.getReadArray(ii)
                self.assertEqual(len(cube._images), 0,
                                 "getReadArray must not populate the cache")
                for hdu in cube.fits[1:]:
                    self.assertFalse(
                        hdu._data_loaded,
                        f"{hdu.name} data retained by the FITS reader")

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
                np.testing.assert_array_equal(cube.getReadArray(2),
                                              self._expected(2))

    def testGetItemStillCaches(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                image = cube[3]
                self.assertIsInstance(image, ImageF)
                self.assertIs(cube[3], image)
                np.testing.assert_array_equal(image.array, self._expected(3))
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
                    np.testing.assert_array_equal(stack[ii], self._expected(ii))
                stack = cube.getImageCube(nreads=2)
                self.assertEqual(stack.shape, (2, *self.DIMS))

    def testReadAllCaches(self):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self._write(filename)
            with ImageCube.fromFile(filename) as cube:
                cube.readAll()
                self.assertEqual(len(cube._images), self.NUM)
                for ii in range(self.NUM):
                    np.testing.assert_array_equal(cube[ii].array,
                                                  self._expected(ii))

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


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
