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


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
