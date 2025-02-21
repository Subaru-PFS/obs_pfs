import sys
from typing import TYPE_CHECKING
import unittest

import numpy as np

from lsst.daf.base import PropertyList
from lsst.afw.image import ImageF, makeExposure, makeMaskedImage
from lsst.obs.pfs.dark import PfsDark
from lsst.obs.pfs.imageCube import ImageCube
import lsst.utils.tests

if TYPE_CHECKING:
    from lsst.afw.image import ExposureF


class PfsDarkTestCase(lsst.utils.tests.TestCase):
    """Tests the functionality of PfsDark"""
    def setUp(self):
        self.metadata = PropertyList()
        self.metadata.set("FOO", "BAR")
        self.metadata.set("BAZ", 42)
        self.dimensions = (10, 10)
        self.exposureValue = 12345.67
        self.cubeSize = 3

    def makeExposure(self, isNir: bool = False) -> "ExposureF":
        """Make an ExposureF for testing"""
        self.metadata.set("W_ARM", 3 if isNir else 1)
        exposure = makeExposure(
            makeMaskedImage(ImageF(np.full(self.dimensions, self.exposureValue, dtype=np.float32)))
        )
        exposure.setMetadata(self.metadata)
        return exposure

    def makeImageCube(self) -> ImageCube:
        """Make an ImageCube for testing"""
        self.metadata.set("W_ARM", 3)  # NIR arm
        cube = ImageCube.empty(self.metadata)
        for ii in range(self.cubeSize):
            cube[ii] = ImageF(np.full(self.dimensions, float(ii), dtype=np.float32))
        return cube

    def assertMetadata(self, metadata: "PropertyList"):
        """Assert that the metadata is as expected"""
        self.assertEqual(metadata.get("FOO"), self.metadata.get("FOO"))
        self.assertEqual(metadata.get("BAZ"), self.metadata.get("BAZ"))

    def assertExposure(self, exposure: "ExposureF"):
        """Assert that the ExposureF is as expected"""
        self.assertMetadata(exposure.getMetadata())
        self.assertFloatsEqual(exposure.image.array, self.exposureValue)

    def assertImageCube(self, cube: ImageCube):
        """Assert that the ImageCube is as expected"""
        self.assertMetadata(cube.metadata)
        for ii in range(self.cubeSize):
            self.assertFloatsEqual(cube[ii].array, ii)

    def testFromFile(self):
        """Test PfsDark.fromFile"""
        exposure = self.makeExposure()
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            exposure.writeFits(filename)
            dark = PfsDark.fromFile(filename)
            self.assertEqual(dark.path, filename)
            self.assertMetadata(dark.metadata)
            self.assertExposure(dark.getCcdDark())
            self.assertRaises(RuntimeError, dark.getNirDark)

    def testFromExposure(self):
        """Test PfsDark.fromExposure"""
        exposure = self.makeExposure()
        dark = PfsDark.fromExposure(exposure)
        self.assertMetadata(dark.metadata)
        self.assertExposure(dark.getCcdDark())
        self.assertRaises(RuntimeError, dark.getNirDark)

        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            dark.writeFits(filename)
            new = PfsDark.fromFile(filename)
            self.assertEqual(new.path, filename)
            self.assertMetadata(new.metadata)
            self.assertExposure(new.getCcdDark())
            self.assertRaises(RuntimeError, new.getNirDark)

    def testFromImageCube(self):
        """Test PfsDark.fromImageCube"""
        cube = self.makeImageCube()
        dark = PfsDark.fromImageCube(cube)
        self.assertImageCube(dark.getNirDark())
        self.assertRaises(RuntimeError, dark.getCcdDark)

    def testNirExposure(self):
        """Test PfsDark.nirExposure"""
        exposure = self.makeExposure(True)  # NIR exposure: for backward compatibility
        dark = PfsDark.fromExposure(exposure)
        self.assertRaises(RuntimeError, dark.getCcdDark)  # It's expecting an ImageCube
        self.assertExposure(dark.getCcdDark(True))  # Force it to return an ExposureF
        self.assertRaises(RuntimeError, dark.getNirDark)  # There is no NIR data


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
