import sys
import unittest

import numpy as np

from lsst.obs.pfs.nirBadRefPixels import NirBadRefPixels
import lsst.utils.tests


class NirBadRefPixelsTestCase(lsst.utils.tests.TestCase):
    """Tests the functionality of NirBadRefPixels"""

    def checkRoundTrip(self, pixels, detectorName):
        calib = NirBadRefPixels.fromList(pixels, detectorName)
        self.assertEqual(calib.pixels.dtype, np.int32)
        self.assertFloatsEqual(calib.pixels, np.array(pixels, dtype=np.int32))
        if detectorName is not None:
            self.assertEqual(calib.metadata.get("DETNAME"), detectorName)

        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            calib.writeFits(filename)
            new = NirBadRefPixels.readFits(filename)
            self.assertFloatsEqual(new.pixels, np.array(pixels, dtype=np.int32))
            if detectorName is not None:
                self.assertEqual(new.metadata.get("DETNAME"), detectorName)

    def testBasic(self):
        """Round-trip a non-empty pixel list"""
        self.checkRoundTrip([230, 301, 342, 4033], "n1")

    def testEmpty(self):
        """Round-trip an empty pixel list (e.g. n2)"""
        self.checkRoundTrip([], "n2")


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
