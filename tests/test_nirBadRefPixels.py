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

    def testProvenanceInHeader(self):
        """Scan provenance survives the FITS round-trip as header cards."""
        provenance = dict(
            generatedBy="tester",
            generatedAt="2026-07-06T00:00:00+00:00",
            arm="n",
            visits=[144587, 144588, 144589],
            visitDates={144587: "2026-06-15", 144588: "2026-06-16"},
            threshold=2.0,
            maxThreshold=4.0,
            minVisits=3,
            nreads=None,
        )
        calib = NirBadRefPixels.fromList([230, 301], "n1", provenance=provenance)
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            calib.writeFits(filename)
            new = NirBadRefPixels.readFits(filename)
        md = new.metadata
        self.assertEqual(md.get("DETNAME"), "n1")
        self.assertEqual(md.get("CALGENBY"), "tester")
        self.assertEqual(md.get("CALVIS"), "144587,144588,144589")
        self.assertEqual(md.get("CALNVIS"), 3)
        self.assertEqual(md.get("CALDBEG"), "2026-06-15")
        self.assertEqual(md.get("CALDEND"), "2026-06-16")
        self.assertEqual(md.get("CALMINVI"), 3)
        self.assertEqual(md.get("CALNRD"), "all")


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
