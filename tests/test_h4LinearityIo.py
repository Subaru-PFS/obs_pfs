import os
import tempfile
import unittest

import numpy as np
from astropy.io import fits

import lsst.utils.tests
from lsst.obs.pfs.h4Linearity.io import isH4LinearityFile
from lsst.obs.pfs.h4Linearity.models import MODEL_REGISTRY


class IsH4LinearityFileTestCase(lsst.utils.tests.TestCase):
    """`isH4LinearityFile` discriminates h4Linearity files (carrying a
    valid ``MODEL`` keyword) from non-h4Linearity FITS files and from
    nonexistent / unreadable paths.
    """

    def _writeFits(self, header_pairs, dirpath):
        path = os.path.join(dirpath, "test.fits")
        hdr = fits.Header()
        for k, v in header_pairs:
            hdr[k] = v
        fits.PrimaryHDU(header=hdr).writeto(path)
        return path

    def testKnownModelReturnsTrue(self):
        knownModel = next(iter(MODEL_REGISTRY))
        with tempfile.TemporaryDirectory() as d:
            path = self._writeFits([("MODEL", knownModel)], d)
            self.assertTrue(isH4LinearityFile(path))

    def testUnknownModelReturnsFalse(self):
        with tempfile.TemporaryDirectory() as d:
            path = self._writeFits([("MODEL", "not-a-model-name-XYZ")], d)
            self.assertFalse(isH4LinearityFile(path))

    def testMissingModelKeyReturnsFalse(self):
        """A FITS file without a MODEL keyword (legacy nirLinearity style)."""
        with tempfile.TemporaryDirectory() as d:
            path = self._writeFits([("OBJECT", "anything")], d)
            self.assertFalse(isH4LinearityFile(path))

    def testNonexistentPathReturnsFalse(self):
        self.assertFalse(isH4LinearityFile("/nonexistent/path/that/does/not/exist.fits"))

    def testNonFitsFileReturnsFalse(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "garbage.fits")
            with open(path, "wb") as f:
                f.write(b"this is not a FITS file")
            self.assertFalse(isH4LinearityFile(path))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
