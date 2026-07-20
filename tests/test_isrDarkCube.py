"""Tests for how ``PfsIsrTask`` reads the NIR dark cube back out.

The dark is stored in electrons and subtracted from the raw ramp in ADU, so both
accessors divide by the gain the dark was taken with. Neither may scale the
`ImageCube`'s cached images in place.
"""

import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

import lsst.utils.tests

from lsst.obs.pfs.imageCube import ImageCube

from testUtils import HAS_DRP_STELLA, requireDrpStella

if HAS_DRP_STELLA:
    # isrTask imports pfs.drp.stella.crosstalk.
    from lsst.obs.pfs.isrTask import PfsIsrTask

GAIN = 2.0
LEVEL = 10.0
NREADS = 3


@requireDrpStella
class DarkCubeAccessTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="isrDarkCube-")
        self.path = os.path.join(self.root, "nirDark.fits")
        data = np.full((NREADS, 2, 2), LEVEL, dtype="f4")
        ImageCube.fromCube(data, dict(GAIN=GAIN)).writeFits(self.path)
        # The helpers under test touch no task state, so skip __init__.
        self.task = PfsIsrTask.__new__(PfsIsrTask)

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def testGetDarkReadBacksOutGain(self):
        dark = ImageCube.readFits(self.path)
        self.assertFloatsAlmostEqual(self.task.getDarkRead(dark, 0), LEVEL / GAIN)

    def testGetDarkReadDoesNotCorruptTheCube(self):
        """Reading the same read twice must give the same answer."""
        dark = ImageCube.readFits(self.path)
        first = self.task.getDarkRead(dark, 0).copy()
        second = self.task.getDarkRead(dark, 0)
        self.assertFloatsAlmostEqual(second, first)
        # ...and the cached image is still in electrons, unscaled.
        self.assertFloatsAlmostEqual(dark[0].array, LEVEL)

    def testGetDarkReadReturnsNewArray(self):
        """The caller must not be handed the cube's cached image."""
        dark = ImageCube.readFits(self.path)
        read = self.task.getDarkRead(dark, 0)
        read += 1000.0
        self.assertFloatsAlmostEqual(dark[0].array, LEVEL)

    def testGetDarkReadWithoutGain(self):
        path = os.path.join(self.root, "noGain.fits")
        ImageCube.fromCube(np.full((NREADS, 2, 2), LEVEL, dtype="f4"), {}).writeFits(path)
        dark = ImageCube.readFits(path)
        read = self.task.getDarkRead(dark, 0)
        self.assertFloatsAlmostEqual(read, LEVEL)
        read += 1000.0  # still a copy, even when the gain is 1
        self.assertFloatsAlmostEqual(dark[0].array, LEVEL)

    def testGetDarkReadPreservesDtype(self):
        dark = ImageCube.readFits(self.path)
        self.assertEqual(self.task.getDarkRead(dark, 0).dtype, np.float32)

    def testGetDarkCubeBacksOutGainAndIsRepeatable(self):
        dark = ImageCube.readFits(self.path)
        self.assertFloatsAlmostEqual(self.task.getDarkCube(dark, NREADS), LEVEL / GAIN)
        self.assertFloatsAlmostEqual(self.task.getDarkCube(dark, NREADS), LEVEL / GAIN)

    def testGetDarkCubeTakesLeadingReads(self):
        """A shorter exposure uses the leading planes of the dark."""
        dark = ImageCube.readFits(self.path)
        self.assertEqual(self.task.getDarkCube(dark, NREADS - 1).shape[0], NREADS - 1)

    def testGetDarkCubeTooShortRaises(self):
        """An exposure with more reads than the dark must not pass silently."""
        dark = ImageCube.readFits(self.path)
        with self.assertRaises(KeyError):
            self.task.getDarkCube(dark, NREADS + 1)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
