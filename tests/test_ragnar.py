import os
import tempfile
import unittest

import numpy as np
import pandas as pd

import lsst.utils.tests

from testUtils import HAS_DRP_STELLA, requireDrpStella

if HAS_DRP_STELLA:
    from lsst.obs.pfs.ragnar import SIGMA_WAV_THRESHOLD, readSweepFile
    from pfs.drp.stella.referenceLine import ReferenceLineSource, ReferenceLineStatus


@requireDrpStella
class RagnarTestCase(lsst.utils.tests.TestCase):
    """Tests of converting a RAGNAR sweep file to a ReferenceLineSet"""

    def setUp(self):
        self.tempDir = tempfile.mkdtemp()

    def tearDown(self):
        for filename in os.listdir(self.tempDir):
            os.remove(os.path.join(self.tempDir, filename))
        os.rmdir(self.tempDir)

    def writeSweep(self, withArm: bool) -> str:
        """Write a small synthetic sweep file

        Parameters
        ----------
        withArm : `bool`
            Include an ``arm`` column, as in a unified sweep file?

        Returns
        -------
        filename : `str`
            Path to the sweep file.
        """
        data = pd.DataFrame(
            dict(
                wav=[557.734, 589.158, 761.5, 834.5, 834.5],
                intensity=[100.0, 20.0, 5.0, 50.0, 49.0],
                sigma_intensity=[1.0, 2.0, 3.0, 4.0, 5.0],
                label=["OI(5577)", "NaI(D1)", "O2(A)|(9-5)P1e(2.5)", "(6-2)P1e(3.5)", "(6-2)P1e(3.5)"],
                clean=[True, True, True, False, True],
                robust=[True, True, True, True, True],
                sigma_wav=[0.0005, 0.00138, 0.0005, 0.0005, 2 * SIGMA_WAV_THRESHOLD],
            )
        )
        if withArm:
            data["arm"] = ["b", "b", "r", "r", "m"]
        filename = os.path.join(self.tempDir, f"sweep_{withArm}.pkl")
        data.to_pickle(filename)
        return filename

    def testUnified(self):
        """Test reading a unified sweep file, selecting an arm"""
        filename = self.writeSweep(True)
        lines = readSweepFile(filename, "r")
        np.testing.assert_array_equal(lines.wavelength, [761.5, 834.5])
        self.assertEqual(list(lines.description), ["O2|OH", "OH"])
        self.assertEqual(list(lines.status), [ReferenceLineStatus.GOOD, ReferenceLineStatus.BLEND])  # unclean
        np.testing.assert_array_equal(lines.intensityErr, [3.0, 4.0])
        self.assertTrue(np.all(lines.source == ReferenceLineSource.RAGNAR))

        lines = readSweepFile(filename, "m")
        self.assertEqual(list(lines.status), [ReferenceLineStatus.BLEND])  # sigma_wav too large

        with self.assertRaises(RuntimeError):
            readSweepFile(filename)  # Arm required
        with self.assertRaises(RuntimeError):
            readSweepFile(filename, "n")  # Arm not present

    def testPerArm(self):
        """Test reading an old per-arm sweep file"""
        filename = self.writeSweep(False)
        lines = readSweepFile(filename)
        self.assertEqual(len(lines), 5)
        self.assertEqual(list(lines.description[:2]), ["OI", "NaI"])
        # The NaI D lines pass the relaxed sigma_wav threshold
        self.assertEqual(lines.status[1], ReferenceLineStatus.GOOD)
        with self.assertRaises(RuntimeError):
            readSweepFile(filename, "b")  # No arm column


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
