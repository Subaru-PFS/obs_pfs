import os
import sys
import unittest

import numpy as np

from lsst.obs.pfs.pfsMapper import PfsMapper
import lsst.utils.tests


class MapperTestCase(lsst.utils.tests.TestCase):
    """A test case for the mapper and butler"""

    def setUp(self):
        self.mapper = PfsMapper(root=os.path.dirname(__file__))
        self.expId = 12345
        self.arms = [None, 'b', 'r', 'n', 'm']
        self.minArmNum = 1
        self.maxArmNum = 4
        self.minSpectrographNum = 1
        self.maxSpectrographNum = 4

    def testDetectorId(self):
        # detector ids for the individual spectrographs and arms
        blues = [0, 3, 6, 9]
        reds = [1, 4, 7, 10]
        nirs = [2, 5, 8, 11]
        ones = [0, 1, 2]
        twos = [3, 4, 5]
        threes = [6, 7, 8]
        fours = [9, 10, 11]

        for armNum in np.arange(self.minArmNum, self.maxArmNum + 1, 1):
            for spectrograph in np.arange(self.minSpectrographNum, self.maxSpectrographNum + 1, 1):
                arm = self.arms[armNum]
                dataId = {'expId': self.expId, 'spectrograph': spectrograph, 'arm': arm}
                detectorId = self.mapper.getDetectorId(dataId)
                if arm in ['b']:
                    self.assertTrue(detectorId in blues)
                elif arm in ['r', 'm']:
                    self.assertTrue(detectorId in reds)
                elif arm in ['n']:
                    self.assertTrue(detectorId in nirs)
                else:
                    self.assertTrue(False)

                if spectrograph in [1]:
                    self.assertTrue(detectorId in ones)
                elif spectrograph in [2]:
                    self.assertTrue(detectorId in twos)
                elif spectrograph in [3]:
                    self.assertTrue(detectorId in threes)
                elif spectrograph in [4]:
                    self.assertTrue(detectorId in fours)
                else:
                    self.assertTrue(False)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    from argparse import ArgumentParser
    parser = ArgumentParser(__file__)
    parser.add_argument("--display", help="Display backend")
    parser.add_argument("--log", action="store_true", help="Activate verbose logging")
    args, argv = parser.parse_known_args()
    display = args.display
    if args.log:
        import lsst.log  # noqa
        lsst.log.setLevel("pfs.drp.stella.FiberTrace", lsst.log.TRACE)
    unittest.main(failfast=True, argv=[__file__] + argv)
