#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python FiberTrace.py
or
   python
   >>> import FiberTrace; FiberTrace.run()
"""
import os
import unittest
import lsst.daf.persistence as dafPersist
import lsst.utils
import lsst.utils.tests as tests

class DitheredFlatsTestCase(tests.TestCase):
    """A test case for measuring FiberTrace quantities"""

    def setUp(self):
        self.xOffsetHdrKeyWordSims = 'sim.slit.xoffset'
        self.xOffsetHdrKeyWordLAM = 'W_FCA_SHIFT'
        drpStellaDataDir = lsst.utils.getPackageDir("drp_stella_data")
        self.butler = dafPersist.Butler(os.path.join(drpStellaDataDir,"tests/data/PFS/"))
        
    def tearDown(self):
        del self.butler

    def testSimulatedDitheredFlat(self):
        """Test that we can read the xOffset"""
        dataId = dict(field="FLAT", visit=29, spectrograph=2, arm="r")
        simFlat = self.butler.get("postISRCCD", dataId, immediate=True)
        xOffset = simFlat.getMetadata().get(self.xOffsetHdrKeyWordSims)
        self.assertAlmostEqual(xOffset, 0.)

    def testLAMExposure(self):
        """Test that we can read the xOffset"""
        dataId = dict(field="OBJECT", visit=300, spectrograph=1, arm="r")
        lamExp = self.butler.get("raw", dataId, immediate=True)
        xOffset = lamExp.getMetadata().get(self.xOffsetHdrKeyWordLAM)
        self.assertAlmostEqual(xOffset, 0.)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(DitheredFlatsTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--verbose", '-v', type=int, default=0, help="Verbosity level")
    args = parser.parse_args()
    verbose = args.verbose
    run(True)
