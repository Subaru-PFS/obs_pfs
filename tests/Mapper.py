#!/usr/bin/env python
"""
Tests for the mapper/butler

Run with:
   python Mapper.py
or
   python
   >>> import Mapper; Mapper.run()
"""

import os
import unittest
import shutil
import numpy as np
import lsst.utils
import lsst.utils.tests as tests
from lsst.obs.pfs.ingest import PfsParseTask

class MapperTestCase(tests.TestCase):
    """A test case for the mapper and butler"""

    def setUp(self):
        drpStellaDataDir = lsst.utils.getPackageDir("drp_stella_data")
        self.butlerDir = os.path.join(drpStellaDataDir,"tests/temp")
        if os.path.exists(self.butlerDir):
            shutil.rmtree(self.butlerDir)
        os.makedirs(self.butlerDir)
        os.makedirs(os.path.join(self.butlerDir,'CALIB'))

        """Write mapper"""
        mapper = "lsst.obs.pfs.PfsMapper"
        text_file = open(os.path.join(self.butlerDir,"_mapper"), "w")
        text_file.write(mapper)
        text_file.close()

        self.rawFiles = os.path.join(drpStellaDataDir,"tests/data/raw/*.fits")
        self.regFile = os.path.join(self.butlerDir,'registry.sqlite3')
        self.dataId = {'visit': '4', 'ccd': 1}

    def tearDown(self):
        """clean up"""
#        if os.path.exists(self.butlerDir):
#            shutil.rmtree(self.butlerDir)

    def testFileNameInterpretation(self):
        pfsParseTask = PfsParseTask(name='name')
        site = 'S'
        category = 'A'
        visit = 4
        arms = ['b', 'r', 'n', 'm']
        minArmNum = 1
        maxArmNum = 4
        minSpectrographNum = 1
        maxSpectrographNum = 4
        for armNum in np.arange(minArmNum, maxArmNum + 1, 1):
            for spectrograph in np.arange(minSpectrographNum, maxSpectrographNum + 1, 1):
                ccd = spectrograph - 1
                if int(armNum) == minArmNum + 1 or int(armNum) == minArmNum + 3:
                    ccd += 4
                elif int(armNum) == minArmNum + 2:
                    ccd += 8
                fileNameFormat = 'PF%1s%1s%06d%1d%1d.fits'
                fileName = fileNameFormat % (site, category, visit, spectrograph, armNum)
                dirName = os.path.dirname( self.rawFiles )
                fileNameWithDir = os.path.join( dirName, fileName )
                info = pfsParseTask.getInfo( fileNameWithDir )
                self.assertTrue(info[0]['site'] == site)
                self.assertTrue(info[0]['category'] == category)
                self.assertTrue(info[0]['visit'] == visit)
                self.assertTrue(info[0]['spectrograph'] == spectrograph)
                self.assertTrue(info[0]['arm'] == arms[armNum-1])
                self.assertTrue(info[0]['filter'] == arms[armNum-1])
                self.assertTrue(info[0]['ccd'] == ccd)

    def testDetectorId(self):
        butler = dafPersist.Butler(self.butlerDir)
        mapper = butler.getMapperClass(self.butlerDir)(root=self.butlerDir)
        camera = mapper.camera
        blues=[0,3,6,9]
        reds=[1,4,7,10]
        nirs=[2,5,8,11]
        ones=[0,1,2]
        twos=[3,4,5]
        threes=[6,7,8]
        fours=[9,10,11]
        for armNum in np.arange(self.minArmNum, self.maxArmNum + 1, 1):
            for spectrograph in np.arange(self.minSpectrographNum, self.maxSpectrographNum + 1, 1):
                dataId = {'visit': self.visit, 'spectrograph': spectrograph, 'arm':self.arms[armNum-1]}
                detectorId = mapper._extractDetectorId(dataId)
                if self.arms[armNum-1] in ['b']:
                    self.assertTrue(detectorId in blues)
                elif self.arms[armNum-1] in ['r', 'm']:
                    self.assertTrue(detectorId in reds)
                elif self.arms[armNum-1] in ['n']:
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

                name = mapper._extractDetectorName(dataId)
                for ccd in camera:
                    ccdNum = ccd.getId()
                    ccdName = ccd.getName()
                    if ccdName in [name]:
                        self.assertEqual(ccdNum, detectorId)

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MapperTestCase)
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
