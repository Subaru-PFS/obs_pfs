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
import shutil
import unittest

import numpy as np

import lsst.daf.persistence as dafPersist
import lsst.log as log
from lsst.obs.pfs.ingest import PfsParseTask
from lsst.obs.pfs.pfsMapper import PfsMapper
import lsst.utils
import lsst.utils.tests as tests

class MapperTestCase(tests.TestCase):
    """A test case for the mapper and butler"""

    def setUp(self):
        self.logger = log.Log.getLogger('MapperTests')
        try:
            self.drpStellaDataDir = lsst.utils.getPackageDir("drp_stella_data")
        except Exception, e:
            self.drpStellaDataDir = None
            self.logger.warn('cannot run tests as drp_stella_data is not setup')
        if self.drpStellaDataDir is not None:
            self.butlerDir = os.path.join(self.drpStellaDataDir,"tests/temp")
            if os.path.exists(self.butlerDir):
                shutil.rmtree(self.butlerDir)
            os.makedirs(self.butlerDir)
            os.makedirs(os.path.join(self.butlerDir,'CALIB'))

            # Write mapper
            mapper = "lsst.obs.pfs.PfsMapper"
            text_file = open(os.path.join(self.butlerDir,"_mapper"), "w")
            text_file.write(mapper)
            text_file.close()

            self.rawFiles = os.path.join(self.drpStellaDataDir,"tests/data/raw/*.fits")
            self.regFile = os.path.join(self.butlerDir,'registry.sqlite3')

            self.site = 'S'
            self.category = 'A'
            self.visit = 4
            self.arms = [None, 'b', 'r', 'n', 'm']
            self.minArmNum = 1
            self.maxArmNum = 4
            self.minSpectrographNum = 1
            self.maxSpectrographNum = 4

    def tearDown(self):
        """clean up"""
        if self.drpStellaDataDir is not None and os.path.exists(self.butlerDir):
            shutil.rmtree(self.butlerDir)

    def testFileNameInterpretation(self):
        if self.drpStellaDataDir is not None:
            pfsParseTask = PfsParseTask(name='name')
            mapper = PfsMapper(
                root=os.path.join(self.drpStellaDataDir,'tests/data/PFS')
            )
            for armNum in np.arange(self.minArmNum, self.maxArmNum + 1, 1):
                for spectrograph in np.arange(self.minSpectrographNum, self.maxSpectrographNum + 1, 1):
                    arm = self.arms[armNum]
                    dataId = {'spectrograph': spectrograph, 'arm': arm}
                    ccd = mapper.getDetectorId(dataId)
                    fileNameFormat = 'PF%1s%1s%06d%1d%1d.fits'
                    fileName = fileNameFormat % (self.site, self.category, self.visit, spectrograph, armNum)
                    dirName = os.path.dirname( self.rawFiles )
                    fileNameWithDir = os.path.join( dirName, fileName )
                    info = pfsParseTask.getInfo( fileNameWithDir )
                    self.assertEqual(info[0]['site'], self.site)
                    self.assertEqual(info[0]['category'], self.category)
                    self.assertEqual(info[0]['visit'], self.visit)
                    self.assertEqual(info[0]['spectrograph'], spectrograph)
                    self.assertEqual(info[0]['arm'], arm)
                    self.assertEqual(info[0]['filter'], arm)
                    self.assertEqual(info[0]['ccd'], ccd)

    def testDetectorId(self):
        if self.drpStellaDataDir is not None:
            butler = dafPersist.Butler(self.butlerDir)
            mapper = butler.getMapperClass(self.butlerDir)(root=self.butlerDir)

            #detector ids for the individual spectrographs and arms
            blues=[0,3,6,9]
            reds=[1,4,7,10]
            nirs=[2,5,8,11]
            ones=[0,1,2]
            twos=[3,4,5]
            threes=[6,7,8]
            fours=[9,10,11]

            for armNum in np.arange(self.minArmNum, self.maxArmNum + 1, 1):
                for spectrograph in np.arange(self.minSpectrographNum, self.maxSpectrographNum + 1, 1):
                    arm = self.arms[armNum]
                    dataId = {'visit': self.visit, 'spectrograph': spectrograph, 'arm':arm}
                    detectorId = mapper.getDetectorId(dataId)
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
