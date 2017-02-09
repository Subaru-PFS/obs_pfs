#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python PfsFiberTrace.py
or
   python
   >>> import FiberTrace; FiberTrace.run()
"""
import numpy as np
import os
import unittest
import lsst.daf.persistence as dafPersist
import lsst.utils
import lsst.utils.tests as tests
import pfs.drp.stella as drpStella
import pfs.drp.stella.utils as drpStellaUtils
from pfs.datamodel.pfsFiberTrace import PfsFiberTrace
from pfs.drp.stella.datamodelIO import PfsFiberTraceIO
from pfs.drp.stella.utils import makeFiberTraceSet

class PfsFiberTraceTestCase(tests.TestCase):
    """A test case for comparing a reconstructed FiberTraceSet to the original"""

    def setUp(self):
        drpStellaDataDir = lsst.utils.getPackageDir("drp_stella_data")
        self.butlerDir = os.path.join(drpStellaDataDir,"tests/data/PFS")
        self.butler = dafPersist.Butler(self.butlerDir)
        self.dataId = dict(visit=104, spectrograph=1, arm='r', dateObs='2016-11-11')
        self.flat = self.butler.get("postISRCCD", self.dataId, immediate=True)

        self.ftffc = drpStella.FiberTraceFunctionFindingControl()
        self.ftffc.fiberTraceFunctionControl.order = 5
        self.ftffc.fiberTraceFunctionControl.xLow = -5
        self.ftffc.fiberTraceFunctionControl.xHigh = 5
        self.ftffc.fiberTraceFunctionControl.nRows = self.flat.getHeight()
        self.ftffc.fiberTraceFunctionControl.interpolation = "POLYNOMIAL"
        self.ftffc.apertureFWHM = 2.6
        self.ftffc.signalThreshold = 10
        self.ftffc.nTermsGaussFit = 3
        self.ftffc.saturationLevel = 65550.
        self.ftffc.minLength = 3880
        self.ftffc.maxLength = 3930
        self.ftffc.nLost = 20

        self.ftpfc = drpStella.FiberTraceProfileFittingControl()

        # This particular flatfile has 12 FiberTraces scattered over the whole CCD
        # If in the future the test data change we need to change these numbers
        self.nFiberTraces = 11
        self.fileNameFormat = "pfsFiberTrace-%10s-0-%1d%1s.fits"

    def tearDown(self):
        del self.flat
        del self.ftffc

    def testPfsFiberTrace(self):
        """Test that we can create a pfsFiberTrace"""

        fiberTraceSet = drpStella.findAndTraceAperturesF(self.flat.getMaskedImage(), self.ftffc)

        """Check that we found self.nFiberTraces FiberTraces"""
        self.assertEqual(fiberTraceSet.size(), self.nFiberTraces)

        """Check that we can set the FiberTraceProfileFittingControl ftpfc"""
        self.assertTrue(fiberTraceSet.setFiberTraceProfileFittingControl(self.ftpfc))

        """Check that we can calculate the spatial profile for all the FiberTraces in fiberTraceSet"""
        self.assertTrue(fiberTraceSet.calcProfileAllTraces())

        """Check that we can write a PfsFiberTrace from the fiberTraceSet"""
        dataId = dict(calibDate=self.dataId['dateObs'], spectrograph=self.dataId['spectrograph'], arm=self.dataId['arm'])
        pfsFiberTrace = drpStellaUtils.createPfsFiberTrace(dataId, fiberTraceSet, self.flat.getHeight())
        print pfsFiberTrace.obsDate, pfsFiberTrace.spectrograph, pfsFiberTrace.arm
        fileName = os.path.join('.',PfsFiberTrace.fileNameFormat % (self.dataId['dateObs'], self.dataId['spectrograph'], self.dataId['arm']))
        pfsFiberTraceIO = PfsFiberTraceIO(pfsFiberTrace)
        pfsFiberTraceIO.writeFits(fileName)

        """Check that we can read a PfsFiberTrace back in"""
        pfsFiberTraceNew = PfsFiberTrace(self.dataId['dateObs'], self.dataId['spectrograph'], self.dataId['arm'])
        pfsFiberTraceNew.read()
        print 'pfsFiberTraceNew = ',pfsFiberTraceNew

        self.assertEqual(pfsFiberTrace.obsDate, pfsFiberTraceNew.obsDate)
        self.assertEqual(pfsFiberTrace.spectrograph, pfsFiberTraceNew.spectrograph)
        self.assertEqual(pfsFiberTrace.arm, pfsFiberTraceNew.arm)
        self.assertAlmostEqual(pfsFiberTrace.fwhm, pfsFiberTraceNew.fwhm)
        self.assertAlmostEqual(pfsFiberTrace.threshold, pfsFiberTraceNew.threshold)
        self.assertEqual(pfsFiberTrace.nTerms, pfsFiberTraceNew.nTerms)
        self.assertAlmostEqual(pfsFiberTrace.saturationLevel, pfsFiberTraceNew.saturationLevel)
        self.assertEqual(pfsFiberTrace.minLength, pfsFiberTraceNew.minLength)
        self.assertEqual(pfsFiberTrace.maxLength, pfsFiberTraceNew.maxLength)
        self.assertEqual(pfsFiberTrace.nLost, pfsFiberTraceNew.nLost)
        self.assertEqual(pfsFiberTrace.traceFunction, pfsFiberTraceNew.traceFunction)
        self.assertEqual(pfsFiberTrace.order, pfsFiberTraceNew.order)
        self.assertAlmostEqual(pfsFiberTrace.xLow, pfsFiberTraceNew.xLow)
        self.assertAlmostEqual(pfsFiberTrace.xHigh, pfsFiberTraceNew.xHigh)
        self.assertEqual(pfsFiberTrace.nCutLeft, pfsFiberTraceNew.nCutLeft)
        self.assertEqual(pfsFiberTrace.nCutRight, pfsFiberTraceNew.nCutRight)
        self.assertEqual(pfsFiberTrace.interpol, pfsFiberTraceNew.interpol)
        self.assertEqual(pfsFiberTrace.swathLength, pfsFiberTraceNew.swathLength)
        self.assertEqual(pfsFiberTrace.overSample, pfsFiberTraceNew.overSample)
        self.assertAlmostEqual(pfsFiberTrace.maxIterSF, pfsFiberTraceNew.maxIterSF)
        self.assertAlmostEqual(pfsFiberTrace.maxIterSig, pfsFiberTraceNew.maxIterSig)
        self.assertAlmostEqual(pfsFiberTrace.lambdaSF, pfsFiberTraceNew.lambdaSF)
        self.assertAlmostEqual(pfsFiberTrace.lambdaSP, pfsFiberTraceNew.lambdaSP)
        self.assertAlmostEqual(pfsFiberTrace.lambdaWing, pfsFiberTraceNew.lambdaWing)
        self.assertAlmostEqual(pfsFiberTrace.lSigma, pfsFiberTraceNew.lSigma)
        self.assertAlmostEqual(pfsFiberTrace.uSigma, pfsFiberTraceNew.uSigma)

        """Create new FiberTraceSet"""
        ftsNew = makeFiberTraceSet(pfsFiberTraceNew, self.flat.getMaskedImage())

        self.assertEqual(fiberTraceSet.size(), ftsNew.size())
        for iFt in range(ftsNew.size()):
            ft = fiberTraceSet.getFiberTrace(iFt)
            ftNew = ftsNew.getFiberTrace(iFt)
            ftf = ft.getFiberTraceFunction()
            ftfNew = ftNew.getFiberTraceFunction()
            ftfc = ftf.fiberTraceFunctionControl
            ftfcNew = ftfNew.fiberTraceFunctionControl
            self.assertEqual(ft.getITrace(), ftNew.getITrace())
            self.assertEqual(ft.isTraceSet(), ftNew.isTraceSet())
            self.assertEqual(ft.isProfileSet(), ftNew.isProfileSet())
            self.assertEqual(ft.isFiberTraceProfileFittingControlSet(), ftNew.isFiberTraceProfileFittingControlSet())
            self.assertEqual(ft.getWidth(), ftNew.getWidth())
            self.assertEqual(ft.getHeight(), ftNew.getHeight())
            self.assertEqual(ftf.yCenter, ftfNew.yCenter)
            self.assertAlmostEqual(ftf.xCenter, ftfNew.xCenter)
            self.assertEqual(ftf.yLow, ftfNew.yLow)
            self.assertEqual(ftf.yHigh, ftfNew.yHigh)

            self.assertEqual(ftfc.interpolation, ftfcNew.interpolation)
            self.assertEqual(ftfc.order, ftfcNew.order)
            self.assertAlmostEqual(ftfc.xLow, ftfcNew.xLow)
            self.assertAlmostEqual(ftfc.xHigh, ftfcNew.xHigh)
            self.assertEqual(ftfc.nPixCutLeft, ftfcNew.nPixCutLeft)
            self.assertEqual(ftfc.nPixCutRight, ftfcNew.nPixCutRight)
            self.assertEqual(ftfc.nRows, ftfcNew.nRows)

            coeffs = ft.getTraceCoefficients()
            coeffsNew = ftNew.getTraceCoefficients()
            self.assertEqual(coeffs.shape[0], coeffsNew.shape[0])
            for iCoeff in range(coeffs.shape[0]):
                self.assertAlmostEqual(coeffs[iCoeff], coeffsNew[iCoeff], places=3)

            xCenters = ft.getXCenters()
            xCentersNew = ftNew.getXCenters()
            self.assertEqual(xCenters.shape[0], xCentersNew.shape[0])
            for iY in range(xCenters.shape[0]):
                self.assertAlmostEqual(xCenters[iY], xCentersNew[iY], places=3)

            profile = ft.getProfile().getArray()
            profileNew = ftNew.getProfile().getArray()
            self.assertEqual(profile.shape[0], profileNew.shape[0])
            self.assertEqual(profile.shape[1], profileNew.shape[1])
            for iX in range(profile.shape[0]):
                for iY in range(profile.shape[1]):
                    if np.abs(profile[iX,iY] - profileNew[iX,iY]) > 0.1:
                        print 'FiberTrace ',iFt,': large difference detected: profile[',iX,',',iY,'] = ',profile[iX,iY],', profileNew[',iX,',',iY,'] = ',profileNew[iX,iY],'      xCenters[',iX,'] = ',xCenters[iX],', xCentersNew[',iX,'] = ',xCentersNew[iX]
                    self.assertAlmostEqual(profile[iX,iY], profileNew[iX,iY], places=3)

            image = ft.getImage().getArray()
            imageNew = ftNew.getImage().getArray()
            self.assertEqual(image.shape[0], imageNew.shape[0])
            self.assertEqual(image.shape[1], imageNew.shape[1])
            for iX in range(image.shape[0]):
                for iY in range(image.shape[1]):
                    if np.abs(image[iX,iY] - imageNew[iX,iY]) > 0.1:
                        print 'FiberTrace ',iFt,': large difference detected: image[',iX,',',iY,'] = ',image[iX,iY],', imageNew[',iX,',',iY,'] = ',imageNew[iX,iY],'      xCenters[',iX,'] = ',xCenters[iX],', xCentersNew[',iX,'] = ',xCentersNew[iX]
                    self.assertAlmostEqual(image[iX,iY], imageNew[iX,iY], places=3)

            self.assertEqual(ft.getITrace(), ftNew.getITrace())
            self.assertAlmostEqual(ft.getFiberTraceFunction().xCenter, ftNew.getFiberTraceFunction().xCenter)
            self.assertEqual(ft.getFiberTraceFunction().yCenter, ftNew.getFiberTraceFunction().yCenter)
            self.assertEqual(ft.getFiberTraceFunction().yLow, ftNew.getFiberTraceFunction().yLow)
            self.assertEqual(ft.getFiberTraceFunction().yHigh, ftNew.getFiberTraceFunction().yHigh)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(PfsFiberTraceTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--display", '-d', default=False, action="store_true", help="Activate display?")
    parser.add_argument("--verbose", '-v', type=int, default=0, help="Verbosity level")
    args = parser.parse_args()
    display = args.display
    verbose = args.verbose
    run(True)
