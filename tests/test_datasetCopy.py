"""Tests for copying/linking datasets between collections (``datasetCopy``).

A dataset's RUN is immutable, so promoting a calibration into a shared collection
means reading it and writing it out again. These tests cover that copy, the
no-copy alternative, inheriting the source certification's validity period, and
the guards around resuming an interrupted copy.
"""

import os
import shutil
import sys
import tempfile
import unittest

import astropy.time
import numpy as np

import lsst.utils.tests
from lsst.daf.butler import CollectionType, Timespan
from lsst.daf.butler._exceptions import MissingDatasetTypeError
from lsst.daf.butler.registry import ConflictingDefinitionError

from lsst.obs.pfs import datasetCopy
from lsst.obs.pfs.imageCube import ImageCube

from testUtils import closeButler, loadScript

BEGIN = astropy.time.Time("2026-07-08T12:00:00", format="isot", scale="tai")


class DatasetCopyTestCase(lsst.utils.tests.TestCase):
    """Round-trip through a throwaway Gen3 butler."""

    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="datasetCopy-repo-")
        self.butler = loadScript("makePfsTestRepo").makePfsTestRepo(
            os.path.join(self.root, "repo"))
        self.srcRun = "u/test/src/put"
        self.butler.registry.registerRun(self.srcRun)
        self.dataIds = [dict(instrument="PFS", arm="n", spectrograph=s) for s in (1, 2)]
        self.data = {}
        self.srcRefs = []
        for i, dataId in enumerate(self.dataIds):
            data = np.full((2, 3, 3), i + 1, dtype="f4")
            self.data[dataId["spectrograph"]] = data
            cube = ImageCube.fromCube(data, dict(W_H4IRPN=1))
            self.srcRefs.append(self.butler.put(cube, "nirDark", dataId, run=self.srcRun))

    def tearDown(self):
        closeButler(self)
        shutil.rmtree(self.root, ignore_errors=True)

    def certifySource(self, calib="u/test/src/nirDarkGen.20260708a"):
        self.butler.registry.registerCollection(calib, CollectionType.CALIBRATION)
        self.butler.registry.certify(calib, self.srcRefs, Timespan(BEGIN, None))
        return calib

    def testCopyMakesIndependentDatasets(self):
        outputRun = "PFS/calib/test/nirDarkGen.20260708a/put"
        newRefs = datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun],
                                       outputRun=outputRun)
        self.assertEqual(len(newRefs), 2)
        for ref in newRefs:
            self.assertEqual(ref.run, outputRun)
        # New dataset IDs: these are copies, not the originals re-pointed.
        self.assertFalse({r.id for r in newRefs} & {r.id for r in self.srcRefs})
        # And the pixels survived the round trip.
        for dataId in self.dataIds:
            cube = self.butler.get("nirDark", dataId, collections=outputRun)
            self.assertFloatsAlmostEqual(cube.getImageCube(),
                                         self.data[dataId["spectrograph"]])

    def testCopySurvivesSourceRemoval(self):
        """The point of copying: the destination must not depend on the source."""
        outputRun = "PFS/calib/test/independent/put"
        datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun], outputRun=outputRun)
        self.butler.removeRuns([self.srcRun])
        for dataId in self.dataIds:
            cube = self.butler.get("nirDark", dataId, collections=outputRun)
            self.assertFloatsAlmostEqual(cube.getImageCube(),
                                         self.data[dataId["spectrograph"]])

    def testTimespanInheritedFromSourceCertification(self):
        calib = self.certifySource()
        outputRun = "PFS/calib/test/inherit/put"
        destCalib = "PFS/calib/test/inherit"
        datasetCopy.transfer(self.butler, ["nirDark"], [calib], outputRun=outputRun,
                             calibCollection=destCalib)
        registry = self.butler.registry
        self.assertEqual(registry.getCollectionType(destCalib), CollectionType.CALIBRATION)
        # Valid after the inherited begin date, and open-ended.
        far = astropy.time.Time("2099-01-01T00:00:00", format="isot", scale="tai")
        for dataId in self.dataIds:
            self.assertIsNotNone(registry.findDataset("nirDark", dataId, collections=destCalib,
                                                      timespan=Timespan(far, far)))
        # Not valid before it was taken: the inherited begin date was respected.
        before = astropy.time.Time("2020-01-01T00:00:00", format="isot", scale="tai")
        for dataId in self.dataIds:
            self.assertIsNone(registry.findDataset("nirDark", dataId, collections=destCalib,
                                                   timespan=Timespan(before, before)))

    def testExplicitBeginDateOverridesSource(self):
        calib = self.certifySource()
        destCalib = "PFS/calib/test/explicit"
        datasetCopy.transfer(self.butler, ["nirDark"], [calib],
                             outputRun="PFS/calib/test/explicit/put",
                             calibCollection=destCalib, beginDate="2030-01-01T00:00:00")
        registry = self.butler.registry
        # 2027 falls after the source begin date but before the explicit one.
        mid = astropy.time.Time("2027-01-01T00:00:00", format="isot", scale="tai")
        self.assertIsNone(registry.findDataset("nirDark", self.dataIds[0],
                                               collections=destCalib,
                                               timespan=Timespan(mid, mid)))

    def testUncertifiedSourceCannotInheritTimespan(self):
        """Inheriting a validity period from a plain RUN is not possible."""
        with self.assertRaises(RuntimeError) as cm:
            datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun],
                                 outputRun="PFS/calib/test/nospan/put",
                                 calibCollection="PFS/calib/test/nospan")
        self.assertIn("not certified", str(cm.exception))

    def testNoCopyLeavesDatasetsInSourceRun(self):
        calib = self.certifySource()
        destCalib = "PFS/calib/test/linked"
        newRefs = datasetCopy.transfer(self.butler, ["nirDark"], [calib],
                                       outputRun=None, calibCollection=destCalib)
        for ref in newRefs:
            self.assertEqual(ref.run, self.srcRun)
        self.assertEqual({r.id for r in newRefs}, {r.id for r in self.srcRefs})

    def testRepeatedCopyConflictsUnlessSkipping(self):
        outputRun = "PFS/calib/test/twice/put"
        datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun], outputRun=outputRun)
        with self.assertRaises(ConflictingDefinitionError):
            datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun], outputRun=outputRun)
        # Resuming an interrupted copy is allowed, and maps to the existing datasets.
        newRefs = datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun],
                                       outputRun=outputRun, skipExisting=True)
        self.assertEqual(len(newRefs), 2)
        for ref in newRefs:
            self.assertEqual(ref.run, outputRun)

    def testDryRunWritesNothing(self):
        outputRun = "PFS/calib/test/dry/put"
        self.assertEqual(
            datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun],
                                 outputRun=outputRun, dryRun=True),
            [])
        with self.assertRaises(Exception):
            self.butler.registry.getCollectionType(outputRun)

    def testMissingDatasetsRaise(self):
        with self.assertRaises(RuntimeError) as cm:
            datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun],
                                 outputRun="PFS/calib/test/none/put",
                                 where="instrument='PFS' AND spectrograph = 3")
        self.assertIn("no datasets", str(cm.exception))

    def testWhereRestrictsCopy(self):
        outputRun = "PFS/calib/test/where/put"
        newRefs = datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun],
                                       outputRun=outputRun,
                                       where="instrument='PFS' AND spectrograph = 2")
        self.assertEqual(len(newRefs), 1)
        self.assertEqual(newRefs[0].dataId["spectrograph"], 2)

    def testCopiedPixelsSurviveStorageClassRoundTrip(self):
        """Regression: copying must not go through ImageCube.get/put.

        A fetched ImageCube is lazily backed by its file and holds no images, so
        writeFits would persist a header with no image HDUs -- losing every read
        without raising. Copying the artifact sidesteps the storage class.
        """
        outputRun = "PFS/calib/test/pixels/put"
        datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun], outputRun=outputRun)
        for dataId in self.dataIds:
            cube = self.butler.get("nirDark", dataId, collections=outputRun)
            self.assertEqual(cube.nreads, 2)
            self.assertFloatsAlmostEqual(cube.getImageCube(),
                                         self.data[dataId["spectrograph"]])
            self.assertEqual(cube.metadata.get("W_H4IRPN"), 1)

    def testDisassembledCompositeRejected(self):
        """A dataset with no single artifact cannot be copied silently."""
        ref = self.srcRefs[0]

        class FakeUris:
            primaryURI = None
            componentURIs = {"image": "file:///nowhere.fits"}

        original = self.butler.getURIs
        self.butler.getURIs = lambda _ref: FakeUris()
        try:
            with self.assertRaises(NotImplementedError):
                datasetCopy.artifactUri(self.butler, ref)
        finally:
            self.butler.getURIs = original

    def testScratchCubesRefusedByDefault(self):
        """The bulky rawISRCube intermediates must not be promoted by accident."""
        with self.assertRaises(ValueError) as cm:
            datasetCopy.transfer(self.butler, ["rawISRCube"], [self.srcRun],
                                 outputRun="PFS/calib/test/scratch/put")
        self.assertIn("intermediate product", str(cm.exception))
        # Nothing was even queried, let alone written.
        with self.assertRaises(Exception):
            self.butler.registry.getCollectionType("PFS/calib/test/scratch/put")

    def testScratchCubesCopyableWhenAsked(self):
        """The refusal is a guard, not a prohibition: allowIntermediates lifts it.

        The throwaway repo has no ``rawISRCube`` dataset type (it needs visit
        dimensions this repo lacks), so the call still fails -- but on the missing
        type, having got past the intermediate guard.
        """
        with self.assertRaises(MissingDatasetTypeError):
            datasetCopy.transfer(self.butler, ["rawISRCube"], [self.srcRun],
                                 outputRun="u/test/scratchCopy/put",
                                 allowIntermediates=True)

    def testIntermediateRefusedEvenAlongsideWantedType(self):
        """Naming a good type does not smuggle an intermediate along with it."""
        with self.assertRaises(ValueError) as cm:
            datasetCopy.transfer(self.butler, ["nirDark", "rawISRCube"], [self.srcRun],
                                 outputRun="PFS/calib/test/mixed/put")
        self.assertIn("rawISRCube", str(cm.exception))

    def testChainedInputCopiesOnlyNamedDatasetType(self):
        """A chain holding other products must not leak them into the copy.

        This is what keeps a scratchCubes collection out of a promotion: only the
        dataset types named on the command line are ever queried.
        """
        otherRun = "u/test/other"
        self.butler.registry.registerRun(otherRun)
        self.butler.put(ImageCube.fromCube(np.zeros((2, 3, 3), dtype="f4"), dict(W_H4IRPN=4)),
                        "nirDark_irp4", self.dataIds[0], run=otherRun)
        chain = "u/test/everything"
        self.butler.registry.registerCollection(chain, CollectionType.CHAINED)
        self.butler.registry.setCollectionChain(chain, [self.srcRun, otherRun])

        outputRun = "PFS/calib/test/fromchain/put"
        newRefs = datasetCopy.transfer(self.butler, ["nirDark"], [chain], outputRun=outputRun)
        self.assertEqual({ref.datasetType.name for ref in newRefs}, {"nirDark"})
        self.assertEqual(
            self.butler.query_datasets("nirDark_irp4", collections=[outputRun],
                                       explain=False, limit=1),
            [])

    def testOutputRunMustBeRun(self):
        chain = "u/test/chain"
        self.butler.registry.registerCollection(chain, CollectionType.CHAINED)
        with self.assertRaises(ValueError) as cm:
            datasetCopy.transfer(self.butler, ["nirDark"], [self.srcRun], outputRun=chain)
        self.assertIn("CHAINED", str(cm.exception))


class CopyDatasetsScriptTestCase(lsst.utils.tests.TestCase):
    """Argument guards in the ``copyDatasets`` wrapper."""

    def setUp(self):
        self.script = loadScript("copyDatasets")
        self.origArgv = sys.argv

    def tearDown(self):
        sys.argv = self.origArgv

    def run_main(self, argv):
        sys.argv = ["copyDatasets"] + argv
        self.script.main()

    def testOutputRunRequired(self):
        with self.assertRaises(SystemExit):
            self.run_main(["/repo", "--input", "in", "--dataset-type", "nirDark"])

    def testNoCopyRejectsOutputRun(self):
        with self.assertRaises(SystemExit):
            self.run_main(["/repo", "--input", "in", "--dataset-type", "nirDark",
                           "--no-copy", "--certify", "calib", "--output-run", "run"])

    def testNoCopyRequiresCertify(self):
        with self.assertRaises(SystemExit):
            self.run_main(["/repo", "--input", "in", "--dataset-type", "nirDark",
                           "--no-copy"])

    def testDatasetTypeHasNoDefault(self):
        """Nothing is copied unless a dataset type is named: no wildcard, no default."""
        with self.assertRaises(SystemExit):
            self.run_main(["/repo", "--input", "in", "--output-run", "out"])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
