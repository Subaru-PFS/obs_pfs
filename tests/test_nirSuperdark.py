"""Tests for the combined-NIR-dark core (``lsst.obs.pfs.nirSuperdark``).

Exercises the read-by-read median combine (``makeMasterDark``), its input-
consistency guards, and the butler round-trip of ``saveNirDark`` including the
IRP-ratio-driven dataset-type selection (``nirDark`` vs ``nirDark_irp4``) and the
IRP/timing metadata cards.
"""

import datetime
import os
import shutil
import sys
import tempfile
import unittest
import warnings

import astropy.time
import astropy.units as u
import erfa
import numpy as np

import lsst.resources
import lsst.utils.tests
from lsst.daf.butler import CollectionType, Timespan

from lsst.obs.pfs.imageCube import ImageCube

from testUtils import HAS_DRP_STELLA, closeButler, loadScript, requireDrpStella

if HAS_DRP_STELLA:
    # nirSuperdark imports pfs.drp.stella.calibs.setCalibHeader.
    from lsst.obs.pfs import nirSuperdark


class FakeButler:
    """Minimal butler stand-in for ``makeMasterDark``.

    ``makeMasterDark`` only calls ``getURI('rawISRCube', dataId, visit=...)`` and
    opens the result with ``fitsio``; a ``ResourcePath`` to an on-disk cube (what
    the real butler returns) is all that is needed.
    """

    def __init__(self, paths):
        self.paths = paths  # visit -> file path

    def getURI(self, datasetType, dataId, visit):
        assert datasetType == "rawISRCube"
        return lsst.resources.ResourcePath(self.paths[visit])


def writeCube(path, data, irpN=1, irpOffset=1, frameTime=1.5, nred=None):
    """Write a synthetic rawISRCube FITS file with H4 header cards.

    A rawISRCube has one plane fewer than the ramp has reads, so ``W_H4NRED`` is
    one more than ``len(data)`` unless overridden to test the guard.
    """
    if nred is None:
        nred = len(data) + 1
    metadata = dict(W_H4IRPN=irpN, W_H4IRPO=irpOffset,
                    W_H4NRED=nred, W_H4FRMT=frameTime)
    ImageCube.fromCube(data, metadata).writeFits(path)


@requireDrpStella
class MakeMasterDarkTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="nirSuperdark-")

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def makeInputs(self, offsetsPerDark, nreads=3, shape=(4, 5), irpN=1):
        """Write darks that share a common ramp plus a constant per-dark offset.

        Dark ``v`` at every read is ``base[r] + offsetsPerDark[v]``, where
        ``base[r]`` is a flat image of value ``10*(r+1)`` (the common flux ramp).
        """
        butler = FakeButler({})
        for v, off in enumerate(offsetsPerDark):
            data = np.empty((nreads,) + shape, dtype="f4")
            for r in range(nreads):
                data[r] = 10 * (r + 1) + off
            path = os.path.join(self.root, f"cube_{v}.fits")
            writeCube(path, data, irpN=irpN)
            butler.paths[v] = path
        return butler, list(range(len(offsetsPerDark)))

    def testCombinePreservesCommonRamp(self):
        """The common flux ramp survives; only inter-dark offsets are removed."""
        offsets = [0.0, 10.0, 20.0]  # median 10 -> preserved in the result
        butler, visits = self.makeInputs(offsets, nreads=3, shape=(4, 5))
        newCube, perReadOffsets, rampInfo = nirSuperdark.makeMasterDark(
            butler, dict(instrument="PFS", arm="n", spectrograph=1), visits)

        self.assertEqual(newCube.shape, (3, 4, 5))
        for r in range(3):
            # result == base + median(offsets) == 10*(r+1) + 10
            self.assertFloatsAlmostEqual(newCube[r], 10 * (r + 1) + 10.0)
            # offsets recorded relative to the per-read median
            np.testing.assert_allclose(perReadOffsets[r], [-10.0, 0.0, 10.0])

        self.assertEqual(rampInfo["W_H4IRPN"], 1)
        # One more read than there are planes: the first read is the reference.
        self.assertEqual(rampInfo["W_H4NRED"], 4)
        self.assertEqual(rampInfo["W_H4FRMT"], 1.5)

    def testIdenticalDarksCombineToInput(self):
        butler, visits = self.makeInputs([0.0, 0.0], nreads=4, shape=(3, 3))
        newCube, offsets, _ = nirSuperdark.makeMasterDark(
            butler, dict(instrument="PFS", arm="n", spectrograph=1), visits)
        for r in range(4):
            self.assertFloatsAlmostEqual(newCube[r], 10 * (r + 1))
            np.testing.assert_allclose(offsets[r], [0.0, 0.0])

    def testMismatchedReadCountRaises(self):
        butler = FakeButler({})
        for v, nreads in ((0, 3), (1, 2)):
            data = np.ones((nreads, 3, 3), dtype="f4")
            path = os.path.join(self.root, f"cube_{v}.fits")
            writeCube(path, data)
            butler.paths[v] = path
        with self.assertRaises(ValueError):
            nirSuperdark.makeMasterDark(
                butler, dict(instrument="PFS", arm="n", spectrograph=1), [0, 1])

    def testPlaneCountMustMatchReadCount(self):
        """A rawISRCube must hold exactly W_H4NRED-1 planes."""
        butler = FakeButler({})
        data = np.ones((3, 3, 3), dtype="f4")
        path = os.path.join(self.root, "cube_0.fits")
        writeCube(path, data, nred=3)  # 3 planes claiming 3 reads; should be 4
        butler.paths[0] = path
        with self.assertRaises(ValueError) as cm:
            nirSuperdark.makeMasterDark(
                butler, dict(instrument="PFS", arm="n", spectrograph=1), [0])
        self.assertIn("W_H4NRED", str(cm.exception))

    def testMismatchedIrpRaises(self):
        butler = FakeButler({})
        for v, irpN in ((0, 1), (1, 4)):
            data = np.ones((3, 3, 3), dtype="f4")
            path = os.path.join(self.root, f"cube_{v}.fits")
            writeCube(path, data, irpN=irpN)
            butler.paths[v] = path
        with self.assertRaises(ValueError):
            nirSuperdark.makeMasterDark(
                butler, dict(instrument="PFS", arm="n", spectrograph=1), [0, 1])


@requireDrpStella
class DatasetTypeForIrpTestCase(lsst.utils.tests.TestCase):
    def testMapping(self):
        self.assertEqual(nirSuperdark.datasetTypeForIrp(1), "nirDark")
        self.assertEqual(nirSuperdark.datasetTypeForIrp(4), "nirDark_irp4")
        self.assertEqual(nirSuperdark.datasetTypeForIrp(2), "nirDark_irp2")


@requireDrpStella
class SaveNirDarkTestCase(lsst.utils.tests.TestCase):
    """Round-trip ``saveNirDark`` through a throwaway Gen3 butler."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(prefix="nirSuperdark-repo-")
        repo = os.path.join(cls.root, "repo")
        cls.butler = loadScript("makePfsTestRepo").makePfsTestRepo(repo)

    @classmethod
    def tearDownClass(cls):
        closeButler(cls)
        shutil.rmtree(cls.root, ignore_errors=True)

    START = datetime.datetime(2026, 7, 8, 12, 0, 0)

    def roundTrip(self, irpN, expectedType, saveDir=None, runColl=None,
                  calibCollection=None, endDate=None):
        # A dedicated run per case so the presence check below is unambiguous,
        # and so no two cases put the same dataset into the same run.
        if runColl is None:
            runColl = f"u/test/PIPE2D-1664/{expectedType}"
        self.butler.registry.registerRun(runColl)
        self.runColl = runColl
        if calibCollection is not None:
            nirSuperdark.ensureCalibCollection(self.butler, calibCollection)
        dataId = dict(instrument="PFS", arm="n", spectrograph=1)
        data = np.arange(3 * 4 * 5, dtype="f4").reshape(3, 4, 5)
        rampInfo = dict(W_H4IRPN=irpN, W_H4IRPO=1, W_H4NRED=3, W_H4FRMT=1.5)
        nirSuperdark.saveNirDark(
            self.butler, dataId, runColl, [11, 22, 33], data,
            start=self.START,
            readNoise=5.0, gain=2.5, rampInfo=rampInfo, saveDir=saveDir,
            calibCollection=calibCollection, endDate=endDate)
        cube = self.butler.get(expectedType, dataId, collections=runColl)
        self.assertEqual(cube.nreads, 3)
        self.assertFloatsAlmostEqual(cube.getImageCube(), data)
        # IRP/timing cards recorded for cadence matching
        self.assertEqual(cube.metadata.get("W_H4IRPN"), irpN)
        self.assertEqual(cube.metadata.get("W_H4IRPO"), 1)  # 1..W_H4IRPN, never 0
        self.assertEqual(cube.metadata.get("W_H4FRMT"), 1.5)
        self.assertEqual(cube.metadata.get("W_H4NRED"), 3)
        # calib bookkeeping and dark descriptors
        self.assertEqual(cube.metadata.get("OBSTYPE"), "dark")
        self.assertEqual(cube.metadata.get("GAIN"), 2.5)
        return dataId

    def testNirDarkDefaultType(self):
        self.roundTrip(1, "nirDark")

    def testNirDarkIrp4Type(self):
        dataId = self.roundTrip(4, "nirDark_irp4")
        # An IRP4 combine must not masquerade as the default nirDark.
        self.assertIsNone(
            self.butler.registry.findDataset("nirDark", dataId, collections=self.runColl))

    def testSavedOutsideButlerBeforePut(self):
        """The cube is written to disk so a butler failure cannot discard it."""
        saveDir = tempfile.mkdtemp(prefix="nirSuperdark-save-", dir=self.root)
        self.roundTrip(1, "nirDark", saveDir=saveDir,
                       runColl="u/test/PIPE2D-1664/nirDark-saved")
        path = os.path.join(saveDir, "nirDark-n1.fits")
        self.assertTrue(os.path.exists(path))
        onDisk = ImageCube.readFits(path)
        self.assertEqual(onDisk.nreads, 3)
        self.assertEqual(onDisk.metadata.get("W_H4IRPN"), 1)

    def testCertifiedIntoCalibCollection(self):
        """The dark is certified so ISR can select it by observation date."""
        calib = "u/test/PIPE2D-1664/certified/nirDarkGen.20260709a"
        dataId = self.roundTrip(1, "nirDark", runColl=f"{calib}/put",
                                calibCollection=calib)
        registry = self.butler.registry
        # Valid at and after the start date...
        after = astropy.time.Time(self.START, format="datetime", scale="utc") + 1 * u.day
        self.assertIsNotNone(
            registry.findDataset("nirDark", dataId, collections=calib,
                                 timespan=Timespan(after, after)))
        # ...and open-ended, so a far-future exposure still finds it. The date is
        # past the leap-second table, so silence ERFA's incidental "dubious year".
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", erfa.ErfaWarning)
            far = astropy.time.Time("2099-01-01T00:00:00", format="isot", scale="utc")
        self.assertIsNotNone(
            registry.findDataset("nirDark", dataId, collections=calib,
                                 timespan=Timespan(far, far)))
        # Not valid before it was taken.
        before = astropy.time.Time("2020-01-01T00:00:00", format="isot", scale="utc")
        self.assertIsNone(
            registry.findDataset("nirDark", dataId, collections=calib,
                                 timespan=Timespan(before, before)))

    def testEndDateClosesValidity(self):
        """A regenerated dark for an older run must not stay valid for newer data."""
        calib = "u/test/PIPE2D-1888/closed/nirDark.20260709a"
        end = datetime.datetime(2026, 8, 1)
        dataId = self.roundTrip(1, "nirDark", runColl=f"{calib}/put",
                                calibCollection=calib, endDate=end)
        registry = self.butler.registry

        def at(when):
            t = astropy.time.Time(when, format="datetime", scale="utc")
            return registry.findDataset("nirDark", dataId, collections=calib,
                                        timespan=Timespan(t, t))

        self.assertIsNotNone(at(datetime.datetime(2026, 7, 20)))
        self.assertIsNone(at(datetime.datetime(2026, 8, 2)))


@requireDrpStella
class CollectionNameTestCase(lsst.utils.tests.TestCase):
    def testCalibCollectionDropsGenSuffix(self):
        """Matches the pipeline calibs, e.g. PFS/calib/.../dark.20260503a."""
        self.assertEqual(
            nirSuperdark.calibCollectionName("u/me/calib", "PIPE2D-1664", "run30",
                                             "nirDark_irp4", "20260709a"),
            "u/me/calib/PIPE2D-1664/run30/nirDark_irp4.20260709a")

    def testGenCollectionAndTimestampedRun(self):
        gen = nirSuperdark.genCollectionName("u/me/calib", "PIPE2D-1664", "run30",
                                             "nirDark_irp4", "20260709a")
        self.assertEqual(gen, "u/me/calib/PIPE2D-1664/run30/nirDark_irp4Gen.20260709a")
        self.assertEqual(nirSuperdark.runName(gen, "20260709T123456Z"),
                         f"{gen}/20260709T123456Z")

    def testGenAndCalibDifferOnlyBySuffix(self):
        args = ("u/me/calib", "PIPE2D-1664", "run30", "nirDark", "20260709a")
        self.assertNotEqual(nirSuperdark.calibCollectionName(*args),
                            nirSuperdark.genCollectionName(*args))

    def testTrailingSlashOnBaseStripped(self):
        """A trailing slash must not become a // in the collection name."""
        args = ("PIPE2D-1664", "irp4", "nirDark_irp4", "20260709a")
        self.assertEqual(nirSuperdark.calibCollectionName("u/cpl/calib/", *args),
                         nirSuperdark.calibCollectionName("u/cpl/calib", *args))
        self.assertEqual(nirSuperdark.genCollectionName("u/cpl/calib///", *args),
                         nirSuperdark.genCollectionName("u/cpl/calib", *args))
        self.assertNotIn("//", nirSuperdark.calibCollectionName("u/cpl/calib/", *args))

    def testEmptyBaseRejected(self):
        for bad in ("", "/", "///"):
            with self.assertRaises(ValueError):
                nirSuperdark.collectionRoot(bad)

    def testDefaultIteration(self):
        self.assertRegex(nirSuperdark.defaultIteration(), r"^\d{8}a$")

    def testDefaultTimestampIsUtc(self):
        self.assertRegex(nirSuperdark.defaultTimestamp(), r"^\d{8}T\d{6}Z$")

    def testScratchCollectionLayout(self):
        self.assertEqual(nirSuperdark.scratchCollectionName("u/cpl/calib/", "PIPE2D-1664",
                                                            "irp4"),
                         "u/cpl/calib/PIPE2D-1664/irp4/scratchCubes")

    def testParseScratchCollection(self):
        expected = ("u/cpl/calib", "PIPE2D-1664", "irp4")
        chain = nirSuperdark.scratchCollectionName(*expected)
        self.assertEqual(nirSuperdark.parseScratchCollection(chain), expected)
        # The timestamped RUN pipetask writes inside the chain parses the same.
        self.assertEqual(nirSuperdark.parseScratchCollection(f"{chain}/20260714T063753Z"),
                         expected)

    def testParseScratchCollectionRejectsOthers(self):
        for name in ("u/cpl/calib/PIPE2D-1664/irp4/other",
                     "PIPE2D-1664/irp4/scratchCubes",
                     "u/cpl/calib/PIPE2D-1664/irp4/scratchCubes/a/b"):
            self.assertIsNone(nirSuperdark.parseScratchCollection(name), name)


class FakeRef:
    def __init__(self, spectrograph, visit):
        self.dataId = dict(instrument="PFS", arm="n", spectrograph=spectrograph, visit=visit)


class FakeRegistry:
    """Records the query and returns rawISRCube refs for the given (spec, visit)s."""

    def __init__(self, pairs):
        self.pairs = pairs
        self.queries = []

    def queryDatasets(self, datasetType, collections, where, findFirst):
        self.queries.append((datasetType, collections, where, findFirst))
        return [FakeRef(s, v) for s, v in self.pairs]


class FakeQueryButler:
    def __init__(self, pairs):
        self.registry = FakeRegistry(pairs)


@requireDrpStella
class VisitSelectionTestCase(lsst.utils.tests.TestCase):
    """Which visits each spectrograph combines, given what the inputs hold."""

    def testAvailableVisitsGroupedBySpectrograph(self):
        butler = FakeQueryButler([(1, 12), (1, 10), (2, 10)])
        available = nirSuperdark.availableVisits(butler, "in", [1, 2, 3])
        self.assertEqual(available, {1: [10, 12], 2: [10], 3: []})
        datasetType, collections, where, findFirst = butler.registry.queries[0]
        self.assertEqual((datasetType, collections, findFirst), ("rawISRCube", "in", True))
        self.assertIn("spectrograph IN (1, 2, 3)", where)
        self.assertIn("arm = 'n'", where)

    def testDefaultIsEverythingAvailable(self):
        self.assertEqual(
            nirSuperdark.resolveVisits(None, {1: [10, 11], 2: [10, 12]}),
            {1: [10, 11], 2: [10, 12]})

    def testRequestedVisitsKept(self):
        self.assertEqual(
            nirSuperdark.resolveVisits([10, 11, 11], {1: [10, 11, 12]}), {1: [10, 11]})

    def testMissingVisitsRaiseByDefault(self):
        """The error names the missing visits, so the fix is obvious."""
        with self.assertRaises(ValueError) as cm:
            nirSuperdark.resolveVisits([10, 11, 12], {1: [10, 12], 2: [10, 11, 12]})
        self.assertIn("n1", str(cm.exception))
        self.assertIn("11", str(cm.exception))
        self.assertIn("--skip-missing", str(cm.exception))

    def testMissingVisitsSkippedOnRequest(self):
        with self.assertLogs(nirSuperdark.logger, "WARNING") as cm:
            visits = nirSuperdark.resolveVisits([10, 11, 12], {1: [10, 12], 2: [10, 11, 12]},
                                                skipMissing=True)
        self.assertEqual(visits, {1: [10, 12], 2: [10, 11, 12]})
        self.assertIn("11", "\n".join(cm.output))

    def testTooFewVisitsRaise(self):
        with self.assertRaises(ValueError) as cm:
            nirSuperdark.resolveVisits(None, {1: [10, 11], 2: [10]})
        self.assertIn("n2", str(cm.exception))
        with self.assertRaises(ValueError):
            nirSuperdark.resolveVisits([10, 11], {1: [10]}, skipMissing=True)


@requireDrpStella
class EnsureOutputsTestCase(lsst.utils.tests.TestCase):
    """``ensureOutputs`` must catch, up front, what would fail the final put."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(prefix="nirSuperdark-outputs-")
        repo = os.path.join(cls.root, "repo")
        cls.butler = loadScript("makePfsTestRepo").makePfsTestRepo(repo)

    @classmethod
    def tearDownClass(cls):
        closeButler(cls)
        shutil.rmtree(cls.root, ignore_errors=True)

    def testCollectionsCreated(self):
        calib = "u/test/PIPE2D-1664/fresh/nirDark.20260709a"
        gen = "u/test/PIPE2D-1664/fresh/nirDarkGen.20260709a"
        runColl = nirSuperdark.runName(gen, "20260709T010101Z")
        nirSuperdark.ensureOutputs(self.butler, runColl, calib, "nirDark")
        registry = self.butler.registry
        self.assertEqual(registry.getCollectionType(runColl), CollectionType.RUN)
        self.assertEqual(registry.getCollectionType(calib), CollectionType.CALIBRATION)

    def testChainedCollectionRaises(self):
        chain = "u/test/PIPE2D-1664/chain"
        self.butler.registry.registerCollection(chain, CollectionType.CHAINED)
        with self.assertRaises(ValueError) as cm:
            nirSuperdark.ensureRunCollection(self.butler, chain)
        self.assertIn("CHAINED", str(cm.exception))

    def testGenCollectionChainsRunFirst(self):
        gen = "u/test/PIPE2D-1664/gen/nirDarkGen.20260709a"
        run = f"{gen}/20260709T010101Z"
        for name in (run, "u/test/PIPE2D-1664/inputs"):
            self.butler.registry.registerRun(name)
        nirSuperdark.ensureGenCollection(self.butler, gen, [run, "u/test/PIPE2D-1664/inputs"])
        registry = self.butler.registry
        self.assertEqual(registry.getCollectionType(gen), CollectionType.CHAINED)
        self.assertEqual(list(registry.getCollectionChain(gen)),
                         [run, "u/test/PIPE2D-1664/inputs"])

    def testRerunPrependsNewRunAndKeepsOld(self):
        """A second attempt must not drop the previous run from the chain."""
        gen = "u/test/PIPE2D-1664/rerun/nirDarkGen.20260709a"
        first, second = f"{gen}/20260709T010101Z", f"{gen}/20260709T020202Z"
        for name in (first, second):
            self.butler.registry.registerRun(name)
        nirSuperdark.ensureGenCollection(self.butler, gen, [first])
        nirSuperdark.ensureGenCollection(self.butler, gen, [second])
        self.assertEqual(list(self.butler.registry.getCollectionChain(gen)),
                         [second, first])
        # Idempotent: repeating the second call changes nothing.
        nirSuperdark.ensureGenCollection(self.butler, gen, [second])
        self.assertEqual(list(self.butler.registry.getCollectionChain(gen)),
                         [second, first])

    def testUnknownDatasetTypeRegistered(self):
        """A newly-introduced nirDark_irp<N> need not pre-exist in the repo."""
        registry = self.butler.registry
        datasetType = "nirDark_irp8"
        with self.assertRaises(KeyError):
            registry.getDatasetType(datasetType)
        nirSuperdark.ensureDatasetType(self.butler, datasetType)
        registered = registry.getDatasetType(datasetType)
        self.assertEqual(registered.storageClass_name, "ImageCube")
        self.assertTrue(registered.isCalibration())
        # Idempotent: a second call on an existing repo must not raise.
        nirSuperdark.ensureDatasetType(self.butler, datasetType)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
