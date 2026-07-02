"""Integration tests for the script-installed NIR calibration products.

Builds a throwaway butler repository, exercises the ``makePfsCalibs`` install
logic, and confirms that validity-date selection and CALIBRATION run chaining
work for the ``defects``, ``nirLinearity`` and ``badRefPixels`` products.
"""

import importlib.util
import os
import shutil
import sys
import tempfile
import unittest

import astropy.time

from lsst.daf.base import PropertyList
from lsst.daf.butler import CollectionType, Timespan
from lsst.obs.pfs.nirLinearity import NirLinearity
import lsst.utils.tests

OBS_PFS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def loadScript(name):
    """Import one of the ``bin.src`` scripts as a module."""
    path = os.path.join(OBS_PFS_DIR, "bin.src", f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def utc(iso):
    return astropy.time.Time(iso, format="isot", scale="utc")


class MakePfsCalibsTestCase(lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.install = loadScript("makePfsCalibs")
        cls.root = tempfile.mkdtemp(prefix="pfsCalibs-")
        repo = os.path.join(cls.root, "repo")
        cls.butler = loadScript("makePfsTestRepo").makePfsTestRepo(repo)
        cls.drpData = cls.install.getPackageDir("drp_pfs_data")

        cls.collection = "u/test/PIPE2D-1856/test"
        cls.defectsColl = f"{cls.collection}/defectsGen.99999999a"
        cls.badRefColl = f"{cls.collection}/badRefPixelsGen.99999999a"
        for coll in (cls.defectsColl, cls.badRefColl):
            cls.butler.registry.registerCollection(coll, CollectionType.CALIBRATION)

        cls.install.installDefects(cls.butler, cls.drpData, cls.defectsColl, "PFS")
        cls.install.installBadRefPixels(
            cls.butler, cls.drpData, cls.badRefColl, "PFS", begin=None, end=None
        )

        # Synthetic nirLinearity: two versions for n2 split at the boundary T.
        cls.linColl = f"{cls.collection}/nirLinearityGen.99999999a"
        cls.butler.registry.registerCollection(cls.linColl, CollectionType.CALIBRATION)
        cls.boundary = utc("2025-02-15T00:00:00")
        dataId = dict(instrument="PFS", arm="n", spectrograph=2)
        for label, timespan in (
            ("v1", Timespan(None, cls.boundary)),
            ("v2", Timespan(cls.boundary, None)),
        ):
            metadata = PropertyList()
            metadata.set("VERSION", label)
            calib = NirLinearity.empty((4, 4), 3, metadata)
            run = f"{cls.linColl}/{label}"
            cls.butler.registry.registerRun(run)
            ref = cls.butler.put(calib, "nirLinearity", dataId, run=run)
            cls.butler.registry.certify(cls.linColl, [ref], timespan)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)

    def findVersion(self, datasetType, dataId, collection, when):
        timespan = Timespan(utc(when), utc(when))
        ref = self.butler.registry.findDataset(
            datasetType, dataId, collections=collection, timespan=timespan
        )
        return None if ref is None else os.path.basename(ref.run)

    def testDefectsDateSelection(self):
        """n2 defects (versions 1970, 2023-07-01, 2025-02-15) select by date."""
        dataId = dict(instrument="PFS", detector=5)  # n2
        self.assertEqual(
            self.findVersion("defects", dataId, self.defectsColl, "2019-01-01T00:00:00"),
            "19700101T000000",
        )
        self.assertEqual(
            self.findVersion("defects", dataId, self.defectsColl, "2024-01-01T00:00:00"),
            "20230701T000000",
        )
        self.assertEqual(
            self.findVersion("defects", dataId, self.defectsColl, "2025-06-01T00:00:00"),
            "20250215T000000",
        )

    def testBadRefPixelsReadback(self):
        """badRefPixels round-trips through the butler, incl. the empty n2."""
        now = Timespan(utc("2026-01-01T00:00:00"), utc("2026-01-01T00:00:00"))
        n1 = self.butler.get(
            "badRefPixels", dict(instrument="PFS", arm="n", spectrograph=1),
            collections=self.badRefColl, timespan=now,
        )
        self.assertGreater(n1.pixels.size, 0)
        self.assertEqual(n1.metadata.get("DETNAME"), "n1")
        n2 = self.butler.get(
            "badRefPixels", dict(instrument="PFS", arm="n", spectrograph=2),
            collections=self.badRefColl, timespan=now,
        )
        self.assertEqual(n2.pixels.size, 0)

    def testNirLinearityDateSelection(self):
        """The two synthetic nirLinearity versions select across the boundary."""
        dataId = dict(instrument="PFS", arm="n", spectrograph=2)
        self.assertEqual(
            self.findVersion("nirLinearity", dataId, self.linColl, "2024-01-01T00:00:00"), "v1"
        )
        self.assertEqual(
            self.findVersion("nirLinearity", dataId, self.linColl, "2025-06-01T00:00:00"), "v2"
        )


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
