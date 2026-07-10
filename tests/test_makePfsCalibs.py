"""Integration tests for the script-installed calibration products.

Builds a throwaway butler repository, exercises the ``makePfsCalibs`` install
logic, and confirms that validity-date selection and CALIBRATION run chaining
work for the ``defects``, ``badRefPixels`` and ``h4Linearity`` products.
"""

import importlib.util
import os
import shutil
import sys
import tempfile
import unittest

import astropy.time
import numpy as np

from lsst.daf.butler import CollectionType, Timespan
from lsst.obs.pfs.h4Linearity.models.polynomial import PolynomialModel
from lsst.obs.pfs.h4Linearity.types import Diagnostics, LinearityCorrection
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


def makeCorrection(value, shape=(4, 4)):
    """Build a tiny degree-1 h4Linearity `LinearityCorrection`.

    The first-order coefficient is filled with ``value`` so distinct
    versions are trivially distinguishable after a round-trip.
    """
    height, width = shape
    coefficients = np.zeros((2, height, width), dtype=np.float32)
    coefficients[1] = value
    diagnostics = Diagnostics(
        residualRms=np.zeros(shape, dtype=np.float32),
        maxAbsResidual=np.zeros(shape, dtype=np.float32),
        nPointsUsed=np.full(shape, 11, dtype=np.int32),
        monotonic=np.ones(shape, dtype=bool),
        conditionNumber=np.ones(shape, dtype=np.float32),
        summary={},
    )
    return LinearityCorrection(
        model=PolynomialModel(order=1),
        coefficients=coefficients,
        fitMin=np.zeros(shape, dtype=np.float32),
        fitMax=np.full(shape, 100.0, dtype=np.float32),
        badPixelMask=np.zeros(shape, dtype=np.uint8),
        diagnostics=diagnostics,
    )


@unittest.skipIf("DRP_PFS_DATA_DIR" not in os.environ, "drp_pfs_data is not setup")
class MakePfsCalibsTestCase(lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.install = loadScript("makePfsCalibs")
        cls.root = tempfile.mkdtemp(prefix="pfsCalibs-")
        repo = os.path.join(cls.root, "repo")
        cls.butler = loadScript("makePfsTestRepo").makePfsTestRepo(repo)
        cls.drpData = cls.install.getPackageDir("drp_pfs_data")

        cls.collection = "u/test/PIPE2D-1856/test"
        cls.iteration = "99999999a"
        cls.timestamp = "99999999T000000Z"
        # DMTN-222 layout: the pipeline selects from ``{product}.{iteration}``
        # (CALIBRATION); the runs live under ``{product}Gen.{iteration}`` (CHAINED).
        cls.defectsColl = f"{cls.collection}/defects.{cls.iteration}"
        cls.defectsGen = f"{cls.collection}/defectsGen.{cls.iteration}"
        cls.badRefColl = f"{cls.collection}/badRefPixels.{cls.iteration}"
        cls.badRefGen = f"{cls.collection}/badRefPixelsGen.{cls.iteration}"
        for coll in (cls.defectsColl, cls.badRefColl):
            cls.butler.registry.registerCollection(coll, CollectionType.CALIBRATION)

        cls.install.installDefects(cls.butler, cls.drpData, cls.defectsColl, cls.defectsGen, "PFS")
        cls.install.installBadRefPixels(
            cls.butler, cls.drpData, cls.badRefColl, cls.badRefGen, cls.timestamp,
            "PFS", begin=None, end=None,
        )

        # h4Linearity: two synthetic versions for n2 split at the boundary T.
        cls.h4LinColl = f"{cls.collection}/h4Linearity.{cls.iteration}"
        cls.butler.registry.registerCollection(cls.h4LinColl, CollectionType.CALIBRATION)
        cls.boundary = utc("2025-02-15T00:00:00")
        cls.corrValues = {"v1": 1.0, "v2": 2.0}
        dataId = dict(instrument="PFS", arm="n", spectrograph=2)
        for label, timespan in (
            ("v1", Timespan(None, cls.boundary)),
            ("v2", Timespan(cls.boundary, None)),
        ):
            run = f"{cls.h4LinColl}/{label}"
            cls.butler.registry.registerRun(run)
            ref = cls.butler.put(makeCorrection(cls.corrValues[label]), "h4Linearity", dataId, run=run)
            cls.butler.registry.certify(cls.h4LinColl, [ref], timespan)

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

    def testH4LinearityDateSelection(self):
        """The two h4Linearity versions select across the boundary."""
        dataId = dict(instrument="PFS", arm="n", spectrograph=2)
        self.assertEqual(
            self.findVersion("h4Linearity", dataId, self.h4LinColl, "2024-01-01T00:00:00"), "v1"
        )
        self.assertEqual(
            self.findVersion("h4Linearity", dataId, self.h4LinColl, "2025-06-01T00:00:00"), "v2"
        )

    def testH4LinearityRoundTrip(self):
        """h4Linearity round-trips through the H4Linearity storage class."""
        now = Timespan(utc("2025-06-01T00:00:00"), utc("2025-06-01T00:00:00"))
        calib = self.butler.get(
            "h4Linearity", dict(instrument="PFS", arm="n", spectrograph=2),
            collections=self.h4LinColl, timespan=now,
        )
        self.assertIsInstance(calib, LinearityCorrection)
        self.assertFloatsEqual(calib.coefficients[1], self.corrValues["v2"])

    def testDatasetTypeSpecsMatchRegistry(self):
        """The ``--register-dataset-types`` specs must match what the
        instrument's ``register()`` produced (so they cannot drift)."""
        for product, (storageClass, dimensions) in self.install.DATASET_TYPES.items():
            datasetType = self.butler.registry.getDatasetType(product)
            self.assertEqual(datasetType.storageClass_name, storageClass)
            self.assertEqual(set(datasetType.dimensions.names), set(dimensions))
            self.assertTrue(datasetType.isCalibration)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
