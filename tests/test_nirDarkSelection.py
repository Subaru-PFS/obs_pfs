"""Tests for the NIR dark ratio selection.

Because the registry does not (yet) carry the IRP ratio dimension, the choice
between ``nirDark`` (ratio 1) and ``nirDark_irp4`` (ratio 4) cannot be made at
qgraph-build time. Instead:

- ``lookupNirDark`` resolves *both* dataset types whenever their datasets
  exist (it does not consult ``irp_ratio``), and
- ``selectNirDark`` runs during runtime initialization and picks the
  ratio-matched dark from the loaded ramp's ``W_H4IRPN`` (``PfsRaw.irpN``),
  handing ``run()`` a single, ratio-agnostic ``nirDark``.
"""
import unittest

import lsst.utils.tests
from lsst.obs.pfs import isrTask as pfsIsrTask


class _DatasetType:
    def __init__(self, name):
        self.name = name


class _DataId:
    """A dataId whose visit record carries no ``irp_ratio`` (registry lacks
    the dimension)."""

    def __init__(self):
        self.records = {"visit": object()}
        self.timespan = None

    def hasRecords(self):
        return True


class _Registry:
    def __init__(self, ref):
        self._ref = ref

    def findDataset(self, datasetType, collections, dataId, timespan):
        return self._ref


class LookupNirDarkTestCase(lsst.utils.tests.TestCase):
    def testResolvesIrp4WithoutIrpRatio(self):
        # The regression case: no irp_ratio in the record, yet nirDark_irp4
        # must still resolve so the runtime picker can choose it.
        ref = object()
        out = pfsIsrTask.lookupNirDark(
            _DatasetType("nirDark_irp4"), _Registry(ref), _DataId(), ["c"])
        self.assertEqual(out, [ref])

    def testResolvesNirDark(self):
        ref = object()
        out = pfsIsrTask.lookupNirDark(
            _DatasetType("nirDark"), _Registry(ref), _DataId(), ["c"])
        self.assertEqual(out, [ref])

    def testEmptyWhenNoDataset(self):
        out = pfsIsrTask.lookupNirDark(
            _DatasetType("nirDark_irp4"), _Registry(None), _DataId(), ["c"])
        self.assertEqual(out, [])


class _Raw:
    def __init__(self, isNir, irpN):
        self._isNir = isNir
        self._irpN = irpN

    def isNir(self):
        return self._isNir

    @property
    def irpN(self):
        return self._irpN


class SelectNirDarkTestCase(lsst.utils.tests.TestCase):
    def testPicksIrp4ForRatio4(self):
        nir1, nir4 = object(), object()
        inputs = {"ccdExposure": _Raw(True, 4), "nirDark": nir1, "nirDarkIrp4": nir4}
        pfsIsrTask.selectNirDark(inputs)
        self.assertIs(inputs["nirDark"], nir4)
        self.assertNotIn("nirDarkIrp4", inputs)

    def testPicksIrp1ForRatio1(self):
        nir1, nir4 = object(), object()
        inputs = {"ccdExposure": _Raw(True, 1), "nirDark": nir1, "nirDarkIrp4": nir4}
        pfsIsrTask.selectNirDark(inputs)
        self.assertIs(inputs["nirDark"], nir1)
        self.assertNotIn("nirDarkIrp4", inputs)

    def testIrp4RampWithoutIrp4DarkYieldsNone(self):
        # Only the irp1 dark resolved: the irp4 ramp must NOT silently fall
        # back to it (that is the wrong-cadence bug); it gets None, which the
        # downstream "No NIR dark cube found" check turns into a clear error.
        nir1 = object()
        inputs = {"ccdExposure": _Raw(True, 4), "nirDark": nir1, "nirDarkIrp4": None}
        pfsIsrTask.selectNirDark(inputs)
        self.assertIsNone(inputs["nirDark"])

    def testCcdLeavesNoNirDark(self):
        inputs = {"ccdExposure": _Raw(False, 1)}
        pfsIsrTask.selectNirDark(inputs)
        self.assertNotIn("nirDark", inputs)
        self.assertNotIn("nirDarkIrp4", inputs)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
