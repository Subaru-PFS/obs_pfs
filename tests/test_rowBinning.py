import unittest

import numpy as np

import lsst.utils.tests
import lsst.geom
from lsst.afw.image import MaskedImageF, ExposureF
from lsst.ip.isr import Defects

from lsst.obs.pfs.rowBinning import RowBinning, getRowBinning
from lsst.obs.pfs.isrTask import PfsIsrTask


class RowBinningTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.rng = np.random.RandomState(12345)
        self.factor = 5
        self.height = 23  # Not a multiple of factor
        self.width = 7
        self.binning = RowBinning(self.factor, self.height)

    def testGetRowBinning(self):
        """Test reading the binning from the header"""
        self.assertEqual(getRowBinning(dict(W_CDRBIN=20)), 20)
        self.assertEqual(getRowBinning(dict(W_CDRBIN=1)), 1)
        self.assertEqual(getRowBinning(dict(W_CDRBIN=0)), 1)
        self.assertEqual(getRowBinning(dict()), 1)
        self.assertEqual(getRowBinning(None), 1)
        self.assertFalse(RowBinning.fromMetadata(dict(), self.height))
        self.assertTrue(RowBinning.fromMetadata(dict(W_CDRBIN=2), self.height))

    def testGeometry(self):
        """Test the mapping between binned and unbinned rows"""
        binning = self.binning
        self.assertEqual(binning.numBinned, 5)
        self.assertFloatsEqual(binning.counts, [5, 5, 5, 5, 3])
        self.assertFloatsEqual(binning.starts, [0, 5, 10, 15, 20])
        self.assertFloatsEqual(binning.centers, [2, 7, 12, 17, 21])
        self.assertFloatsEqual(binning.index, np.arange(self.height)//self.factor)
        yBinned = np.linspace(-0.5, 4.5, 11)
        self.assertFloatsAlmostEqual(binning.toBinned(binning.toUnbinned(yBinned)), yBinned, atol=1.0e-12)
        # Edges of binned pixels map to edges of unbinned pixels
        self.assertFloatsAlmostEqual(binning.toUnbinned([-0.5, 0.5]), [-0.5, self.factor - 0.5], atol=1.0e-12)

    def testArrays(self):
        """Test binning of arrays"""
        binning = self.binning
        array = self.rng.normal(size=(self.height, self.width))
        summed = binning.sum(array)
        self.assertEqual(summed.shape, (binning.numBinned, self.width))
        for jj in range(binning.numBinned):
            rows = slice(binning.starts[jj], binning.starts[jj] + binning.counts[jj])
            self.assertFloatsAlmostEqual(summed[jj], array[rows].sum(axis=0), atol=1.0e-12)
        self.assertFloatsAlmostEqual(binning.mean(array), summed/binning.counts[:, np.newaxis], atol=1.0e-12)
        self.assertFloatsAlmostEqual(binning.mean(array.T, axis=1), binning.mean(array).T, atol=1.0e-12)

        mask = np.zeros((self.height, self.width), dtype=np.int32)
        mask[6, 3] = 0x2
        mask[8, 3] = 0x4
        orred = binning.bitwiseOr(mask)
        self.assertEqual(orred[1, 3], 0x6)
        self.assertEqual(np.sum(orred != 0), 1)

        replicated = binning.replicate(array, "mean")
        self.assertEqual(replicated.shape, array.shape)
        self.assertFloatsAlmostEqual(binning.mean(replicated), binning.mean(array), atol=1.0e-12)
        for jj in range(binning.numBinned):
            rows = slice(binning.starts[jj], binning.starts[jj] + binning.counts[jj])
            self.assertTrue(np.all(replicated[rows] == replicated[binning.starts[jj]]))

        with self.assertRaises(RuntimeError):
            binning.sum(array[1:])
        with self.assertRaises(ValueError):
            binning.bin(array, "median")

    def makeMaskedImage(self):
        image = MaskedImageF(self.width, self.height)
        image.image.array[:] = self.rng.normal(size=(self.height, self.width))
        image.variance.array[:] = self.rng.uniform(1, 2, size=(self.height, self.width))
        image.mask.array[:] = 0
        image.mask.array[7, 2] = image.mask.getPlaneBitMask("BAD")
        return image

    def testBinMaskedImage(self):
        """Test collapsing an image to one row per binned row"""
        binning = self.binning
        image = self.makeMaskedImage()
        image.setXY0(lsst.geom.Point2I(3, 0))
        binned = binning.binMaskedImage(image)
        self.assertEqual(binned.getDimensions(), lsst.geom.Extent2I(self.width, binning.numBinned))
        self.assertEqual(binned.getXY0(), image.getXY0())
        self.assertFloatsAlmostEqual(binned.image.array, binning.mean(image.image.array), rtol=1.0e-6)
        self.assertFloatsAlmostEqual(binned.variance.array, binning.mean(image.variance.array), rtol=1.0e-6)
        self.assertEqual(binned.mask.array[1, 2], image.mask.getPlaneBitMask("BAD"))
        self.assertEqual(np.sum(binned.mask.array != 0), 1)

        summed = binning.binMaskedImage(image, "sum")
        self.assertFloatsAlmostEqual(summed.image.array, binning.sum(image.image.array), rtol=1.0e-5)

        exposure = ExposureF(image)
        exposure.getMetadata().set("W_CDRBIN", self.factor)
        binnedExp = binning.binExposure(exposure)
        self.assertEqual(binnedExp.getHeight(), binning.numBinned)
        self.assertEqual(binnedExp.getMetadata().get("W_CDRBIN"), self.factor)

    def testReplicateMaskedImage(self):
        """Test binning a calibration while keeping it full-sized"""
        binning = self.binning
        for how in ("mean", "sum"):
            image = self.makeMaskedImage()
            original = image.clone()
            binning.replicateMaskedImage(image, how)
            self.assertEqual(image.getDimensions(), original.getDimensions())
            expect = binning.bin(original.image.array, how)
            self.assertFloatsAlmostEqual(image.image.array[binning.starts], expect, rtol=1.0e-5)
            if how == "mean":
                expectVar = binning.sum(original.variance.array)/binning.counts[:, np.newaxis]**2
            else:
                expectVar = binning.sum(original.variance.array)
            self.assertFloatsAlmostEqual(image.variance.array[binning.starts], expectVar, rtol=1.0e-5)
            bad = image.mask.getPlaneBitMask("BAD")
            self.assertTrue(np.all(image.mask.array[5:10, 2] == bad))
            self.assertEqual(np.sum(image.mask.array != 0), self.factor)

    def testBinDefects(self):
        """Test growing defects to cover entire groups of binned rows"""
        binning = self.binning
        defects = Defects([
            lsst.geom.Box2I(lsst.geom.Point2I(1, 6), lsst.geom.Point2I(1, 6)),
            lsst.geom.Box2I(lsst.geom.Point2I(3, 4), lsst.geom.Point2I(4, 11)),
            lsst.geom.Box2I(lsst.geom.Point2I(5, 21), lsst.geom.Point2I(5, 22)),
        ])
        task = PfsIsrTask.__new__(PfsIsrTask)  # binDefects doesn't need the task to be configured
        binned = task.binDefects(defects, binning)
        expect = np.zeros((self.height, self.width), dtype=bool)
        expect[5:10, 1] = True
        expect[0:15, 3:5] = True
        expect[20:23, 5] = True
        got = np.zeros_like(expect)
        for defect in binned:
            box = defect.getBBox()
            got[box.getMinY():box.getMaxY() + 1, box.getMinX():box.getMaxX() + 1] = True
        self.assertTrue(np.array_equal(got, expect))


class BinnedOverscanTestCase(lsst.utils.tests.TestCase):
    def testFitOverscanBinned(self):
        """Test that the overscan fit is constant within each binned group"""
        from lsst.obs.pfs.overscan import PfsOverscanCorrectionTask
        rng = np.random.RandomState(12345)
        config = PfsOverscanCorrectionTask.ConfigClass()
        config.fitType = "AKIMA_SPLINE"
        config.order = 10
        task = PfsOverscanCorrectionTask(config=config)

        binning = 20
        firstBinnedRow = 48
        height, width = 4300, 30
        rows = np.arange(height)
        group = np.where(rows < firstBinnedRow, rows,
                         firstBinnedRow + binning*((rows - firstBinnedRow)//binning))
        _, first, inverse = np.unique(group, return_index=True, return_inverse=True)
        perRead = 1000.0 + 0.001*first[:, np.newaxis] + rng.normal(scale=1.0, size=(len(first), width))
        image = MaskedImageF(width, height)
        image.image.array[:] = perRead[inverse]  # each read repeated, as the DAQ writes binned data
        image.mask.array[:] = 0
        image.variance.array[:] = 1.0

        result = task.fitOverscanBinned(image, False, binning=binning, firstBinnedRow=firstBinnedRow,
                                        fitOverscan=task.fitOverscan)
        value = np.asarray(result.overscanValue)
        self.assertEqual(value.shape, (height,))
        for start in range(firstBinnedRow, height, binning):
            self.assertTrue(np.all(value[start:start + binning] == value[start]))
        self.assertFloatsAlmostEqual(value, 1000.0 + 0.001*first[inverse], atol=0.5)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
