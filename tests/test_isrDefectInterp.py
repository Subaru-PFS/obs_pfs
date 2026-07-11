"""``PfsIsrTask`` must honor ``config.doInterpolate`` for defect handling.

The H4 NIR path masks defects as BAD and, when ``doInterpolate`` is True,
interpolates over them (stamping INTRP). With ``doInterpolate`` False it must
mask them but leave the pixel values untouched — the base
``maskAndInterpolateDefects`` interpolates unconditionally, which this checks
is no longer the behavior.
"""
import unittest

import lsst.utils.tests
import lsst.geom as geom
import lsst.afw.image as afwImage
from lsst.ip.isr import Defects
from lsst.obs.pfs import isrTask as pfsIsrTask


def _makeIsrTask(doInterpolate):
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doInterpolate = doInterpolate
    config.doDark = True
    config.doFlat = False
    config.doDefect = True
    config.doSaturationInterpolation = False
    config.h4.quickCDS = False
    config.h4.doIPC = False
    config.h4.doWriteRawCube = True
    config.h4.doLinearize = False
    config.h4.doCR = False
    config.validate()
    return pfsIsrTask.PfsIsrTask(config=config)


def _makeExposureAndDefects():
    exp = afwImage.ExposureF(20, 20)
    exp.image.array[:] = 10.0
    exp.image.array[10, 10] = 9999.0            # the defect pixel's bad value
    box = geom.Box2I(geom.Point2I(10, 10), geom.Extent2I(1, 1))
    return exp, Defects([box])


class DefectInterpolationTestCase(lsst.utils.tests.TestCase):
    def testInterpolatesWhenEnabled(self):
        task = _makeIsrTask(doInterpolate=True)
        exp, defects = _makeExposureAndDefects()
        task._maskDefects(exp, defects)
        m = exp.mask
        self.assertTrue(m.array[10, 10] & m.getPlaneBitMask("BAD"))
        self.assertTrue(m.array[10, 10] & m.getPlaneBitMask("INTRP"))
        self.assertNotAlmostEqual(float(exp.image.array[10, 10]), 9999.0)

    def testMasksWithoutInterpolatingWhenDisabled(self):
        task = _makeIsrTask(doInterpolate=False)
        exp, defects = _makeExposureAndDefects()
        task._maskDefects(exp, defects)
        m = exp.mask
        self.assertTrue(m.array[10, 10] & m.getPlaneBitMask("BAD"))      # still masked
        self.assertFalse(m.array[10, 10] & m.getPlaneBitMask("INTRP"))   # not interpolated
        self.assertEqual(float(exp.image.array[10, 10]), 9999.0)         # value untouched


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
