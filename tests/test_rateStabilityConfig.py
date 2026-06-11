import unittest

import lsst.utils.tests
from lsst.obs.pfs.isrTask import H4Config


class H4ConfigRateStabilityTestCase(lsst.utils.tests.TestCase):
    """The rate-stability fields exist on H4Config with the documented
    defaults and accept overrides."""

    def testDefaults(self):
        cfg = H4Config()
        self.assertTrue(cfg.doRateStability)
        self.assertAlmostEqual(cfg.rateStabilityThreshold, 0.20)
        self.assertAlmostEqual(cfg.rateStabilityRateFloorADU, 5.0)
        self.assertEqual(cfg.rateStabilityMinDeltasPerSegment, 3)
        self.assertFalse(cfg.maskRateUnstable)

    def testFieldsAreSettable(self):
        cfg = H4Config()
        cfg.doRateStability = False
        cfg.rateStabilityThreshold = 0.30
        cfg.rateStabilityRateFloorADU = 10.0
        cfg.rateStabilityMinDeltasPerSegment = 5
        cfg.maskRateUnstable = True
        cfg.validate()


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
