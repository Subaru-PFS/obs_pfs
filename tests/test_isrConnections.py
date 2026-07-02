"""Connection-contract tests for ``PfsIsrTask``.

The NIR calibrations (``h4Linearity``, ``badRefPixels``) must be *optional*
prerequisite inputs (``minimum=0``). That is what lets data be reduced when no
such calibration is available for the observation date / detector — e.g. older
n2 data taken with the previous detector, for which no h4Linearity solution
exists. In that case the quantum graph still builds, ISR warns, and reduction
proceeds without the correction rather than stopping.
"""
import sys
import unittest

import lsst.utils.tests
from lsst.obs.pfs.isrTask import PfsIsrConnections, PfsIsrTask


class PfsIsrConnectionsTestCase(lsst.utils.tests.TestCase):
    def testLinearityIsOptional(self):
        """With the default config (h4.doLinearize=True) the ``h4Linearity``
        input is present but optional, so a missing calib does not block the
        quantum graph."""
        config = PfsIsrTask.ConfigClass()
        self.assertTrue(config.h4.doLinearize)
        connections = PfsIsrConnections(config=config)
        self.assertIn("linearity", connections.prerequisiteInputs)
        self.assertEqual(connections.linearity.minimum, 0)

    def testLinearityRemovedWhenDisabled(self):
        """The input is dropped entirely when linearization is disabled."""
        config = PfsIsrTask.ConfigClass()
        config.h4.doLinearize = False
        connections = PfsIsrConnections(config=config)
        self.assertNotIn("linearity", connections.prerequisiteInputs)

    def testBadRefPixelsIsOptional(self):
        """``badRefPixels`` is likewise optional under the default config."""
        config = PfsIsrTask.ConfigClass()
        self.assertTrue(config.h4.doIRPbadPixels)
        connections = PfsIsrConnections(config=config)
        self.assertIn("badRefPixels", connections.prerequisiteInputs)
        self.assertEqual(connections.badRefPixels.minimum, 0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
