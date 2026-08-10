"""The ``h4.crMinReads`` knob is deprecated and ignored, but the field is
retained so configs that still set it (e.g. an older drp_stella) load without
error. Setting it must not raise, must not change behaviour, and must emit a
noticeable ``FutureWarning`` (via the field's ``deprecated=`` marker).
"""
from __future__ import annotations

import unittest
import warnings

import lsst.utils.tests

from testUtils import HAS_DRP_STELLA, requireDrpStella

if HAS_DRP_STELLA:
    # isrTask imports pfs.drp.stella.crosstalk.
    from lsst.obs.pfs import isrTask as pfsIsrTask


def _baseConfig():
    config = pfsIsrTask.PfsIsrTask.ConfigClass()
    config.doFlat = False
    config.doDark = False
    config.doDefect = False
    config.doSaturationInterpolation = False
    config.h4.quickCDS = False
    config.h4.doIPC = False
    config.h4.doWriteRawCube = False
    config.h4.doLinearize = False
    return config


@requireDrpStella
class CrMinReadsDeprecationTestCase(lsst.utils.tests.TestCase):
    def test_settingWarnsAndDoesNotRaise(self):
        # The field must exist again (drp_stella compat): setting it emits a
        # FutureWarning, validate() still passes, and the value is preserved.
        config = _baseConfig()
        with self.assertWarns(FutureWarning) as cm:
            config.h4.crMinReads = 4
        self.assertIn("crMinReads", str(cm.warning))
        config.validate()
        self.assertEqual(config.h4.crMinReads, 4)

    def test_defaultIsUnset(self):
        # Left alone the knob is unset, and validating a default config emits no
        # deprecation warning.
        config = _baseConfig()
        self.assertIsNone(config.h4.crMinReads)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            config.validate()
        self.assertFalse([w for w in caught if issubclass(w.category, FutureWarning)])

    def test_taskBuildsWithKnobSet(self):
        # A config that sets the deprecated knob must still build the task.
        config = _baseConfig()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            config.h4.crMinReads = 4
        config.validate()
        pfsIsrTask.PfsIsrTask(config=config)  # must not raise


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(__import__("sys").modules["__main__"])
    unittest.main(failfast=True)
