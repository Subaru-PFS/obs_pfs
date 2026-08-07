import sys
import unittest

import lsst.utils.tests

from lsst.obs.pfs.translator import PfsTranslator


class TranslatorTestCase(lsst.utils.tests.TestCase):
    """Tests for the PFS metadata translator"""

    def testIrpRatio(self):
        """W_H4IRPN is translated to ext_irp, and is None when absent"""
        self.assertIn("irp", PfsTranslator.extensions)

        translator = PfsTranslator({"INSTRUME": "PFS", "W_H4IRPN": 4})
        self.assertEqual(translator.to_ext_irp(), 4)

        translator = PfsTranslator({"INSTRUME": "PFS", "W_H4IRPN": 1})
        self.assertEqual(translator.to_ext_irp(), 1)

        # CCD exposures lack the card and must yield None (no default).
        translator = PfsTranslator({"INSTRUME": "PFS"})
        self.assertIsNone(translator.to_ext_irp())


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
