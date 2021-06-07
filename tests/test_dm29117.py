import sys
import unittest

from lsst.obs.pfs import PfsMapper  # noqa: imported for side-effect
from lsst.daf.base import PropertyList

import lsst.utils.tests


class DM29117TestCase(lsst.utils.tests.TestCase):
    """A test case for DM-29117

    DM-29117 (not fixed in LSST 18.1.0) prevents overwriting a PropertyList
    entry with an entry of a different type. This prevents us from fixing wrong
    pfsDesignId values when they would need to be promoted from int to long.
    """
    def testFixed(self):
        """Test that the DM-29117 bug is fixed"""
        name = "foo"
        small = 123  # An int
        big = 0xfeedfacedeadbeef  # Not an int
        pl = PropertyList()
        pl.set(name, small)
        self.assertEqual(pl.get(name), small)
        pl.set(name, big)  # without the DM-29117 fix, would get an exception here
        self.assertEqual(pl.get(name), big)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
