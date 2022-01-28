from lsst.obs.pfs.utils import getLamps, getLampElements
from lsst.daf.base import PropertyList
import lsst.utils.tests


class ReadLampsTestCase(lsst.utils.tests.TestCase):
    """Tests the reading of lamp metadata"""

    def setUp(self):
        pass

    def createLampMetadata(self, activeKeys):
        """Creates a dictionary of lamp keys and boolean values (on or off).
        All lamp are set to 'off', with
        the exception of the lamps corresponding to the input keys.

        Parameters
        ----------
        activeKeys: `list`[`str`]
            header keywords corresponding to the active lamps

        Returns
        -------
        metadata : `dict`[`bool`]
            Header metadata.
        """
        metadata = {
            "W_AITNEO": False,
            "W_AITXEN": False,
            "W_AITHGA": False,
            "W_AITKRY": False,
            "W_AITARG": False,
            "W_AITHGC": False,
            "W_AITQTH": False,
            "W_AITDEU": False,
        }

        for key in activeKeys:
            metadata[key] = True

        return metadata

    def testGetLamps(self):
        """Tests that metadata reflects the active lamps"""
        md = self.createLampMetadata(['W_AITHGC'])
        self.assertEqual(getLamps(md), {'HgCd'})

        md = self.createLampMetadata(['W_AITARG'])
        self.assertEqual(getLamps(md), {'Ar'})

        md = self.createLampMetadata(['W_AITQTH'])
        self.assertEqual(getLamps(md), {'Quartz'})

        md = self.createLampMetadata(['W_AITARG', 'W_AITKRY'])
        lampSet = getLamps(md)
        for lamp in ['Ar', 'Kr']:
            self.assertTrue(lamp in lampSet,
                            f'Lamp {lamp} is not in active lamp set {lampSet}')

    def checkLampElements(self, keys, expectedElements):
        """Checks that expected lamp elements are retrieved
        for the supplied lamp header keys.

        Parameters
        ----------
        keys: `list`[`str`]
            header keywords corresponding to the available lamps
        expectedElements: `set`[`str`]
            The expected lamp elements
        """

        md = self.createLampMetadata(keys)
        pList = PropertyList()
        for key, value in md.items():
            pList.addBool(key, value)
        actualElements = set(getLampElements(pList))
        print(f'{expectedElements}, {actualElements}')
        print(expectedElements.difference(actualElements))
        self.assertTrue(actualElements == expectedElements,
                        f'for keys {keys} '
                        f'actual elements found: {actualElements}, '
                        f'expected: {expectedElements}')

    def testGetLampElements(self):
        """Tests that metadata reflects the active lamp elements"""

        self.checkLampElements(['W_AITHGC'], {'Hg', 'Cd', 'Ar'})
        self.checkLampElements(['W_AITHGA'], {'Hg', 'Ar'})
        self.checkLampElements(['W_AITKRY'], {'Kr'})
        self.checkLampElements(['W_AITNEO'], {'Ne'})
