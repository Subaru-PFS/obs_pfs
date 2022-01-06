"""
Utilities for working with PFS data
"""


def getLamps(md):
    """Return a list of lit lamps, according to the header

    Parameters
    ----------
    md : `dict`
        Header metadata.

    Returns
    -------
    lamps : `set` of `str`
        Set of lit lamps.
    """
    menu = {
        # Modern "short" headers
        "W_AITNEO": "Ne",
        "W_AITXEN": "Xe",
        "W_AITHGA": "HgAr",
        "W_AITKRY": "Kr",
        "W_AITARG": "Ar",
        "W_AITHGC": "HgCd",
        "W_AITQTH": "Quartz",
        "W_AITDEU": "D",
        # Old "long" headers (2017)
        "W_AIT_SRC_Ne": "Ne",
        "W_AIT_SRC_HgAr": "HgAr",
        "W_AIT_SRC_Xe": "Xe",
        "W_AIT_SRC_Qth": "Quartz",
    }
    return {menu[key] for key in menu if md.get(key, False)}


def getLampElements(md):
    """Return a set of the elements found in the lamps that are on

    @param md: dafBase.PropertyList containing the header
    """

    md = md.toDict()
    menu = {
        # Modern "short" headers
        "W_AITNEO": ["Ne"],
        "W_AITXEN": ["Xe"],
        "W_AITHGA": ["Hg", "Ar"],
        "W_AITKRY": ["Kr"],
        "W_AITARG": ["Ar"],
        "W_AITHGC": ["Hg", "Cd", "Ar"],
        # Old "long" headers (2017)
        "W_AIT_SRC_Ne": ["Ne"],
        "W_AIT_SRC_HgAr": ["Hg", "Ar"],
        "W_AIT_SRC_Xe": ["Xe"],
    }
    elements = [menu[key] for key in menu if md.get(key, False)]
    return {el for sublist in elements for el in sublist}
