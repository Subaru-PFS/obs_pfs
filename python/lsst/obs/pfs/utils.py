"""
Utilities for working with PFS data
"""

import os


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
    # The lamp headers aren't perfectly reliable, so let's first check some others
    dataType = (md.get("DATA-TYP") or "").lower().strip()
    if dataType == "object":
        return set()
    if dataType == "flat":
        return set(["Quartz"])

    menu = {
        # The W_ENIISA headers were used for a short, but regrettable, time in early 2024
        "W_ENIISA": "Ar_eng",
        "W_ENIISC": "Quartz_eng",
        "W_ENIISH": "HgAr_eng",
        "W_ENIISK": "Kr_eng",
        "W_ENIISN": "Ne_eng",
        # THREE (count them, 3) characters for the elements
        "W_ENIARG": "Ar_eng",
        "W_ENIQTH": "Quartz_eng",
        "W_ENIHGA": "HgAr_eng",
        "W_ENIKRY": "Kr_eng",
        "W_ENINEO": "Ne_eng",
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

    lamps = getLamps(md)

    menu = {
        # Engineering fibres ("IIS")
        "Quartz_eng": [],
        "Ar_eng": ["Ar_eng"],
        "HgAr_eng": ["Hg_eng", "Ar_eng"],
        "Kr_eng": ["Kr_eng"],
        "Ne_eng": ["Ne_eng"],
        # Science fibres
        "Quartz": [],
        "Ar": ["Ar"],
        "HgAr": ["Hg", "Ar"],
        "HgCd": ["Hg", "Cd", "Ar"],
        "Kr": ["Kr"],
        "Ne": ["Ne"],
        "Xe": ["Xe"],
    }
    elements = [menu[key] for key in lamps]
    return {el for sublist in elements for el in sublist}


def getCalibPath(refOrButler):
    """Attempt to figure out the calibration path

    This only works with the Gen2 middleware, and involves digging in the
    internals.

    Parameters
    ----------
    refOrButler : `ButlerDataRef` or `Butler`
        Data reference or data butler.

    Returns
    -------
    path : `str`
        Path for calibs, or ``UNKNOWN``.
    """
    from lsst.daf.persistence import ButlerDataRef
    if isinstance(refOrButler, ButlerDataRef):
        butler = refOrButler.getButler()
    else:
        butler = refOrButler
    try:
        return os.path.realpath(butler._repos.inputs()[-1].cfg.mapperArgs["calibRoot"])
    except Exception:
        return "UNKNOWN"
