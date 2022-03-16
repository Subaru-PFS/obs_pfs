import os
import re
from dataclasses import dataclass

__all__ = ("ParsedRawFilename", "parseRawFilename")


@dataclass
class ParsedRawFilename:
    """Results of parsing a raw filename

    Parameters
    ----------
    site : `str`
        Site at which raw file was created. May be one of:
        - ``J``: JHU
        - ``L``: LAM
        - ``S``: Subaru
        - ``F``: Fake, produced by the simulator
    category : `str`
        Category of file. May be one of:
        - ``A``: Raw data from the spectrographs.
        - ``B``: Up-the-ramp files from the NIR arms.
        - ``C``: Metrology camera.
        - ``D``: Auto-guider.
    exposure : `int`
        Exposure number. Also known as the "visit" number.
    arm : `str`
        The spectrograph arm. For spectrograph files, this is one of ``brnm``.
        For non-spectrograph files, this may be ``g`` for the guider.
    armNum : `int`
        Arm number. An index corresponding to the spectrograph arm.
    spectrograph : `int`
        Spectrograph module number. For spectrograph data, this is in the range
        1 to 4.
    sequence : `int`
        Sequence number for non-spectrograph data.
    """
    site: str
    category: str
    exposure: int
    arm: str
    armNum: int
    spectrograph: int
    sequence: int = None


def parseRawFilename(filename: str) -> ParsedRawFilename:
    """Parse the filename of a raw PFS file

    Parameters
    ----------
    filename : `str`
        Filename to be parsed.

    Returns
    -------
    result : `ParsedFilename`
        Parse results.

    Raises
    ------
    ValueError
        If the filename doesn't match the expected format.
    RuntimeError
        If the spectrograph number is out of range.
    """
    sites = '[JLXIASPF]'  # Possible site names
    categories = '[ABCDS]'  # Possible category names
    basename = os.path.basename(filename)
    matches = re.search(r"^PF(%s)(%s)-?(\d{6})([brnm0-9])([0-9])\.fits" % (sites, categories), basename)
    if not matches:
        raise ValueError("Unable to interpret filename: %s" % (filename,))
    site, category, exposure, first, second = matches.groups()
    exposure = int(exposure)
    if category == "A":  # spectrograph
        if first in "brnm":
            arm = first
            armNum = dict(b=1, r=2, n=3, m=4)[arm]
            spectrograph = int(second)
        else:
            spectrograph = int(first)
            armNum = int(second)
            arm = "brnm"[armNum - 1]
        if spectrograph < 1 or spectrograph > 4:
            raise RuntimeError(f"Unrecognised spectrograph ({spectrograph}) for base name {basename}")
        return ParsedRawFilename(site=site, category=category, exposure=exposure,
                                 spectrograph=spectrograph, arm=arm, armNum=armNum)
    elif category == "D":  # guider
        sequence = 10*int(first) + int(second)
        return ParsedRawFilename(site=site, category=category, exposure=exposure, spectrograph=0, arm="g",
                                 armNum=0, sequence=sequence)
