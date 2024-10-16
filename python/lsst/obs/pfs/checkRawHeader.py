from typing import Union

import astropy.io.fits
from astropy import log
from astropy.time import Time
from lsst.obs.pfs.parseFilename import parseRawFilename

__all__ = ("checkRawHeader",)


POD = Union[bool, str, int, float]


def checkKeyword(header: astropy.io.fits.Header, keyword: str, expected: POD, comment: str = None,
                 allowFix: bool = False) -> bool:
    """Check that the keyword exists in the header

    Parameters
    ----------
    header : `astropy.io.fits.Header`
        FITS header; may be updated if the ``keyword`` is not present.
    keyword : `str`
        Keyword of interest
    expected : `bool`, `str`, `int` or `float`
        Expected value of ``keyword``.
    comment : `str`
        Comment to give keyword in FITS header if updating.
    allowFix : `bool`, optional
        Allow the header to be updated with the expected value?

    Returns
    -------
    modified : `bool`
        Did we modify the header?

    Raises
    ------
    ValueError
        If the keyword is missing or has the incorrect value, and
        ``allowFix=False``.
    """
    if keyword in header:
        value = header[keyword]
        if value == expected:
            return False
        if allowFix:
            header[keyword] = (expected, comment)
            log.info(f"Fixed value of {keyword} = {repr(value)} --> {repr(expected)}")
            return True
        raise ValueError(f"Value mismatch for {keyword}: got {repr(value)} but expected {repr(expected)}")
    if allowFix:
        header[keyword] = (expected, comment)
        log.info(f"Added value of {keyword} = {repr(expected)}")
        return True
    raise ValueError(f"No header keyword {keyword}")


def checkRawHeader(filename: str, allowFix: bool = False):
    """Check that a PFS raw header includes the appropriate keywords

    These keywords include:
    - ``TELESCOP``: Telescope
    - ``INSTRUME``: Instrument
    - ``DETECTOR``: Detector name
    - ``DET-ID``: Detector number
    - ``W_VISIT``: PFS exposure visit number
    - ``W_ARM``: Spectrograph arm 1=b, 2=r, 3=n, 4=medRed
    - ``W_SPMOD``: Spectrograph module. 1-4 at Subaru
    - ``W_SITE``: PFS DAQ location: S=Subaru, J=JHU, L=LAM, A=ASIAA, F=simulator
    - ``MJD-STR``: Start time of exposure

    Parameters
    ----------
    filename : `str`
        Name of file to check.
    allowFix : `bool`, optional
        Allow fixing the header if it is non-conformant; by default ``False``.
    """
    data = parseRawFilename(filename)
    log.info(f"Checking {filename}")
    with astropy.io.fits.open(filename, "update" if allowFix else "readonly", save_backup=True) as fits:
        modified = False
        header = fits[0].header
        if header["NAXIS"] == 0:
            # Probably compressed
            modified |= checkKeyword(header, "INHERIT", True, "Inherit from previous header", allowFix)
            for key in ("TELESCOP", "INSTRUME", "DETECTOR", "W_VISIT", "W_ARM", "W_SPMOD", "W_SITE"):
                if key in header:
                    del header[key]
                    modified = True

            header = fits[1].header
            if header["NAXIS"] != 2:
                raise RuntimeError(f"Can't find main header for {filename}")

        try:
            modified |= checkKeyword(header, "TELESCOP", "Subaru 8.2m", "Telescope", allowFix)
            modified |= checkKeyword(header, "INSTRUME", "PFS", "Instrument", allowFix)
            modified |= checkKeyword(header, "W_SITE", data.site,
                                     "PFS DAQ site: S=Subaru, J=JHU L=LAM, F=sim", allowFix)
            modified |= checkKeyword(header, "W_VISIT", data.exposure, "PFS exposure visit number", allowFix)
            if data.category != "D":
                modified |= checkKeyword(header, "DETECTOR", f"{data.arm}{data.spectrograph}",
                                         "Name of the detector", allowFix)
                if data.arm == "m":
                    detId = -3*(data.spectrograph - 1) - 1
                else:
                    detId = 3*(data.spectrograph - 1) + data.armNum - 1
                modified |= checkKeyword(header, "DET-ID", detId, "Detector ID", allowFix)
                modified |= checkKeyword(header, "W_ARM", data.armNum,
                                         "Spectrograph arm: 1=b, 2=r, 3=n, 4=medRed", allowFix)
                modified |= checkKeyword(header, "W_SPMOD", data.spectrograph, "Spectrograph module: 1-4",
                                         allowFix)
            if "MJD-STR" not in header:
                mjd = Time(header["DATE-OBS"], scale="utc").mjd
                modified |= checkKeyword(header, "MJD-STR", mjd, "Start time of exposure", allowFix)
        except ValueError as exc:
            raise ValueError(f"Bad header for {filename}") from exc

        if modified:
            log.warning(f"Updated {filename}")
