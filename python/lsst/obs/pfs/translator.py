import os
import re
from types import SimpleNamespace

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

from astro_metadata_translator import SubaruTranslator, cache_translation
from astro_metadata_translator.translators.helpers import altaz_from_degree_headers

from lsst.utils import getPackageDir

__all__ = ["PfsTranslator"]


class PfsTranslator(SubaruTranslator):
    """Metadata translator for PFS FITS headers"""

    name = "PFS"
    """Name of this translation class"""

    supported_instrument = "PFS"
    """Supports the PFS instrument."""

    default_search_path = os.path.join(getPackageDir("obs_pfs"), "corrections")
    """Default search path to use to locate header correction files."""

    default_resource_root = os.path.join(getPackageDir("obs_pfs"), "corrections")
    """Default resource path root to use to locate header correction files."""

    _const_map = {"boresight_rotation_angle": "unknown",
                  "boresight_rotation_coord": "unknown",
                  "telescope": "Subaru",
                  }
    """Constant mappings"""

    _trivial_map = {
        # Existing
        "dark_time": ("DARKTIME", dict(unit=u.s)),
        "detector_group": "SPECNUM",  # Spectrograph
        "detector_name": "ARM",
        "exposure_time": ("EXPTIME", dict(unit=u.s)),
        "object": "IMAGETYP",  # XXX to be updated; hopefully to OBJECT
        "observation_type": "IMAGETYP",
        "physical_filter": "ARM",
        # Guesses for what we will use in the future
        "boresight_airmass": ("AIRMASS", dict(default=np.nan)),
        "pressure": ("OUT-PRS", dict(unit=u.hPa, default=np.nan)),
        "relative_humidity": ("OUT-HUM", dict(default=np.nan)),
        "science_program": ("PROP-ID", dict(default=0)),
        "temperature": ("OUT-TMP", dict(unit=u.K, default=np.nan)),
    }
    """One-to-one mappings"""

    _armToIndex = {"blue": 1,
                   "red": 2,
                   "nir": 3,
                   "medium": 4,
                   }

    def __init__(self, header, filename=None):
        super().__init__(header, filename=filename)
        if filename is not None:
            self._parsedData = self.parseFilename()

    @classmethod
    def can_translate(cls, header, filename=None):
        """Indicate whether this translation class can translate the
        supplied header.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        if "INSTRUME" in header:  # Not in header as of August 2019
            return header["INSTRUME"] == "PFS"
        if "HIERARCH revision.FEE" in header and header["HIERARCH revision.FEE"].startswith("PFS"):
            return True
        return False

    def parseFilename(self):
        if self.filename is None:
            raise ValueError("No filename to parse")
        sites = '[JLXIASPF]'  # Possible site names
        categories = '[ABCDS]'  # Possible category names
        matches = re.search(r"^PF(%s)(%s)-?(\d{6})(\d)(\d)\.fits$" % (sites, categories),
                            os.path.basename(self.filename))
        if not matches:
            raise ValueError("Unable to interpret filename: %s" % (self.filename,))
        site, category, visit, spectrograph, arm = matches.groups()
        return SimpleNamespace(site=site, category=category, visit=int(visit),
                               spectrograph=int(spectrograph), arm=int(arm))

    @cache_translation
    def to_altaz_begin(self):
        return altaz_from_degree_headers(self, (("ALTITUDE", "AZIMUTH"),),
                                         self.to_datetime_begin())

    @cache_translation
    def to_datetime_begin(self):
        # Neven writes:
        # > by comparing `DATE-OBS` for exposures with different `EXPTIME`, I am
        # > concluding that `DATE-OBS` refers to the end of the exposure
        value = self.to_datetime_begin() - self._header["EXPTIME"]
        self._used_these_cards("EXPTIME")
        return value

    @cache_translation
    def to_datetime_end(self):
        value = self._from_fits_date_string(self._header["DATE-OBS"], scale="utc")
        self._used_these_cards("DATE-OBS")
        return value

    @cache_translation
    def to_detector_exposure_id(self):
        return 100*self.to_exposure_id() + self.to_detector_num()

    @cache_translation
    def to_detector_num(self):
        return 10*self.to_detector_group() + self._armToIndex(self.to_detector_name())

    @cache_translation
    def to_detector_serial(self):
        value = self._header["HIERARCH serial.CCD0"] + "+" + self._header["HIERARCH serial.CCD1"]
        self._used_these_cards("HIERARCH serial.CCD0", "HIERARCH serial.CCD1")
        return value

    @cache_translation
    def to_exposure_id(self):
        return self.to_observation_id()

    @cache_translation
    def to_observation_id(self):
        if "W_VISIT" in self._header:
            self._used_these_cards("W_VISIT")
            return self._header["W_VISIT"]
        if self.filename is not None:
            return self._parsedData.visit
        raise ValueError("Can't determine visit number")

    @cache_translation
    def to_tracking_radec(self):
        # Guess for what we'll use in the future
        ra = self._header["RA2000"] if "RA2000" in self._header else np.nan
        dec = self._header["DEC2000"] if "DEC2000" in self._header else np.nan
        radec = SkyCoord(ra, dec, frame="fk5", unit=(u.hourangle, u.deg),
                         obstime=self.to_datetime_begin(), location=self.to_location()).transform_to("icrs")
        self._used_these_cards("RA2000", "DEC2000")
        return radec

    @cache_translation
    def to_visit_id(self):
        return self.to_observation_id()

    @cache_translation
    def to_instrument(self):
        if self._parsedData is None:
            self._parsedData = self.parseFilename(self.filename)
        return "PFS" if self._parsedData.site == "S" else f"PFS-{self._parsedData.site}"

    def search_paths(self):
        """Search paths to use when searching for header fix up correction
        files.

        This method is necessary for LSST 18.1.0. It can be removed once we
        upgrade past ``w.2019.25``.

        Returns
        -------
        paths : `list`
            Directory paths to search. Can be an empty list if no special
            directories are defined.

        Notes
        -----
        Uses the classes ``default_search_path`` property if defined.
        """
        return [self.default_search_path]
