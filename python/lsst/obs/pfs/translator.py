import os

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, AltAz

from astro_metadata_translator import SubaruTranslator, cache_translation, PropertyDefinition
from astro_metadata_translator.translators.helpers import altaz_from_degree_headers

from lsst.utils import getPackageDir

from .utils import getLamps

__all__ = ["PfsTranslator"]


class PfsTranslator(SubaruTranslator):
    """Metadata translator for PFS FITS headers"""

    name = "PFS"
    """Name of this translation class"""

    supported_instrument = "PFS"
    """Supports the PFS instrument."""

    default_search_path = [os.path.join(getPackageDir("obs_pfs"), "corrections")]
    """Default search path to use to locate header correction files."""

    default_resource_root = os.path.join(getPackageDir("obs_pfs"), "corrections")
    """Default resource path root to use to locate header correction files."""

    _const_map = {"boresight_rotation_angle": Angle(0.0*u.deg),
                  "boresight_rotation_coord": "sky",
                  "telescope": "Subaru",
                  "observation_type": "science",
                  }
    """Constant mappings"""

    _trivial_map = {
        # Existing
        "exposure_id": "W_VISIT",
        "observation_id": "W_VISIT",
        "visit_id": "W_VISIT",
        "detector_group": "W_SPMOD",
        "exposure_time": ("EXPTIME", dict(unit=u.s)),
        "detector_serial": "DETECTOR",
        "object": "IMAGETYP",  # XXX to be updated; hopefully to OBJECT
        "ext_spectrograph": "W_SPMOD",
        "ext_pfs_design_id": "W_PFDSGN",
        "ext_dither": "W_ENFCAZ",
        "ext_shift": "W_ENFCAY",
        "ext_focus": "W_ENFCAX",
        # Guesses for what we will use in the future
        "boresight_airmass": ("AIRMASS", dict(default=np.nan)),
        "pressure": ("OUT-PRS", dict(unit=u.hPa, default=np.nan)),
        "relative_humidity": ("OUT-HUM", dict(default=np.nan)),
        "science_program": ("PROP-ID", dict(default=0)),
        "temperature": ("OUT-TMP", dict(unit=u.K, default=np.nan)),
    }
    """One-to-one mappings"""

    extensions = dict(
        pfs_design_id=PropertyDefinition("Top-end configuration ID", "int", int),
        spectrograph=PropertyDefinition("Spectrograph module [1-4]", "int", int),
        arm=PropertyDefinition("Arm of spectrograph [brnm]", "str", str),
        dither=PropertyDefinition("Slit offset in spatial dimension", "float", float),
        shift=PropertyDefinition("Slit offset in spectral dimension", "float", float),
        focus=PropertyDefinition("Focus offset", "float", float),
        lamps=PropertyDefinition("Lamps that are in use", "str", str),
    )

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
        return "INSTRUME" in header and header["INSTRUME"] == "PFS"

    @cache_translation
    def to_altaz_begin(self):
        if "ALTITUDE" not in self._header or "AZIMUTH" not in self._header:
            return AltAz(0.0*u.deg, 0.0*u.deg, obstime=self.to_datetime_begin(),
                         location=self.to_location())
        return altaz_from_degree_headers(self, (("ALTITUDE", "AZIMUTH"),),
                                         self.to_datetime_begin())

    @cache_translation
    def to_datetime_begin(self):
        # Neven writes:
        # > by comparing `DATE-OBS` for exposures with different `EXPTIME`, I am
        # > concluding that `DATE-OBS` refers to the end of the exposure
        value = self.to_datetime_end() - self._header["EXPTIME"]*u.second
        self._used_these_cards("EXPTIME")
        return value

    @cache_translation
    def to_datetime_end(self):
        value = self._from_fits_date_string(self._header["DATE-OBS"], scale="utc")
        self._used_these_cards("DATE-OBS")
        return value

    @cache_translation
    def to_detector_name(self):
        arm = self.to_ext_arm()
        spectrograph = self._header["W_SPMOD"]
        self._used_these_cards("W_SPMOD")
        return arm + str(spectrograph)

    @cache_translation
    def to_physical_filter(self):
        return self.to_ext_arm()

    @cache_translation
    def to_detector_exposure_id(self):
        return 100*self.to_exposure_id() + self.to_detector_num()

    @cache_translation
    def to_detector_num(self):
        arm = self._header["W_ARM"]
        spectrograph = self._header["W_SPMOD"]
        self._used_these_cards("W_ARM", "W_SPMOD")
        return 4*(spectrograph - 1) + arm - 1

    @cache_translation
    def to_tracking_radec(self):
        # Guess for what we'll use in the future
        ra = self._header["RA2000"] if "RA2000" in self._header else 0.0
        dec = self._header["DEC2000"] if "DEC2000" in self._header else 0.0
        radec = SkyCoord(ra, dec, frame="fk5", unit=(u.hourangle, u.deg),
                         obstime=self.to_datetime_begin(), location=self.to_location()).transform_to("icrs")
        self._used_these_cards("RA2000", "DEC2000")
        return radec

    @cache_translation
    def to_instrument(self):
        site = self._header["W_SITE"]
        self._used_these_cards("W_SITE")
        return "PFS" if site == "S" else f"PFS-{site}"

    @cache_translation
    def to_dark_time(self):
        if "DARKTIME" in self._header:
            darkTime = self._header["DARKTIME"]
            self._used_these_cards("DARKTIME")
            return darkTime*u.s
        return self.to_exposure_time()

    @cache_translation
    def to_ext_lamps(self):
        return ",".join(sorted(getLamps(self._header)))

    @cache_translation
    def to_ext_arm(self):
        armNum = self._header["W_ARM"]
        self._used_these_cards("W_ARM")
        return {1: "b", 2: "r", 3: "n", 4: "m"}[armNum]
