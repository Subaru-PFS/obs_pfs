import os
import pathlib

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time

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

    # Provided by superclasses:
    # * location (SubaruTranslator)
    # * observation_counter (SubaruTranslator)
    # * observing_day (MetadataTranslator)

    _const_map = {
        "boresight_rotation_coord": "mount",
        "group_counter_end": 0,
        "group_counter_start": 0,
        "telescope": "Subaru",
    }
    """Constant mappings"""

    _trivial_map = {
        # Existing
        "boresight_airmass": ("AIRMASS", dict(default=np.nan)),
        "boresight_rotation_angle": ("INSROT", dict(unit=u.deg, default=0.0)),
        "detector_group": "W_SPMOD",
        "detector_serial": "DETECTOR",
        "detector_unique_name": "DETECTOR",
        "exposure_group": ("W_PFDSNM", dict(default="")),
        "exposure_time": ("EXPTIME", dict(unit=u.s)),
        "focus_z": ("FOC-VAL", dict(unit=u.um, default=np.nan)),
        "pressure": ("OUT-PRS", dict(unit=u.hPa, default=np.nan)),
        "relative_humidity": ("OUT-HUM", dict(default=np.nan)),
        "science_program": ("PROP-ID", dict(default="unknown")),
        "temperature": ("OUT-TMP", dict(unit=u.K, default=np.nan)),
        "ext_spectrograph": "W_SPMOD",
        "ext_pfs_design_id": "W_PFDSGN",
        "ext_dither": "W_ENFCAZ",
        "ext_shift": "W_ENFCAY",
        "ext_focus": "W_ENFCAX",
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
        return (("INSTRUME" in header and header["INSTRUME"] == "PFS") or
                ("SLIT" in header and header["SLIT"] == "PFS"))

    def from_header_string(self, *keywords):
        """Return the first non-empty value from the header

        Parameters
        ----------
        *keywords : `str`
            Keywords to search for in the header

        Returns
        -------
        value : `str`
            First non-empty value found in the header, or `None` if none found.
        """
        for key in keywords:
            if key not in self._header:
                continue
            value = self._header[key].strip()
            if not value:
                continue
            self._used_these_cards(key)
            return self._header[key]
        return None

    @cache_translation
    def to_visit_id(self):
        if self._header["W_SITE"] == "J":
            path = pathlib.Path(self.filename)
            visit_id = int(path.stem[4:-2], base=10)
        else:
            visit_id = self._header["W_VISIT"]
        return visit_id

    @cache_translation
    def to_exposure_id(self):
        return self.to_visit_id()

    @cache_translation
    def to_observation_id(self):
        if self._header["W_SITE"] == "J":
            path = pathlib.Path(self.filename)
            visit_id = path.stem[4:-2]
        else:
            visit_id = str(self._header["W_VISIT"])
        return visit_id



    @cache_translation
    def to_object(self):
        # IMAGETYP is to support the simulator, which doesn't set OBJECT
        value = self.from_header_string("OBJECT", "IMAGETYP")
        return value if value else "UNKNOWN"

    @cache_translation
    def to_altaz_begin(self):
        if "ALTITUDE" not in self._header or "AZIMUTH" not in self._header or self._header["ALTITUDE"] == 'no available value':
            return AltAz(90.0*u.deg, 0.0*u.deg, obstime=self.to_datetime_begin(),
                         location=self.to_location())
        return altaz_from_degree_headers(self, (("ALTITUDE", "AZIMUTH"),),
                                         self.to_datetime_begin())

    @cache_translation
    def to_datetime_begin(self):
        value = Time(self._header["MJD-STR"], format="mjd", scale=self._header.get("TIMESYS", "utc").lower())
        self._used_these_cards("MJD-STR", "TIMESYS")
        return value

    @cache_translation
    def to_datetime_end(self):
        return self.to_datetime_begin() + self.to_exposure_time()

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
        armNum = self._header["W_ARM"]
        spectrograph = self._header["W_SPMOD"]
        self._used_these_cards("W_ARM", "W_SPMOD")
        if armNum == 4:  # arm=m
            return - 3*(spectrograph - 1) - 1
        return 3*(spectrograph - 1) + armNum - 1

    @cache_translation
    def to_tracking_radec(self):
        # Guess for what we'll use in the future
        ra = self._header["RA2000"] if "RA2000" in self._header else 0.0
        dec = self._header["DEC2000"] if "DEC2000" in self._header else 0.0
        if ra == 'no available value':
            ra = 0.0
        if dec == 'no available value':
            dec = 0.0
        radec = SkyCoord(ra, dec, frame="fk5", unit=(u.hourangle, u.deg),
                         obstime=self.to_datetime_begin(), location=self.to_location()).transform_to("icrs")
        self._used_these_cards("RA2000", "DEC2000")
        return radec

    @cache_translation
    def to_instrument(self):
        site = self._header["W_SITE"]
        self._used_these_cards("W_SITE")
        return "PFS" if site in {"S", "J"} else f"PFS-{site}"

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

    @cache_translation
    def to_can_see_sky(self):
        domeShutter = self._header.get("W_TSHUTR", "closed")
        self._used_these_cards("W_TSHUTR")
        return domeShutter == "open"

    @cache_translation
    def to_has_simulated_content(self):
        site = self._header["W_SITE"]
        self._used_these_cards("W_SITE")
        return site == "F"

    @cache_translation
    def to_observation_reason(self):
        obsType = self.to_observation_type()
        return dict(science="science", UNKNOWN="UNKNOWN").get(obsType, "calibration")

    @cache_translation
    def to_observation_type(self):
        value = self.from_header_string("DATA-TYP", "IMAGETYP")
        if value is None:
            return "UNKNOWN"
        if value == "OBJECT":
            return "science"
        return value.lower()
