import os
import re

import lsst.afw.image as afwImage
from lsst.obs.pfs.pfsMapper import PfsMapper
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask

__all__ = ["PfsParseConfig", "PfsParseTask"]


class PfsParseConfig(ParseConfig):
    def setDefaults(self):
        ParseConfig.setDefaults(self)
        self.translators["field"] = "translate_field"
        self.translators["pfsConfigId"] = "translate_pfsConfigId"
        self.translators["slitOffset"] = "translate_slitOffset"


class PfsParseTask(ParseTask):
    ConfigClass = PfsParseConfig
    nSpectrograph = 4  # Number of spectrographs
    arms = ['b', 'r', 'n', 'm']  # Possible spectrograph arm names
    sites = '[JLXIASPF]'  # Possible site names
    categories = '[ABCDS]'  # Possible category names

    def getInfo(self, filename):
        """Get information about the image from the filename and its contents

        The information we want isn't in the headers, but in the filname.

        Parameters
        ----------
        filename : `str`
            Name of file to inspect

        Returns
        -------
        info : `dict`
            File properties from the PHU.
        infoList : `list` of `dict`
            File properties from the extensions.
        """
        self.log.debug('interpreting filename <%s>' % filename)

        path, name = os.path.split(filename)

        matches = re.search(r"^PF(%s)(%s)-?(\d{6})(\d)(\d)\.fits$" % (self.sites, self.categories), name)
        if not matches:
            raise RuntimeError("Unable to interpret filename: %s" % filename)
        site, category, visit, spectrograph, armNum = matches.groups()

        armNum = int(armNum)
        spectrograph = int(spectrograph)
        visit = int(visit, base=10)

        # spectrograph 2 was called 9 at JHU, at least in the early days, so if we find
        # a spectrograph 9 we re-assign its number to 2
        if spectrograph == 9:
            spectrograph = 2
            self.log.info("Visit %06d has spectrograph == 9" % visit)

        if spectrograph < 1 or spectrograph > self.nSpectrograph + 1:
            raise RuntimeError('spectrograph (=%d) out of bounds 1..%d' % (spectrograph, self.nSpectrograph))

        try:
            arm = self.arms[armNum - 1]   # 1-indexed
        except IndexError as exc:
            raise IndexError('armNum=%d out of bounds 1..%d' % (armNum, max(self.arms.keys()))) from exc

        ccd = PfsMapper.computeDetectorId(spectrograph, arm)
        self.log.debug('site = <%s>, category = <%s>, visit = <%s>, spectrograph = <%d>, armNum = <%d>, '
                       'arm = <%s>, ccd = <%d>' %
                       (site, category, visit, spectrograph, armNum, arm, ccd))

        info = dict(site=site, category=category, visit=visit, filter=arm, arm=arm,
                    spectrograph=spectrograph, ccd=ccd)

        if os.path.exists(filename):
            header = afwImage.readMetadata(filename)
            info = self.getInfoFromMetadata(header, info=info)
        return info, [info]

    def translate_field(self, md):
        """Get 'field' from IMAGETYP

        This is temporary, until an OBJECT (or similar) header can be provided,
        but it's better than setting everything to the same thing.
        """
        field = md.get("IMAGETYP").strip()
        if field in ("#", ""):
            field = "UNKNOWN"
        self.log.debug('PfsParseTask.translate_field: field = %s' % field)
        return re.sub(r'\W', '_', field).upper()

    def translate_pfsConfigId(self, md):
        """Get 'pfsConfigId' from metadata

        Fall back to the value to 0x0 if it's not present.
        """
        key = "PFSCONFIGID"
        if md.exists(key):
            return md.get(key)
        return 0x0

    def translate_slitOffset(self, md):
        """Get slitOffset from header metadata

        We try mutliple header keywords before resorting to an assumed zero
        offset.
        """
        for key in ("sim.slit.xoffset",  # From the simulator
                    "W_FCA_DITHER"  # From real data at LAM
                    ):
            if md.exists(key):
                value = md.get(key)
                return value if isinstance(value, float) else 0.0
        return 0.0


class PfsCalibsParseTask(CalibsParseTask):
    ConfigClass = PfsParseConfig
    calibTypes = ["Bias", "Dark", "DetectorMap", "FiberTrace", "FiberFlat"]

    def _translateFromCalibId(self, field, md):
        data = md.get("CALIB_ID")
        match = re.search(".*%s=(\S+)" % field, data)
        return match.groups()[0]

    def translate_ccd(self, md):
        return int(self._translateFromCalibId("ccd", md))

    def translate_arm(self, md):
        return self._translateFromCalibId("arm", md)

    def translate_calibDate(self, md):
        return self._translateFromCalibId("calibDate", md)

    def translate_spectrograph(self, md):
        return int(self._translateFromCalibId("spectrograph", md))

    def translate_visit0(self, md):
        # We don't have a way of setting this yet; it's always zero
        return 0
