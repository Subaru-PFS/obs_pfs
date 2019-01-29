import os
import re

import lsst.afw.image as afwImage
from lsst.obs.pfs.pfsMapper import PfsMapper
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig, IngestTask, IngestConfig
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask

from pfs.datamodel.pfsConfig import PfsConfig

__all__ = ["PfsParseConfig", "PfsParseTask", "PfsIngestTask"]


class PfsParseConfig(ParseConfig):
    """Configuration for PfsParseTask"""
    def setDefaults(self):
        ParseConfig.setDefaults(self)
        self.translators["field"] = "translate_field"
        self.translators["pfiDesignId"] = "translate_pfiDesignId"
        self.translators["slitOffset"] = "translate_slitOffset"


class PfsParseTask(ParseTask):
    """Parse a PFS FITS image to determine butler keyword-value pairs"""
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
        site, category, expId, spectrograph, armNum = matches.groups()

        armNum = int(armNum)
        spectrograph = int(spectrograph)
        expId = int(expId, base=10)

        # spectrograph 2 was called 9 at JHU, at least in the early days, so if we find
        # a spectrograph 9 we re-assign its number to 2
        if spectrograph == 9:
            spectrograph = 2
            self.log.info("Visit %06d has spectrograph == 9" % expId)

        if spectrograph < 1 or spectrograph > self.nSpectrograph + 1:
            raise RuntimeError('spectrograph (=%d) out of bounds 1..%d' % (spectrograph, self.nSpectrograph))

        try:
            arm = self.arms[armNum - 1]   # 1-indexed
        except IndexError as exc:
            raise IndexError('armNum=%d out of bounds 1..%d' % (armNum, max(self.arms.keys()))) from exc

        ccd = PfsMapper.computeDetectorId(spectrograph, arm)
        self.log.debug('site = <%s>, category = <%s>, expId = <%s>, spectrograph = <%d>, armNum = <%d>, '
                       'arm = <%s>, ccd = <%d>' %
                       (site, category, expId, spectrograph, armNum, arm, ccd))

        info = dict(site=site, category=category, expId=expId, visit=expId, filter=arm, arm=arm,
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

    def translate_pfiDesignId(self, md):
        """Get 'pfiDesignId' from metadata

        Fall back to the value to 0x0 if it's not present.
        """
        key = "W_PFDSGN"
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
    """Parse a PFS calib image for butler keyword-value pairs"""
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
        try:
            return int(self._translateFromCalibId("visit0", md))
        except Exception:
            return 0


class PfsIngestConfig(IngestConfig):
    """Configuration for PfsIngestTask"""
    def setDefaults(self):
        super().setDefaults()
        self.parse.retarget(PfsParseTask)
        self.register.columns = {'site': 'text',  # J: JHU, L: LAM, X: Subaru offline, I: IPMU, A: ASIAA,
                                                  # S: Summit, P: Princeton, F: simulation (fake)
                                 'category': 'text',  # A: science, B: NTR, C: Meterology, D: HG
                                 'field': 'text',  # Observation name
                                 'expId': 'int',  # Exposure identifier; better alias for "visit"
                                 'visit': 'int',  # Required because hard-coded in LSST's CameraMapper
                                 'ccd': 'int',  # [0-11]
                                 'filter': 'text',  # b: blue, r: red, n: nir, m: medium resolution red
                                 'arm': 'text',  # b: blue, r: red, n: nir, m: medium resolution red
                                 'spectrograph': 'int',  # Spectrograph module: 1-4
                                 'dateObs': 'text',  # Date of observation
                                 'expTime': 'double',  # Exposure time
                                 'dataType': 'text',  # Type of exposure
                                 'taiObs': 'text',  # Time of observation
                                 'pfiDesignId': 'int',  # Configuration of the top-end
                                 'slitOffset': 'double',  # Horizontal slit offset
                                 }
        self.register.unique = ['site', 'category', 'expId', 'visit', 'filter', 'arm', 'spectrograph',
                                'pfiDesignId']
        self.register.visit = ['expId', 'visit', 'field', 'filter', 'spectrograph', 'arm', 'dateObs',
                               'taiObs', 'pfiDesignId', 'slitOffset']

        self.parse.translation = {'dataType': 'IMAGETYP',
                                  'expTime': 'EXPTIME',
                                  'dateObs': 'DATE-OBS',
                                  'taiObs': 'DATE-OBS',
                                  }
        self.parse.defaults = {'ccdTemp': "0",  # Added in commissioning run 3
                               }
        self.parse.translators.update(field='translate_field',
                                      dateObs='translate_date',
                                      taiObs='translate_date')


class PfsIngestTask(IngestTask):
    """Ingest PFS images and configs into the data repository"""
    ConfigClass = PfsIngestConfig

    def ingestPfsConfig(self, dirName, fileInfo, args):
        """Ingest a PfsConfig file

        Parameters
        ----------
        dirName : `str`
            Directory in which the file resides.
        fileInfo : `dict`
            Key-value pairs defining the file.
        args : `argparse.Namespace`
            Parsed command-line arguments.
        """
        fileName = PfsConfig.fileNameFormat % (fileInfo["pfiDesignId"], fileInfo["expId"])
        infile = os.path.join(dirName, fileName)
        outfile = args.butler.get("pfsConfig_filename", fileInfo)[0]
        if not os.path.exists(outfile):
            self.ingest(infile, outfile, mode=args.mode, dryrun=args.dryrun)

    def runFile(self, infile, registry, args):
        """Examine and ingest a single PFS image file

        Also ingests the required PfsConfig file, which should be in the same
        directory as the raw image file.

        Parameters
        ----------
        infile : `str`
            Name of file to process
        registry : `lsst.pipe.tasks.ingest.RegistryContext` or `None`
            Registry into which to record file info.
        args : `argparse.Namespace`
            Parsed command-line arguments.

        Returns
        -------
        hduInfoList : `list` of `dict`
            Information from each of the HDUs.
        """
        try:
            hduInfoList = super().runFile(infile, registry, args)
            self.ingestPfsConfig(os.path.dirname(infile), hduInfoList[0], args)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            raise
        return hduInfoList
