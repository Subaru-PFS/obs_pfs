import os
import re

import lsst.afw.image as afwImage
from lsst.obs.pfs.pfsMapper import PfsMapper
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig, IngestTask, IngestConfig, IngestArgumentParser
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask
from lsst.pex.config import Field

from pfs.datamodel.pfsConfig import PfsConfig, PfsDesign

__all__ = ["PfsParseConfig", "PfsParseTask", "PfsIngestTask"]


class PfsParseConfig(ParseConfig):
    """Configuration for PfsParseTask"""
    pfsDesignId = Field(dtype=int, default=0x0, doc="Default value for pfsDesignId")

    def setDefaults(self):
        ParseConfig.setDefaults(self)
        self.translators["field"] = "translate_field"
        self.translators["pfsDesignId"] = "translate_pfsDesignId"
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

    def translate_pfsDesignId(self, md):
        """Get 'pfsDesignId' from metadata

        Fall back to the value to 0x0 if it's not present.
        """
        key = "W_PFDSGN"
        if md.exists(key):
            value = md.get(key)
            if isinstance(value, int):
                return value
        self.log.warn("No value set for pfsDesignId; using default (0x%016x)" % (self.config.pfsDesignId,))
        return self.config.pfsDesignId

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


def setIngestConfig(config):
    """Set the configuration for ingestion

    This has been factored out so it can be used from
    ``PfsIngestConfig.setDefaults`` (run via ``ingestPfsImages.py``) or from
    the obs package override for ``ingestImages.py``.

    Parameters
    ----------
    config : subclass of `lsst.pipe.tasks.IngestConfig`
        Configuration to set.
    """
    config.parse.retarget(PfsParseTask)
    config.register.columns = {'site': 'text',  # J: JHU, L: LAM, X: Subaru offline, I: IPMU, A: ASIAA,
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
                               'pfsDesignId': 'int',  # Configuration of the top-end
                               'slitOffset': 'double',  # Horizontal slit offset
                               }
    config.register.unique = ['site', 'category', 'expId', 'visit', 'filter', 'arm', 'spectrograph',
                              'pfsDesignId']
    config.register.visit = ['expId', 'visit', 'field', 'filter', 'spectrograph', 'arm', 'dateObs',
                             'taiObs', 'pfsDesignId', 'slitOffset']

    config.parse.translation = {'dataType': 'IMAGETYP',
                                'expTime': 'EXPTIME',
                                'dateObs': 'DATE-OBS',
                                'taiObs': 'DATE-OBS',
                                }
    config.parse.defaults = {'ccdTemp': "0",  # Added in commissioning run 3
                             }
    config.parse.translators.update(field='translate_field',
                                    dateObs='translate_date',
                                    taiObs='translate_date')


class PfsIngestArgumentParser(IngestArgumentParser):
    """ArgumentParser for PFS ingestion

    Sets PFS-specific command-line arguments.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_argument("--pfsConfigDir",
                          help="Directory with pfsConfig/pfsDesign files (default: with images)")


class PfsIngestConfig(IngestConfig):
    """Configuration for PfsIngestTask"""
    def setDefaults(self):
        super().setDefaults()
        setIngestConfig(self)


class PfsIngestTask(IngestTask):
    """Ingest PFS images and configs into the data repository"""
    ConfigClass = PfsIngestConfig
    ArgumentParser = PfsIngestArgumentParser
    _DefaultName = "ingestPfs"

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
        outfile = args.butler.get("pfsConfig_filename", fileInfo)[0]
        if os.path.exists(outfile):
            # Don't clobber one that got put there when we ingested a different spectrograph,arm
            return

        pfsDesignId = fileInfo["pfsDesignId"]
        expId = fileInfo["expId"]
        fileName = PfsConfig.fileNameFormat % (pfsDesignId, expId)
        infile = os.path.join(dirName, fileName)
        if os.path.exists(infile):
            self.ingest(infile, outfile, mode=args.mode, dryrun=args.dryrun)
            return

        # Check for a PfsDesign, and use that instead
        designName = os.path.join(dirName, PfsDesign.fileNameFormat % (pfsDesignId,))
        if not os.path.exists(designName):
            raise RuntimeError("Unable to find PfsConfig or PfsDesign for pfsDesignId=0x%016x" %
                               (pfsDesignId,))
        design = PfsDesign.read(pfsDesignId, dirName)
        keywords = ("pfsDesignId", "raBoresight", "decBoresight",
                    "fiberId", "tract", "patch", "ra", "dec", "catId", "objId",
                    "targetType", "fiberMag", "filterNames", "pfiNominal")
        kwargs = {kk: getattr(design, kk) for kk in keywords}
        kwargs["expId"] = expId
        kwargs["pfiCenter"] = kwargs["pfiNominal"]
        PfsConfig(**kwargs).write(dirName)
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
            pfsConfigDir = args.pfsConfigDir if args.pfsConfigDir is not None else os.path.dirname(infile)
            self.ingestPfsConfig(pfsConfigDir, hduInfoList[0], args)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            raise
        return hduInfoList
