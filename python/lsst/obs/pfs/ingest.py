import os
import re
import shutil
from functools import partialmethod
import datetime
import dateutil.parser
import collections

from astro_metadata_translator import fix_header

from lsst.obs.pfs.pfsMapper import PfsMapper
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig, IngestTask, IngestConfig, IngestArgumentParser
from lsst.pipe.tasks.ingest import RegisterTask, IngestError
from lsst.pipe.tasks.ingestCalibs import IngestCalibsTask, IngestCalibsConfig
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask, CalibsRegisterTask
from lsst.pipe.tasks.ingestPgsql import PgsqlRegisterTask
from lsst.pex.config import Field, ConfigurableField
from lsst.afw.fits import readMetadata
from lsst.obs.pfs.utils import getLamps
import lsst.pipe.tasks.ingest

from pfs.datamodel.pfsConfig import PfsConfig, PfsDesign
from .translator import PfsTranslator

__all__ = ["PfsParseConfig", "PfsParseTask", "PfsIngestTask", "PfsIngestCalibsTask", "PfsPgsqlIngestTask"]


class PfsParseConfig(ParseConfig):
    """Configuration for PfsParseTask"""
    pfsDesignId = Field(dtype=int, default=0x0, doc="Default value for pfsDesignId")

    def setDefaults(self):
        ParseConfig.setDefaults(self)
        self.translators["field"] = "translate_field"
        self.translators["dataType"] = "translate_dataType"
        self.translators["pfsDesignId"] = "translate_pfsDesignId"
        self.translators["slitOffset"] = "translate_slitOffset"
        self.translators["lamps"] = "translate_lamps"
        self.translators["dateObs"] = "translate_date"
        self.translators["dither"] = "translate_dither"
        self.translators["shift"] = "translate_shift"
        self.translators["focus"] = "translate_focus"
        self.translators["attenuator"] = "translate_attenuator"
        self.translators["photodiode"] = "translate_photodiode"
        self.translators["taiObs"] = "translate_datetime"


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
        site, category, visit, spectrograph, armNum = matches.groups()

        info = dict(site=site, category=category, visit=int(visit, base=10))

        if category in ("A", "B"):
            # Spectrograph
            armNum = int(armNum)
            spectrograph = int(spectrograph)

            # spectrograph 2 was called 9 at JHU, at least in the early days, so if we find
            # a spectrograph 9 we re-assign its number to 2
            if spectrograph == 9:
                spectrograph = 2
                self.log.info("Visit %06d has spectrograph == 9" % visit)

            if spectrograph < 1 or spectrograph > self.nSpectrograph + 1:
                raise RuntimeError('spectrograph (=%d) out of bounds 1..%d' %
                                   (spectrograph, self.nSpectrograph))

            try:
                arm = self.arms[armNum - 1]   # 1-indexed
            except IndexError as exc:
                raise IndexError('armNum=%d out of bounds 1..%d' % (armNum, max(self.arms.keys()))) from exc

            ccd = PfsMapper.computeDetectorId(spectrograph, arm)
            self.log.debug('site = <%s>, category = <%s>, visit = <%s>, spectrograph = <%d>, armNum = <%d>, '
                           'arm = <%s>, ccd = <%d>' %
                           (site, category, visit, spectrograph, armNum, arm, ccd))

            info["filter"] = arm
            info["arm"] = arm
            info["spectrograph"] = spectrograph
            info["ccd"] = ccd
        elif category == "D":
            # Guider
            info["sequence"] = 10*int(spectrograph) + int(armNum)
            # Necessary dummy values
            info["filter"] = "NONE"
            info["arm"] = "x"
            info["ccd"] = 0
            info["spectrograph"] = 0
            info["hdu"] = 0
        else:
            raise RuntimeError(f"Unrecognised file category ({category}) for file {filename}")

        if os.path.exists(filename):
            header = readMetadata(filename)

            # Patch the header
            original = header.toOrderedDict()
            fixed = header.toOrderedDict()
            fix_header(fixed, translator_class=PfsTranslator, filename=filename)
            for kk, vv in fixed.items():
                if kk in original and vv == original[kk]:
                    continue
                header.set(kk, vv)

            info = self.getInfoFromMetadata(header, info=info)
        return info, [info]

    def getInfoFromMetadata(self, md, info=None):
        """Attempt to pull the desired information out of the header

        This is done through two mechanisms:
        * translation: a property is set directly from the relevant header keyword
        * translator: a property is set with the result of calling a method

        The translator methods receive the header metadata and should return the
        appropriate value, or None if the value cannot be determined.

        This method has been copied from LSST's pipe_tasks in order to fix a
        bug in LSST 18.1.0 (PIPE2D-458).

        @param md      FITS header
        @param info    File properties, to be supplemented
        @return info
        """
        if info is None:
            info = {}
        for p, h in self.config.translation.items():
            value = md.get(h, None)
            if value is not None:
                if isinstance(value, str):
                    value = value.strip()
                info[p] = value
            elif p in self.config.defaults:
                info[p] = self.config.defaults[p]
            else:
                self.log.warn("Unable to find value for %s (derived from %s)" % (p, h))
        for p, t in self.config.translators.items():
            func = getattr(self, t)
            try:
                value = func(md)
            except Exception as e:
                self.log.warn("%s failed to translate %s: %s", t, p, e)
                value = None
            if value is not None:
                info[p] = value
        return info

    def getDestination(self, butler, info, filename):
        """Get destination for the file

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
            Data butler.
        info : `dict`
            Data identifier keyword-value pairs.
        filename : `str`
            Input filename.

        Returns
        -------
        destination : `str`
            Destination filename.
        """
        category = info["category"]
        dataset = {"A": "raw",
                   "B": "raw",
                   "D": "guider",
                   }[category]
        raw = butler.get(dataset + "_filename", info)[0]
        # Ensure filename is devoid of cfitsio directions about HDUs
        c = raw.find("[")
        if c > 0:
            raw = raw[:c]
        return raw

    def translate_field(self, md):
        """Get 'field' from metadata

        Get a potentially useful value from typical header keywords. As PFS
        evolves, the header keyword with the value we want changes (IMAGETYP
        for engineering at LAM, DATA-TYP for engineering at Subaru, and then
        OBJECT when we get on the sky).
        """
        for key in ("OBJECT", "DATA-TYP", "IMAGETYP"):
            if not md.exists(key):
                continue
            field = md.get(key).strip()
            if field not in ("#", "##NODATA##", ""):
                break
        else:
            field = "UNKNOWN"
        self.log.debug('PfsParseTask.translate_field: field = %s' % field)
        return re.sub(r'\W', '_', field).upper()

    def translate_dataType(self, md):
        """Get 'dataType' from metadata

        Get a value from typical header keywords. As PFS evolves, the header
        keyword with the value we want changes (IMAGETYP for engineering at LAM,
        DATA-TYP for engineering and on the sky at Subaru.
        """
        for key in ("DATA-TYP", "IMAGETYP"):
            if not md.exists(key):
                continue
            field = md.get(key).strip()
            if field not in ("#", "##NODATA##", ""):
                break
        else:
            field = "UNKNOWN"
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

    def translate_lamps(self, md):
        """Get lamps from header metadata"""
        return ",".join(sorted(getLamps(md.toDict())))

    def getNumeric(self, md, keyword, default=0.0):
        """Get a numerical value from a header which may be a string

        Casting these to numeric types results in an exception, so we set them
        to a default value.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList`
            The FITS header.
        keyword : `str`
            Header keyword from which to get the value.
        default : numeric
            Value to use if the header keyword isn't present, or the value is
            not a numeric type.

        Returns
        -------
        value : numeric
            Value of interest.
        """
        if not md.exists(keyword):
            return default
        value = md.get(keyword)
        try:
            return float(value)
        except Exception:
            return default

    def translate_datetime(self, md):
        """Get an ISO-formatted date+time

        Sometimes the time is included in DATE-OBS, sometimes it isn't.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList`
            Metadata from the header.

        Returns
        -------
        datetime : `str`
            ISO-formatted date+time.
        """
        if "DATE-OBS" not in md and "DATE" in md:
            # Early guider data
            return md.get("DATE")
        dateObs = md.get("DATE-OBS")
        if "T" in dateObs:
            return dateObs
        time = md.get("UT")
        return f"{dateObs}T{time}"

    def translate_date(self, md):
        """Convert a full DATE-OBS to a mere date

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList`
            Metadata from the header.

        Returns
        -------
        datetime : `str`
            ISO-formatted date+time.
        """
        if "DATE-OBS" in md:
            date = md.getScalar("DATE-OBS").strip()
        elif "DATE" in md:
            date = md.getScalar("DATE").strip()
        else:
            raise RuntimeError("Unable to find DATE-OBS or DATE keyword in header")
        c = date.find("T")
        if c > 0:
            date = date[:c]
        return date


# Set up multiple columns to be populated using the "getNumeric" method
for column, keyword, default in (('dither', 'W_ENFCAZ', -9999),
                                 ('shift', 'W_ENFCAY', -9999),
                                 ('focus', 'W_ENFCAX', -9999),
                                 ('attenuator', 'W_AITATT', -9999),
                                 ('photodiode', 'W_AITPHO', -9999),
                                 ):
    setattr(PfsParseTask, "translate_" + column,
            partialmethod(PfsParseTask.getNumeric, keyword=keyword, default=default))


class PfsCalibsParseTask(CalibsParseTask):
    """Parse a PFS calib image for butler keyword-value pairs"""
    ConfigClass = PfsParseConfig

    def getCalibType(self, filename):
        """Return the calibration dataset type

        We inspect the OBSTYPE header.

        Parameters
        ----------
        filename : `str`
            Input filename.

        Returns
        -------
        calibType : `str`
            Calibration dataset type.
        """
        md = readMetadata(filename, self.config.hdu)
        if not md.exists("OBSTYPE"):
            raise RuntimeError("Unable to find the required header keyword OBSTYPE in %s, hdu %d" %
                               (filename, self.config.hdu))
        return md.getScalar("OBSTYPE").strip()

    def _translateFromCalibId(self, field, md):
        data = md.get("CALIB_ID")
        match = re.search(r".*%s=(\S+)" % field, data)
        return match.groups()[0]

    def translate_ccd(self, md):
        return int(self._translateFromCalibId("ccd", md))

    def translate_arm(self, md):
        return self._translateFromCalibId("arm", md)

    def translate_calibDate(self, md):
        return self._translateFromCalibId("calibDate", md)

    def translate_calibTime(self, md):
        """Get time of calibration frame from metadata

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList`
            Metadata from header.

        Returns
        -------
        calibTime : `str`
            Time of calibration frame.
        """
        try:
            return self._translateFromCalibId("calibTime", md)
        except Exception:
            calibDate = self._translateFromCalibId("calibDate", md)
            return calibDate[:calibDate.find("T")] if "T" in calibDate else calibDate

    def translate_spectrograph(self, md):
        return int(self._translateFromCalibId("spectrograph", md))

    def translate_visit0(self, md):
        try:
            return int(self._translateFromCalibId("visit0", md))
        except Exception:
            return 0


class PfsRegisterTask(RegisterTask):
    def addRow(self, conn, info, dryrun=False, create=False, table=None):
        """Add a row to the file table (typically 'raw').

        This override of the LSST implementation removes the visit table
        because PFS doesn't use the visits table, and its creation can be
        expensive. Furthermore, the LSST implementation hard-codes the
        "OR IGNORE" part of the SQL, which doesn't work for PostgreSQL.

        Parameters
        ----------
        conn : `sqlite3.Connection`
            Database connection.
        info : `dict`
            File properties to add to database.
        dryrun : `bool`, optional
            Simulate what would happen?
        create : `bool`, optional
            Create a table if it doesn't exist? This parameter appears to be
            vestigial.
        table : `str`, optional
            Name of table.
        """
        with conn:
            if table is None:
                table = self.config.table
            ignoreClause = ""
            if self.config.ignore:
                ignoreClause = " OR IGNORE"
            sql = "INSERT%s INTO %s (%s) VALUES (" % (ignoreClause, table, ",".join(self.config.columns))
            sql += ",".join([self.placeHolder] * len(self.config.columns)) + ")"
            values = [self.typemap[tt](info[col]) for col, tt in self.config.columns.items()]

            if dryrun:
                print("Would execute: '%s' with %s" % (sql, ",".join([str(value) for value in values])))
            else:
                conn.cursor().execute(sql, values)


def setIngestConfig(config, retarget=True):
    """Set the configuration for ingestion

    This has been factored out so it can be used from
    ``PfsIngestConfig.setDefaults`` (run via ``ingestPfsImages.py``) or from
    the obs package override for ``ingestImages.py``.

    Parameters
    ----------
    config : subclass of `lsst.pipe.tasks.IngestConfig`
        Configuration to set.
    retarget : `bool`, optional
        Perform the retargeting? This should be disabled if the ``config`` is
        already set up to use the appropriate ``parse`` and ``register`` tasks
        (e.g., because it was subclassed over overridden).
    """
    if retarget:
        config.parse.retarget(PfsParseTask)
        config.register.retarget(PfsRegisterTask)
    config.register.columns = {'site': 'text',  # J: JHU, L: LAM, X: Subaru offline, I: IPMU, A: ASIAA,
                                                # S: Summit, P: Princeton, F: simulation (fake)
                               'category': 'text',  # A: science, B: NTR, C: Meterology, D: HG
                               'field': 'text',  # Observation name
                               'visit': 'int',  # Required because hard-coded in LSST's CameraMapper
                               'ccd': 'int',  # [0-11]
                               'filter': 'text',  # b: blue, r: red, n: nir, m: medium resolution red
                               'arm': 'text',  # b: blue, r: red, n: nir, m: medium resolution red
                               'spectrograph': 'int',  # Spectrograph module: 1-4
                               'dateObs': 'text',  # Date of observation; used for filenames
                               'expTime': 'double',  # Exposure time
                               'dataType': 'text',  # Type of exposure
                               'taiObs': 'text',  # Date+time of observation; used for finding calibs
                               'pfsDesignId': 'int',  # Configuration of the top-end
                               'slitOffset': 'double',  # Horizontal slit offset; kept for backwards compat.
                               'dither': 'double',  # Slit offset in spatial dimension
                               'shift': 'double',  # Slit offset in spectral dimension
                               'focus': 'double',  # Focus offset
                               'lamps': 'text',  # Lamps that are lit
                               'attenuator': 'double',  # Attenuator setting
                               'photodiode': 'double',  # Photodiode reading (cd/m^2)
                               }
    config.register.unique = ['site', 'category', 'visit', 'filter', 'arm', 'spectrograph',
                              'pfsDesignId']
    config.register.visit = ['visit', 'field', 'dateObs', 'taiObs', 'pfsDesignId',
                             'slitOffset', 'dither', 'shift', 'focus',
                             'lamps', 'attenuator', 'photodiode',
                             ]

    config.parse.translation = {'expTime': 'EXPTIME',
                                }
    config.parse.defaults = {'ccdTemp': "0",  # Added in commissioning run 3
                             }


class PfsIngestArgumentParser(IngestArgumentParser):
    """ArgumentParser for PFS ingestion

    Sets PFS-specific command-line arguments.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_argument("--pfsConfigDir",
                          help="Base directory with pfsConfig/pfsDesign files (default: with images)")


class PfsIngestConfig(IngestConfig):
    """Configuration for PfsIngestTask"""
    parse = ConfigurableField(target=PfsParseTask, doc="Parse the input headers")
    register = ConfigurableField(target=PfsRegisterTask, doc="Add the file to the registry")

    def setDefaults(self):
        super().setDefaults()
        setIngestConfig(self, False)


class PfsIngestTask(IngestTask):
    """Ingest PFS images and configs into the data repository"""
    ConfigClass = PfsIngestConfig
    ArgumentParser = PfsIngestArgumentParser
    _DefaultName = "ingestPfs"

    def findPfsConfig(self, dirName, fileInfo):
        """Find a PfsConfig file

        We search for the PfsConfig file in the following places, in order:
        * ``dirName/``
        * ``dirName/<dateObs>/pfsConfig/``

        If the PfsConfig file is not found there, we search for a PfsDesign
        file in the ``dirName`` directory, and use that to write a PfsConfig
        file in the ``dirName/<dateObs>/pfsConfig`` directory.

        Parameters
        ----------
        dirName : `str`
            (Base) directory in which the file resides.
        fileInfo : `dict`
            Key-value pairs defining the file.

        Returns
        -------
        path : `str`
            Path to the PfsConfig file.

        Raises
        ------
        RuntimeError
            If no PfsConfig or PfsDesign file can be found.
        """
        pfsDesignId = fileInfo["pfsDesignId"]
        visit = fileInfo["visit"]
        fileName = PfsConfig.fileNameFormat % (pfsDesignId, visit)

        path = os.path.join(dirName, fileName)
        if os.path.exists(path):
            return path
        dateDirName = os.path.join(dirName, fileInfo["dateObs"], "pfsConfig")
        path = os.path.join(dateDirName, fileName)
        if os.path.exists(path):
            return path

        self.log.warn(
            "Unable to find PfsConfig file for pfsDesignId=0x%016x, visit=%d; using pfsDesign",
            pfsDesignId,
            visit,
        )

        designName = os.path.join(dirName, PfsDesign.fileNameFormat % (pfsDesignId,))
        if not os.path.exists(designName):
            raise RuntimeError(f"Unable to find PfsConfig or PfsDesign for pfsDesignId=0x{pfsDesignId:016x}")
        design = PfsDesign.read(pfsDesignId, dirName)
        if not os.path.exists(dateDirName):
            os.makedirs(dateDirName)
        PfsConfig.fromPfsDesign(design, visit, design.pfiNominal).write(dateDirName)
        assert os.path.exists(path), f"Failed to write pfsConfig at expected path: {path}"
        return path

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
        infile = self.findPfsConfig(dirName, fileInfo)
        self.ingest(infile, outfile, mode=args.mode, dryrun=args.dryrun)

    def ingest(self, infile, outfile, mode="move", dryrun=False):
        """Ingest a file into the image repository.

        This is copied from LSST 18.1.0, with the addition of commit
        ``7d27e3e8694b5f62b6ae8ef1ae6a7ad35f2829ee`` to detect when we are
        re-linking the same file.

        Parameters
        ----------
        infile : `str`
            Name of input file.
        outfile : `str`
            Name of output file (file in repository).
        mode : `str`
            Mode of ingest (copy/link/move/skip).
        dryrun : `bool`
            Only report what would occur, rather than actually doing anything?

        Returns
        -------
        success : `bool`
            Whether the file was successfully ingested.
        """
        if mode == "skip":
            return True
        if dryrun:
            self.log.info("Would %s from %s to %s" % (mode, infile, outfile))
            return True
        try:
            outdir = os.path.dirname(outfile)
            if not os.path.isdir(outdir):
                try:
                    os.makedirs(outdir)
                except OSError:
                    # Silently ignore mkdir failures due to race conditions
                    if not os.path.isdir(outdir):
                        raise
            if os.path.lexists(outfile):
                if self.config.clobber:
                    os.unlink(outfile)
                else:
                    raise RuntimeError("File %s already exists; consider --config clobber=True" % outfile)

            if mode == "copy":
                lsst.pipe.tasks.ingest.assertCanCopy(infile, outfile)
                shutil.copyfile(infile, outfile)
            elif mode == "link":
                if os.path.exists(outfile):
                    if os.path.samefile(infile, outfile):
                        self.log.debug("Already linked %s to %s: ignoring" % (infile, outfile))
                    else:
                        self.log.warn("%s already has a file at the target location (%s): ignoring "
                                      "(set clobber=True to overwrite)" % (infile, outfile))
                    return False
                os.symlink(os.path.abspath(infile), outfile)
            elif mode == "move":
                lsst.pipe.tasks.ingestassertCanCopy(infile, outfile)
                os.rename(infile, outfile)
            else:
                raise AssertionError("Unknown mode: %s" % mode)
            self.log.info("%s --<%s>--> %s" % (infile, mode, outfile))
        except Exception as e:
            self.log.warn("Failed to %s %s to %s: %s" % (mode, infile, outfile, e))
            if not self.config.allowError:
                raise
            return False
        return True

    def runFile(self, infile, registry, args, pos):
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
        pos : `int`
            Position number of this file in the input list.

        Returns
        -------
        hduInfoList : `list` of `dict`
            Information from each of the HDUs to be ingested.
        """
        if self.isBadFile(infile, args.badFile):
            self.log.info("Skipping declared bad file %s", infile)
            return None
        try:
            fileInfo, hduInfoList = self.parse.getInfo(infile)
        except Exception as e:
            if not self.config.allowError:
                raise RuntimeError(f"Error parsing {infile}") from e
            self.log.warning("Error parsing %s (%s); skipping", infile, e)
            return None
        if self.isBadId(fileInfo, args.badId.idList):
            self.log.info("Skipping declared bad file %s: %s", infile, fileInfo)
            return None
        if registry is not None and self.register.check(registry, fileInfo):
            if args.ignoreIngested:
                return None
            self.log.warning("%s: already ingested: %s", infile, fileInfo)
        outfile = self.parse.getDestination(args.butler, fileInfo, infile)
        if not self.ingest(infile, outfile, mode=args.mode, dryrun=args.dryrun):
            return None
        if hduInfoList is None:
            return None
        if registry is None:
            return hduInfoList
        for info in hduInfoList:
            try:
                self.register.addRow(registry, info, dryrun=args.dryrun, create=args.create)
            except Exception as exc:
                raise IngestError(f"Failed to register file {infile}", infile, pos) from exc

        if hduInfoList[0]["category"] in ("A", "B"):
            pfsConfigDir = args.pfsConfigDir if args.pfsConfigDir is not None else os.path.dirname(infile)
            self.ingestPfsConfig(pfsConfigDir, hduInfoList[0], args)

        return None  # No further registration should be performed


class PfsPgsqlRegisterTask(PgsqlRegisterTask):
    def addRow(self, conn, info, dryrun=False, create=False, table=None):
        """Add a row to the file table (typically 'raw').

        This override of the LSST implementation removes the visit table
        because PFS doesn't use the visits table, and its creation can be
        expensive. Furthermore, the LSST implementation hard-codes the
        "OR IGNORE" part of the SQL, which doesn't work for PostgreSQL.

        Parameters
        ----------
        conn : `sqlite3.Connection`
            Database connection.
        info : `dict`
            File properties to add to database.
        dryrun : `bool`, optional
            Simulate what would happen?
        create : `bool`, optional
            Create a table if it doesn't exist? This parameter appears to be
            vestigial.
        table : `str`, optional
            Name of table.
        """
        with conn:
            if table is None:
                table = self.config.table
            sql = "INSERT INTO %s (%s) VALUES (" % (table, ",".join(self.config.columns))
            sql += ",".join([self.placeHolder] * len(self.config.columns)) + ")"
            values = [self.typemap[tt](info[col]) for col, tt in self.config.columns.items()]

            if dryrun:
                print("Would execute: '%s' with %s" % (sql, ",".join([str(value) for value in values])))
            else:
                conn.cursor().execute(sql, values)


class PfsPgsqlIngestConfig(PfsIngestConfig):
    """Configuration for PfsPgsqlIngestTask"""
    parse = ConfigurableField(target=PfsParseTask, doc="Parse the input headers")
    register = ConfigurableField(target=PfsPgsqlRegisterTask, doc="Add the file to the registry")

    def setDefaults(self):
        super().setDefaults()
        setIngestConfig(self, False)


class PfsPgsqlIngestTask(PfsIngestTask):
    ConfigClass = PfsPgsqlIngestConfig


class PfsIngestCalibsConfig(IngestCalibsConfig):
    """Configuration for PfsIngestCalibsTask"""
    def setDefaults(self):
        super().setDefaults()
        self.register.retarget(PfsCalibsRegisterTask)


class PfsIngestCalibsTask(IngestCalibsTask):
    ConfigClass = PfsIngestCalibsConfig

    def runFile(self, infile, registry, args):
        """Examine and ingest a single file

        Parameters
        ----------
        infile : `str`
            Name of file to process.
        registry : `sqlite3.Connection`
            Registry database connection.
        args : `argparse.Namespace`
            Parsed command-line arguments.

        Returns
        -------
        calibType : `str`
            Calibration type (e.g., bias, dark, flat), or None.
        hduInfoList : `list` of `dict`
            Parsed information from FITS HDUs, or None.
        """
        fileInfo, hduInfoList = self.parse.getInfo(infile)
        calibType = self.parse.getCalibType(infile)
        if calibType not in self.register.config.tables:
            self.log.warn(str("Skipped adding %s of observation type '%s' to registry "
                              "(must be one of %s)" %
                              (infile, calibType, ", ".join(self.register.config.tables))))
            return None, None
        if args.mode != 'skip':
            outfile = self.parse.getDestination(args.butler, fileInfo, infile)
            ingested = self.ingest(infile, outfile, mode=args.mode, dryrun=args.dryrun)
            if not ingested:
                self.log.warn(str("Failed to ingest %s of observation type '%s'" %
                              (infile, calibType)))
                return None, None

        if self.register.check(registry, fileInfo, table=calibType):
            if args.ignoreIngested:
                return None, None
            self.log.warn("%s: already ingested: %s" % (infile, fileInfo))

        # If we have a bias/dark/flat for 'r' or 'm' arm, ingest the alternate: the pixels are identical
        assert len(hduInfoList) == 1
        hduInfo = hduInfoList[0]
        if calibType in ("bias", "dark", "flat") and hduInfo["arm"] in ("r", "m"):
            copy = hduInfo.copy()
            copy["arm"] = "r" if copy["arm"] == "m" else "m"
            hduInfoList.append(copy)
            copyFile = self.parse.getDestination(args.butler, copy, infile)
            self.ingest(outfile, copyFile, mode="copy" if args.mode == "move" else args.mode,
                        dryrun=args.dryrun)

        return calibType, hduInfoList

    def run(self, args):
        """Ingest all specified files and add them to the registry"""
        calibRoot = args.calib if args.calib is not None else args.output
        filenameList = self.expandFiles(args.files)
        with self.register.openRegistry(calibRoot, create=args.create, dryrun=args.dryrun) as registry:
            for infile in filenameList:
                calibType, hduInfoList = self.runFile(infile, registry, args)
                if hduInfoList is None:
                    continue
                for info in hduInfoList:
                    if self.register.check(registry, info, table=calibType):
                        self.register.updateRow(registry, info, dryrun=args.dryrun, table=calibType)
                    else:
                        self.register.addRow(registry, info, dryrun=args.dryrun,
                                             create=args.create, table=calibType)
            if not args.dryrun:
                self.register.updateValidityRanges(registry, args.validity)
            else:
                self.log.info("Would update validity ranges here, but dryrun")


def _convertToDate(dateString):
    """Convert a string into a date object"""
    return dateutil.parser.parse(dateString)


class PfsCalibsRegisterTask(CalibsRegisterTask):
    def updateRow(self, conn, info, table, dryrun=False):
        """Update a row in the table

        Parameters
        ----------
        conn : `sqlite3.Connection`
            Database connection.
        info : `dict` [`str`: POD]
            File properties to record in database.
        table : `str`
            Name of table in database.
        dryrun : `bool`
            Pretend only?
        """
        if table is None:
            table = self.config.table

        columns = set([col for col in self.config.columns if col in info])
        columns.difference_update(set(self.config.unique))
        values = []
        sql = f"UPDATE {table} SET "
        sql += ", ".join(["%s=%s" % (col, self.placeHolder) for col in columns])
        values += [info[col] for col in columns]
        sql += " WHERE "
        sql += " AND ".join(["%s=%s" % (col, self.placeHolder) for col in self.config.unique])
        values += [info[col] for col in self.config.unique]

        if dryrun:
            self.log.fatal("Would execute: '%s' with %s", sql, ",".join([str(value) for value in values]))
        else:
            conn.cursor().execute(sql, values)

    def fixSubsetValidity(self, conn, table, detectorData, validity):
        """Update the validity ranges among selected rows in the registry.

        For defects and qe_curve, the products are valid from their start date until
        they are superseded by subsequent defect data.
        For other calibration products, the validity ranges are checked and
        if there are overlaps, a midpoint is used to fix the overlaps,
        so that the calibration data with whose date is nearest the date
        of the observation is used.

        @param conn: Database connection
        @param table: Name of table to be selected
        @param detectorData: Values identifying a detector (from columns in self.config.detector)
        @param validity: Validity range (days)
        """
        columns = ", ".join([self.config.calibDate, self.config.validStart, self.config.validEnd])
        sql = "SELECT id, %s FROM %s" % (columns, table)
        sql += " WHERE " + " AND ".join(col + "=?" for col in self.config.detector)
        sql += " ORDER BY " + self.config.calibDate
        cursor = conn.cursor()
        cursor.execute(sql, detectorData)
        rows = cursor.fetchall()

        try:
            valids = collections.OrderedDict([(_convertToDate(row[self.config.calibDate]), [None, None]) for
                                              row in rows])
        except Exception:
            det = " ".join("%s=%s" % (k, v) for k, v in zip(self.config.detector, detectorData))
            # Sqlite returns unicode strings, which cannot be passed through SWIG.
            self.log.warn(str("Skipped setting the validity overlaps for %s %s: missing calibration dates" %
                              (table, det)))
            return
        dates = list(valids.keys())
        if table in self.config.validityUntilSuperseded:
            # A calib is valid until it is superseded
            for thisDate, nextDate in zip(dates[:-1], dates[1:]):
                valids[thisDate][0] = thisDate - datetime.timedelta(microseconds=1)
                valids[thisDate][1] = nextDate - datetime.timedelta(microseconds=2)
            valids[dates[-1]][0] = dates[-1] - datetime.timedelta(microseconds=1)
            valids[dates[-1]][1] = _convertToDate("2037-12-31")  # End of UNIX time
            valids[dates[0]][0] = _convertToDate("1970-01-01")  # Start of UNIX time
        else:
            # A calib is valid within the validity range (in days) specified.
            for dd in dates:
                valids[dd] = [dd - datetime.timedelta(validity), dd + datetime.timedelta(validity)]
            # Fix the dates so that they do not overlap, which can cause the butler to find a
            # non-unique calib.
            midpoints = [t1 + (t2 - t1)//2 for t1, t2 in zip(dates[:-1], dates[1:])]
            for i, (date, midpoint) in enumerate(zip(dates[:-1], midpoints)):
                if valids[date][1] > midpoint:
                    nextDate = dates[i + 1]
                    valids[nextDate][0] = midpoint + datetime.timedelta(microseconds=1)
                    valids[date][1] = midpoint
            del midpoints
        del dates
        # Update the validity data in the registry
        for row in rows:
            calibDate = _convertToDate(row[self.config.calibDate])
            validStart = valids[calibDate][0].isoformat()
            validEnd = valids[calibDate][1].isoformat()
            sql = "UPDATE %s" % table
            sql += " SET %s=?, %s=?" % (self.config.validStart, self.config.validEnd)
            sql += " WHERE id=?"
            conn.execute(sql, (validStart, validEnd, row["id"]))
