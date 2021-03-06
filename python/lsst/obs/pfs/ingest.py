import os
import re
from functools import partialmethod
import datetime
import dateutil.parser
import collections

from astro_metadata_translator import fix_header

import lsst.afw.image as afwImage
from lsst.obs.pfs.pfsMapper import PfsMapper
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig, IngestTask, IngestConfig, IngestArgumentParser
from lsst.pipe.tasks.ingest import RegisterTask
from lsst.pipe.tasks.ingestCalibs import IngestCalibsTask, IngestCalibsConfig
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask, CalibsRegisterTask
from lsst.pex.config import Field
from lsst.afw.fits import readMetadata
from lsst.obs.pfs.utils import getLamps

from pfs.datamodel.pfsConfig import PfsConfig, PfsDesign
from .translator import PfsTranslator

__all__ = ["PfsParseConfig", "PfsParseTask", "PfsIngestTask", "PfsIngestCalibsTask"]


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
        dateObs = md.get("DATE-OBS")
        if "T" in dateObs:
            return dateObs
        time = md.get("UT")
        return f"{dateObs}T{time}"


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
    def addVisits(self, conn, dryrun=False, table=None):
        """Generate the visits table (typically 'raw_visits') from the
        file table (typically 'raw').

        @param conn    Database connection
        @param table   Name of table in database
        """
        if table is None:
            table = self.config.table
        sql = f"INSERT INTO {table}_visit SELECT * FROM (SELECT "
        sql += ",".join(self.config.visit)
        sql += f" FROM {table} GROUP BY visit HAVING ROWID = MIN(ROWID)) AS vv1"  # Take first row as standard
        sql += " WHERE NOT EXISTS "
        sql += f"(SELECT vv2.visit FROM {table}_visit AS vv2 WHERE vv1.visit = vv2.visit)"
        if dryrun:
            self.log.fatal("Would execute: %s", sql)
        else:
            conn.cursor().execute(sql)


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
                                'taiObs': 'DATE-OBS',
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
        visit = fileInfo["visit"]
        fileName = PfsConfig.fileNameFormat % (pfsDesignId, visit)
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
        PfsConfig.fromPfsDesign(design, visit, design.pfiNominal).write(dirName)
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
            if hduInfoList is None:
                return None
            pfsConfigDir = args.pfsConfigDir if args.pfsConfigDir is not None else os.path.dirname(infile)
            self.ingestPfsConfig(pfsConfigDir, hduInfoList[0], args)
        except Exception:
            import traceback
            traceback.print_exc()
            raise
        return hduInfoList


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
        if args.calibType is None:
            calibType = self.parse.getCalibType(infile)
        else:
            calibType = args.calibType
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
                valids[thisDate][0] = thisDate
                valids[thisDate][1] = nextDate - datetime.timedelta(microseconds=1)
            valids[dates[-1]][0] = dates[-1]
            valids[dates[-1]][1] = _convertToDate("2037-12-31")  # End of UNIX time
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
