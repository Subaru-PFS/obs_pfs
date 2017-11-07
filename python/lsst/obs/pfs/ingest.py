import os
import re

import lsst.afw.image as afwImage
from lsst.obs.pfs.pfsMapper import PfsMapper
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask
import lsst.utils

class PfsParseConfig(ParseConfig):
    def setDefaults(self):
        ParseConfig.setDefaults(self)
        self.translators["field"] = "translate_field"

class PfsParseTask(ParseTask):
    ConfigClass = PfsParseConfig

    nSpectrograph = 4
    arms = ['b', 'r', 'n', 'm']
    sites = '[JLXIASPF]'
    categories = '[ABCDS]'
    
    def getInfo(self, filename):
        """Get information about the image from the filename and its contents

        @param filename    Name of file to inspect
        @return File properties; list of file properties for each extension
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
            raise Exception('spectrograph (=%d) out of bounds 1..%d' % (spectrograph, self.nSpectrograph))

        try:
            arm = self.arms[armNum - 1]   # 1-indexed
        except IndexError:
            raise Exception('armNum (=%d) out of bounds 1..%d' % (armNum, max(self.arms.keys())))

        dataId = dict(spectrograph=spectrograph, arm=arm)
        ccd = PfsMapper.computeDetectorId(spectrograph, arm)
        self.log.debug(
            'site = <%s>, category = <%s>, visit = <%s>, spectrograph = <%d>, armNum = <%d>, arm = <%s>, ccd = <%d>'
            % (site,category,visit,spectrograph,armNum,arm,ccd)
        )

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
        """Get 'pfsConfigId' from metadata, setting the value to 0x0 if it's not present
        """
        key = "PFSCONFIGID"

        if md.exists(key):
            return md.get(key)
        else:
            return 0x0

class PfsCalibsParseTask(CalibsParseTask, PfsParseTask):
    ConfigClass = PfsParseConfig

    calibTypes = ["Bias", "Dark", "DetectorMap", "FiberTrace", "FiberFlat"]
    
    def getInfo(self, filename):
        """Get information about the image from the filename and its contents

        @param filename    Name of file to inspect
        @return File properties; list of file properties for each extension
        """
        self.log.debug('interpreting filename <%s>' % filename)

        path, name = os.path.split(filename)
        matches = re.search(r"^pfs(%s)-(\d{4}-\d{2}-\d{2})-(\d)-([brnm])(\d)\.fits$" % \
                            ("|".join(self.calibTypes)), name)
        
        calibDate = None                # one possible label for calibrations (deprecated!)
        visit0 = 0                      # another possible label;  the one in the datamodel
        if matches:
            calibType, calibDate, _, arm, spectrograph = matches.groups()
        else:
            matches = re.search(r"^pfs(%s)-(\d{6})-([brnm])(\d)\.fits$" % ("|".join(self.calibTypes)), name)
            if matches:
                calibType, visit0, arm, spectrograph = matches.groups()
                visit0 = int(visit0)

        if matches is None:
            raise RuntimeError("Unable to interpret %s as a calibration file" % filename)
        
        try:
            armNum = self.arms.index(arm) + 1
        except ValueError:
            raise RuntimeError("found unknown filter \"%s\" in %s; expected one of %s" %
                               (arm, filename, self.arms))

        spectrograph = int(spectrograph)

        if spectrograph < 1 or spectrograph > self.nSpectrograph + 1:
            raise Exception('spectrograph (=%d) out of bounds 1..%d' % (spectrograph, self.nSpectrograph))

        dataId = dict(spectrograph=spectrograph, arm=arm)
        ccd = PfsMapper.computeDetectorId(spectrograph, arm)
        self.log.debug(
            'visit0 = <%d>, spectrograph = <%d>, arm = <%s>, ccd = <%d>' %
            (visit0, spectrograph, arm, ccd)
        )

        info = dict(visit0=visit0, arm=arm, filter=arm, spectrograph=spectrograph, ccd=ccd,
                    calibDate=calibDate)

        if os.path.exists(filename):
            header = afwImage.readMetadata(filename)
            info = self.getInfoFromMetadata(header, info=info)

        return info, [info]
    
