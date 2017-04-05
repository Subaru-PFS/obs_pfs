import re
import os
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig
import lsst.afw.image as afwImage

class PfsParseConfig(ParseConfig):
    def setDefaults(self):
        ParseConfig.setDefaults(self)
        self.translators["field"] = "translate_field"

class PfsParseTask(ParseTask):
    ConfigClass = PfsParseConfig

    def getInfo(self, filename):
        """Get information about the image from the filename and its contents

        @param filename    Name of file to inspect
        @return File properties; list of file properties for each extension
        """
        self.log.debug('interpreting filename <%s>' % filename)
        minSpectrograph = 1
        maxSpectrograph = 4
        minArmNum = 1
        maxArmNum = 4
        maxCCDs = 4
        arms = ['b', 'r', 'n', 'm']
        sites = '[JLXIASPF]'
        categories = '[ABCDS]'
        matches = re.search("PF("+sites+")("+categories+")(\d{6})(\d)(\d).fits", filename)
        if not matches:
            """ Try old format """
            matches = re.search("PF("+sites+")("+categories+")-(\d{6})(\d)(\d).fits", filename)
            if not matches:
                raise RuntimeError("Unable to interpret filename: %s" % filename)
        site, category, visit, spectrograph, armNum = matches.groups()
        if int(spectrograph) == 9:
            spectrograph = '2'
        if int(spectrograph) < minSpectrograph or int(spectrograph) > maxSpectrograph:
            message = 'spectrograph (=',spectrograph,') out of bounds'
            raise Exception(message)
        ccd = int(spectrograph)-1
        if int(armNum) == minArmNum:
            """do nothing"""
        elif (int(armNum) == minArmNum + 1) or (int(armNum) == minArmNum + 3):
            ccd += maxCCDs
        elif int(armNum) == minArmNum + 2:
            ccd += 2*maxCCDs
        else:
            message = 'arm number (=',armNum,') out of bounds [',minArmNum,'...',maxArmNum,']'
            raise Exception(message)

        filter = arms[ int(armNum) - minArmNum ]
        arm = filter
        self.log.debug('site = <%s>, category = <%s>, visit = <%s>, spectrograph = <%s>, armNum = <%s>, filter = <%s>, ccd = %d' % (site, category, visit, spectrograph, armNum, filter, ccd))

        info = dict(site=site, category=category, visit=int(visit, base=10), filter=filter, arm=arm, spectrograph=int(spectrograph), ccd=ccd)
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
