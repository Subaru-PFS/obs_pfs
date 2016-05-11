#ingestImages.py '/Users/azuri/spectra/pfs/PFS/' --mode=link '/Users/azuri/spectra/pfs/raw/2016-01-12/*.fits' --loglevel 'info' -C /Users/azuri/stella-git/obs_pfs/config/pfs/ingest.py
import re
import datetime

from lsst.pipe.tasks.ingest import IngestTask, ParseTask, IngestArgumentParser, ParseConfig
import lsst.afw.image as afwImage

class PfsParseConfig(ParseConfig):
    def setDefaults(self):
        ParseConfig.setDefaults(self)
#        self.translators["date"] = "translate_date"
        self.translators["field"] = "translate_field"
#        self.defaults["filter"] = "NONE"

class PfsParseTask(ParseTask):
    ConfigClass = PfsParseConfig

    def getInfo(self, filename):
        """Get information about the image from the filename and its contents

        @param filename    Name of file to inspect
        @return File properties; list of file properties for each extension
        """
        #matches = re.search("^PFSA(\d{6})(\d)(\d).fits", filename)
        matches = re.search("PF([JLXIASPF])([ABCDS])(\d{6})(\d)(\d).fits", filename)
        if not matches:
            matches = re.search("PF([JLXIASPF])([ABCDS])-(\d{6})(\d)(\d).fits", filename)
            if not matches:
                raise RuntimeError("Unable to interpret filename: %s" % filename)
        site, category, visit, filterInt, spectrograph = matches.groups()
        if int(spectrograph) > 4:
            spectrograph = '4'
        ccd = int(spectrograph)-1
        filter = ''
        if filterInt == '0':
            filter = 'b'
        elif filterInt == '1':
            filter = 'r'
            ccd += 4
        elif filterInt == '2':
            filter = 'n'
            ccd += 8
        else:
            filter = 'm'
            ccd += 4
        arm = filter

        header = afwImage.readMetadata(filename)
        info = dict(site=site, category=category, visit=int(visit), filter=filter, arm=arm, spectrograph=int(spectrograph), ccd=int(ccd))
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
        self.log.info('PfsParseTask.translate_field: field = %s' % field)
        return re.sub(r'\W', '_', field).upper()
