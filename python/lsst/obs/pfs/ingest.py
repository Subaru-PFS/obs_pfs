import re
import os
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig
import lsst.afw.image as afwImage
from lsst.obs.pfs.pfsMapper import PfsMapper
import lsst.utils

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
        arms = ['b', 'r', 'n', 'm']
        sites = '[JLXIASPF]'
        categories = '[ABCDS]'
        matches = re.search("PF("+sites+")("+categories+")(\d{6})(\d)(\d).fits", filename)
        if not matches:
            # Try old format
            matches = re.search("PF("+sites+")("+categories+")-(\d{6})(\d)(\d).fits", filename)
            if not matches:
                raise RuntimeError("Unable to interpret filename: %s" % filename)
        site, category, visit, spectrograph, armNum = matches.groups()

        # spectrograph '2' was called '9' at JHU, at least in the early days, so if we find
        # a spectrograph '9' we re-assign its number to '2'
        if int(spectrograph) == 9:
            spectrograph = '2'

        if int(spectrograph) < minSpectrograph or int(spectrograph) > maxSpectrograph:
            raise Exception('spectrograph (=$s) out of bounds' % spectrograph)

        if int(armNum) < minArmNum or int(armNum) > maxArmNum:
            raise Exception('armNum (=$s) out of bounds' % armNum)

        filter = arms[ int(armNum) - minArmNum ]
        arm = filter

        dataId = {'spectrograph': int(spectrograph), 'arm': arm}
        mapper = PfsMapper(root=lsst.utils.getPackageDir('drp_stella'))
        ccd = mapper.getDetectorId(dataId)
        self.log.info(
            'site = <%s>, category = <%s>, visit = <%s>, spectrograph = <%s>, armNum = <%s>, filter = <%s>, ccd = <%d>'
            % (site,category,visit,spectrograph,armNum,filter,ccd)
        )

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
