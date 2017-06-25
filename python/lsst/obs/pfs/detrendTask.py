from lsst.obs.subaru.isr import SubaruIsrTask
from lsst.pex.config import Config, ConfigurableField
from lsst.pipe.base import CmdLineTask

class DetrendConfig(Config):
    isr = ConfigurableField(target=SubaruIsrTask, doc="Instrumental signature removal")

class DetrendTask(CmdLineTask):
    _DefaultName = "detrend"
    ConfigClass = DetrendConfig

    def __init__(self, *args, **kwargs):
        CmdLineTask.__init__(self, *args, **kwargs)
        self.makeSubtask("isr")

    def run(self, dataRef):
        exp = self.isr.runDataRef(dataRef).exposure  # Should do ISR and CCD assembly
        dataRef.put(exp, "calexp")

    def _getConfigName(self):
        return None
    def _getMetadataName(self):
        return None
    def _getEupsVersionsName(self):
        return None
