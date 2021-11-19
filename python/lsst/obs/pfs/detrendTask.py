from lsst.ip.isr import IsrTask
from lsst.pipe.tasks.repair import RepairTask
from lsst.pex.config import Config, ConfigurableField, Field
from lsst.pipe.base import CmdLineTask


class DetrendConfig(Config):
    isr = ConfigurableField(target=IsrTask, doc="Instrumental signature removal")
    doRepair = Field(dtype=bool, default=True, doc="Repair artifacts?")
    repair = ConfigurableField(target=RepairTask, doc="Task to repair artifacts")
    windowed = Field(dtype=bool, default=False,
                     doc="Reduction of windowed data, for real-time acquisition? Implies "
                     "overscanFitType=MEDIAN")

    def validate(self):
        if self.windowed:
            self.isr.windowed = True
        super().validate()


class DetrendTask(CmdLineTask):
    _DefaultName = "detrend"
    ConfigClass = DetrendConfig

    def __init__(self, *args, **kwargs):
        CmdLineTask.__init__(self, *args, **kwargs)
        self.makeSubtask("isr")
        self.makeSubtask("repair")

    def runDataRef(self, dataRef):
        exposure = self.isr.runDataRef(dataRef).exposure  # Should do ISR and CCD assembly

        if self.config.doRepair:
            #
            # We need a PSF, so get one from the config
            #
            modelPsfConfig = self.config.repair.interp.modelPsf
            psf = modelPsfConfig.apply()
            exposure.setPsf(psf)
            #
            # Do the work
            #
            self.repair.run(exposure)

        dataRef.put(exposure, "calexp")

    def _getConfigName(self):
        return None

    def _getMetadataName(self):
        return None

    def _getEupsVersionsName(self):
        return None
