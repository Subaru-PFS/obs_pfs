# This file is part of obs_pfs
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.daf.persistence.butlerExceptions import NoResults

import lsst.ip.isr as ipIsr
from lsst.ip.isr import RunIsrTask, RunIsrConfig

___all__ = ["IsrTask", "IsrTaskConfig", "RunIsrTask", "RunIsrConfig"]


class PfsIsrTaskConfig(ipIsr.IsrTaskConfig):
    """Configuration parameters for PFS's IsrTask.

    Items are grouped in the order in which they are executed by the task.
    """
    doBrokenRedShutter = pexConfig.Field(dtype=bool, default=False,
                                         doc="Attempt to correct for a broken red shutter?")

    def validate(self):
        super().validate()


class PfsIsrTask(ipIsr.IsrTask):
    """Apply common instrument signature correction algorithms to a raw frame.

    We run the vanilla ISR, with the exception of optionally correcting for
    the effects of a broken red-arm shutter on a PFS spectrograph
    methods have been split into subtasks that can be redirected
    appropriately.

    Parameters
    ----------
    args : `list`
        Positional arguments passed to the Task constructor. None used at this time.
    kwargs : `dict`, optional
        Keyword arguments passed on to the Task constructor. None used at this time.
    """
    ConfigClass = PfsIsrTaskConfig
    _DefaultName = "isr"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):
        """Perform instrument signature removal on a ButlerDataRef of a Sensor.

        This method contains the `CmdLineTask` interface to the ISR
        processing.  All IO is handled here, freeing the `run()` method
        to manage only pixel-level calculations.  The steps performed
        are:
        - Read in necessary detrending/isr/calibration data.
        - Process raw exposure in `run()`.
        - Persist the ISR-corrected exposure as "postISRCCD" if
          config.doWrite=True.

        Parameters
        ----------
        sensorRef : `daf.persistence.butlerSubset.ButlerDataRef`
            DataRef of the detector data to be processed

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with component:
            - ``exposure`` : `afw.image.Exposure`
                The fully ISR corrected exposure.

        Raises
        ------
        RuntimeError
            Raised if a configuration option is set to True, but the
            required calibration data does not exist.

        """
        self.log.info("Performing ISR on sensor %s" % (sensorRef.dataId))

        ccdExposure = sensorRef.get(self.config.datasetType)

        camera = sensorRef.get("camera")
        if camera is None and self.config.doAddDistortionModel:
            raise RuntimeError("config.doAddDistortionModel is True "
                               "but could not get a camera from the butler")
        isrData = self.readIsrData(sensorRef, ccdExposure)

        result = self.run(ccdExposure, camera=camera, **isrData.getDict())

        if self.config.doBrokenRedShutter and sensorRef.dataId["arm"] == 'r':
            #
            # If we're using the `r` arm subtract the next image if it's a bias
            #
            # If the source of photons is stable, this corrects for the extra
            # light accumulated during CCD wipe and readout
            #
            # Don't fix a bias; this is especially important as if we do, it'll lead
            # to an attempt to correct the bias we're using for the shutter correction
            #
            # N.b. we can't rely on the sensorRef.dataId because we adjusted the visit number
            md = sensorRef.get("raw_md")
            if md.get('DATA-TYP', "").upper() != "BIAS":
                sensorRef.dataId["visit"] += 1

                exp2 = None
                try:
                    md = sensorRef.get("raw_md")
                    if md.get('DATA-TYP', "").upper() == "BIAS":
                        self.log.info("Reading visit %d for broken red shutter correction",
                                      sensorRef.dataId["visit"])

                        exp2 = self.runDataRef(sensorRef).exposure
                except NoResults as e:
                    self.log.warn("Unable to read next frame for broken red shutter correction %s",
                                  sensorRef.dataId)
                else:
                    sensorRef.dataId["visit"] -= 1

                if exp2:
                    self.log.info("Correcting for broken red shutter")
                    result.exposure.maskedImage -= exp2.maskedImage

        if self.config.doWrite:
            sensorRef.put(result.exposure, "postISRCCD")
            if result.preInterpolatedExposure is not None:
                sensorRef.put(result.preInterpolatedExposure, "postISRCCD_uninterpolated")
        if result.ossThumb is not None:
            isrQa.writeThumbnail(sensorRef, result.ossThumb, "ossThumb")
        if result.flattenedThumb is not None:
            isrQa.writeThumbnail(sensorRef, result.flattenedThumb, "flattenedThumb")

        return result
