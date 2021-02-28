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
import numpy as np
import numpy.linalg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.daf.persistence.butlerExceptions import NoResults

import lsst.ip.isr as ipIsr
from lsst.ip.isr import isrQa

___all__ = ["IsrTask", "IsrTaskConfig"]


class BrokenShutterConfig(pexConfig.Config):
    """Configuration parameters for working around a broken shutter
    """
    useAnalytic = pexConfig.Field(dtype=bool, default=False,
                                  doc="Use an analytic correction")
    window = pexConfig.Field(dtype=int, default=1,
                             doc="Number of frames +- to search for a bias")

    causalWindow = pexConfig.Field(dtype=bool, default=True,
                                   doc="Only search for later biases?")

    t_wipe = pexConfig.Field(dtype=float, default=5.23,
                             doc="Time taken to wipe the CCD (s)")
    t_read = pexConfig.Field(dtype=float, default=36.79,
                             doc="Time taken to read the CCD (s)")
    t_stare = pexConfig.Field(dtype=float, default=3.06,
                              doc="Time during readout when charge isn't being clocked CCD (s)")

    def validate(self):
        super().validate()
        if self.window > 1 and self.useAnalytic:
            raise ValueError(f"You may not specify both (window == {self.window}) > 1 and useAnalytic")


class PfsIsrTaskConfig(ipIsr.IsrTaskConfig):
    """Configuration parameters for PFS's IsrTask.

    Items are grouped in the order in which they are executed by the task.
    """
    doBrokenRedShutter = pexConfig.Field(dtype=bool, default=False,
                                         doc="Attempt to correct for a broken red shutter?")

    brokenRedShutter = pexConfig.ConfigField(
        dtype=BrokenShutterConfig, doc="Broken shutter related configuration options."
    )

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
            if self.config.doBrokenRedShutter and self.config.brokenRedShutter.useAnalytic:
                #
                # Build a model which generates the as-readout data R given
                # the (unknown) true signal S, with integration time t_exp:
                #    R_i = S_i + (t_wipe*S_red + t_stare*S_i + t_read*S_blue,i)/t_exp
                # where the second term gives the photons resulting from reading
                # the chip with the shutter open, and S_blue/red are the photons bluer (redder)
                # than S_i (i.e. sum_0^i S_i and sum_i^N-1 S_i)
                #
                # We write this as R = S + kern@S, or
                #  S = (1 + kern)^{-1} R
                #
                # N.b. we could cache ikern as it only depends on the exposure time
                # and the config parameters giving times
                #
                t_exp = result.exposure.getInfo().getVisitInfo().getExposureTime()
                t_wipe = self.config.brokenRedShutter.t_wipe
                t_read = self.config.brokenRedShutter.t_read
                t_stare = self.config.brokenRedShutter.t_stare

                N = result.exposure.getHeight()
                kern = np.zeros((N, N))
                kern += t_wipe/N*(1 - np.tri(N, k=1))
                kern += t_read/N*np.tri(N, k=-1)
                np.fill_diagonal(kern, t_stare)
                kern /= t_exp

                np.fill_diagonal(kern, 1 + kern.diagonal())  # i.e. add 1 to the diagonal
                ikern = np.linalg.inv(kern)

                self.log.info("Correcting for broken red shutter analytically")
                result.exposure.image.array[:] = ikern@result.exposure.image.array
            else:
                #
                # If we're using the `r` arm subtract a neighbouring image if it's a bias
                #
                # If the source of photons is stable, this corrects for the extra
                # light accumulated during CCD wipe and readout
                #
                # We search for the bias in a window around the visit v0 in order
                #   v0 + 1, v0 - 1, v0 + 2, v0 - 2, ...
                # (if brokenRedShutter.causalWindow is True, only use v0 + 1, v0 + 2, ...)
                #
                # Don't fix a bias; this is especially important as if we do, it'll lead
                # to an attempt to correct the bias we're using for the shutter correction
                #
                # N.b. we can't rely on the sensorRef.dataId because we adjusted the visit number
                # before possibly calling this routine recursively to fetch the desired bias
                md = sensorRef.get("raw_md")
                if md.get('DATA-TYP', "").upper() != "BIAS":
                    dataId = dict(arm=sensorRef.dataId["arm"],
                                  filter=sensorRef.dataId["arm"],  # yes, this is necessary
                                  spectrograph=sensorRef.dataId["spectrograph"])
                    dataId0 = sensorRef.dataId

                    for dv in sum([[i, -i] for i in range(1, self.config.brokenRedShutter.window + 1)], []):
                        if self.config.brokenRedShutter.causalWindow and dv < 0:
                            continue

                        exp2 = None
                        try:
                            sensorRef.dataId = dataId
                            dataId["visit"] = dataId0["visit"] + dv

                            md = sensorRef.get("raw_md")
                            if md.get('DATA-TYP', "").upper() == "BIAS":
                                self.log.info("Reading visit %d for broken red shutter correction",
                                              sensorRef.dataId["visit"])

                                exp2 = self.runDataRef(sensorRef).exposure
                        except NoResults as e:
                            self.log.warn("Unable to read %d %s%d for broken red shutter correction: %s",
                                          sensorRef.dataId["visit"], sensorRef.dataId["arm"],
                                          sensorRef.dataId["spectrograph"], e)
                        finally:
                            sensorRef.dataId = dataId0

                        if exp2:
                            self.log.info("Correcting for broken red shutter using visit %d", dataId["visit"])
                            result.exposure.maskedImage -= exp2.maskedImage

                            break

                    if exp2 is None:
                        self.log.warn("Unable to find bias frame for broken red shutter correction %s",
                                      sensorRef.dataId)

        if self.config.doWrite:
            sensorRef.put(result.exposure, "postISRCCD")
            if result.preInterpolatedExposure is not None:
                sensorRef.put(result.preInterpolatedExposure, "postISRCCD_uninterpolated")
        if result.ossThumb is not None:
            isrQa.writeThumbnail(sensorRef, result.ossThumb, "ossThumb")
        if result.flattenedThumb is not None:
            isrQa.writeThumbnail(sensorRef, result.flattenedThumb, "flattenedThumb")

        return result
