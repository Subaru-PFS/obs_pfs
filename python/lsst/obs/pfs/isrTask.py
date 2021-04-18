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
import lsst.log
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.daf.persistence.butlerExceptions import NoResults

import lsst.afw.math as afwMath
import lsst.ip.isr as ipIsr
from lsst.ip.isr import isrQa

___all__ = ["IsrTask", "IsrTaskConfig"]


class BrokenShutterConfig(pexConfig.Config):
    """Configuration parameters for working around a broken shutter
    """
    brokenShutterList = pexConfig.ListField(dtype=str, default=["r1"],
                                            doc="List of cameras with broken shutters (e.g. ['r1'])")

    checkParallelOverscan = pexConfig.Field(dtype=bool, default=False,
                                            doc="Use parallel overscan to guess if the shutter is broken?")
    maximumAllowedParallelOverscanFraction = \
        pexConfig.Field(dtype=float, default=0.01,
                        doc="Largest permitted fraction of flux seen "
                        "in parallel overscan if shutter's working")

    useAnalytic = pexConfig.Field(dtype=bool, default=False,
                                  doc="Use an analytic correction")
    window = pexConfig.Field(dtype=int, default=1,
                             doc="Number of frames +- to search for a bias")

    causalWindow = pexConfig.Field(dtype=bool, default=True,
                                   doc="Only search for later biases?")

    doCacheCorrections = pexConfig.Field(dtype=bool, default=False,
                                         doc="Use a cache of analytic correction matrices?")
    cacheSize = pexConfig.Field(dtype=int, default=3,
                                doc="Number of analytic correction matrices to cache")

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
        if not self.doSaturationInterpolation and "SAT" in self.maskListToInterpolate:  # fixed on LSST master
            self.maskListToInterpolate.remove("SAT")
#
# Code to handle the broken shutter matrix matrix inverse cache
#


class BrokenShutterKernelCache:
    __cache = {}
    __maxsize = 0
    log = lsst.log.getLogger("BrokenShutterKernelCache")

    def __init__(self, maxsize=0):
        self.setCacheSize(maxsize)

    @classmethod
    def compute(cls, mat, t_exp):
        """Compute and cache a new kernel

        Kernels are cached and looked up by ``int(t_exp + 0.5)``

        Parameters:
           mat : `np.array`
             The matrix to be inverted to give a new kernel
           t_exp : `float`
             The exposure time associated with `mat`

        Return:
           The new kernel
        """
        key = int(t_exp + 0.5)
        if key not in cls.__cache:
            if cls.__maxsize > 0 and len(cls.__cache) >= cls.__maxsize:  # need to delete an old kernel
                k_min = None                                            # the least popular kernel
                n_min = 0
                for k, v in cls.__cache.items():
                    if k_min is None or v[1] < n_min:
                        k_min = k
                        n_min = v[1]

                if k_min is not None:
                    cls.log.info(f"Clearing cached kernel for key {k_min}")
                    del cls.__cache[k_min]

            cls.log.info(f"Computing ikern for key {key} ({t_exp}s)")
            ikern = np.linalg.inv(mat)
            cls.__cache[key] = [ikern, 0]

        cls.__cache[key][1] += 1
        return cls.__cache[key][0]

    @classmethod
    def setCacheSize(cls, maxsize):
        """Set the maximum number of cached kernels

        N.b. Kernels are not freed until you next
        compute a new one that doesn't fit in the cache;
        see also clear()

        Parameters:
        maxsize : `int`
           The maximum number of saved kernels

        Returns:
           The old value of maxsize
        """
        old = cls.__maxsize
        cls.__maxsize = maxsize

        return old

    @classmethod
    def clear(cls):
        """Empty the cache"""
        cls.__cache = {}

    @classmethod
    def show(cls):
        """Print some statistics about the cache"""
        print("BrokenShutterKernelCache:")
        print(f"Cache size: {'unlimited' if cls.__maxsize <= 0 else cls.__maxsize}")
        print("key   nUsed")
        for k, c in cls.__cache.items():
            print(f"{k:4}  {c[1]}")


try:
    brokenShutterKernelCache
except NameError:
    brokenShutterKernelCache = BrokenShutterKernelCache()


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

        if self.config.doBrokenRedShutter and \
           ccdExposure.getDetector().getName() in self.config.brokenRedShutter.brokenShutterList:
            if not self.config.brokenRedShutter.useAnalytic:
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

    def run(self, exposure, *args, **kwargs):
        doBrokenRedShutter = self.config.doBrokenRedShutter and \
            exposure.getDetector().getName() in self.config.brokenRedShutter.brokenShutterList

        if doBrokenRedShutter and self.config.brokenRedShutter.checkParallelOverscan:
            excessPOscan = []
            meanFlux = []
            for amp in exposure.getDetector():
                # N.b. MEANCLIP is broken for U16 data
                data = afwMath.makeStatistics(exposure[amp.getRawDataBBox()].convertF().image,
                                              afwMath.MEANCLIP).getValue()
                pOscan = afwMath.makeStatistics(exposure[amp.getRawVerticalOverscanBBox()].image,
                                                afwMath.MEDIAN).getValue()
                sOscan = afwMath.makeStatistics(exposure[amp.getRawHorizontalOverscanBBox()].image,
                                                afwMath.MEDIAN).getValue()

                meanFlux.append(data)
                excessPOscan.append(pOscan - sOscan)

            flux = np.mean(meanFlux)
            excess = np.mean(excessPOscan)
            if excess <= self.config.brokenRedShutter.maximumAllowedParallelOverscanFraction*flux:
                doBrokenRedShutter = False

            self.log.info("Checking parallel overscan: %g %s %g*%g %s",
                          excess, (">" if doBrokenRedShutter else "<="),
                          self.config.brokenRedShutter.maximumAllowedParallelOverscanFraction, flux,
                          ("" if doBrokenRedShutter else "; assuming not broken"))

        results = super().run(exposure, *args, **kwargs)
        exposure = results.exposure

        if doBrokenRedShutter and self.config.brokenRedShutter.useAnalytic:
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
            # and the config parameters specifying times for the readout phases
            #
            t_exp = exposure.getInfo().getVisitInfo().getExposureTime()
            t_wipe = self.config.brokenRedShutter.t_wipe
            t_read = self.config.brokenRedShutter.t_read
            t_stare = self.config.brokenRedShutter.t_stare

            if t_exp == 0:          # we can't correct a bias
                self.log.debug("Not correcting bias for broken red shutter analytically")
            else:
                brokenShutterKernelCache.setCacheSize(self.config.brokenRedShutter.cacheSize)

                N = exposure.getHeight()
                kern = np.zeros((N, N))
                kern += t_wipe/N*(1 - np.tri(N, k=1))  # i.e. upper tri
                kern += t_read/N*np.tri(N, k=-1)
                np.fill_diagonal(kern, t_stare)
                kern /= t_exp

                np.fill_diagonal(kern, 1 + kern.diagonal())  # i.e. add 1 to the diagonal

                ikern = brokenShutterKernelCache.compute(kern, t_exp)

                if not self.config.brokenRedShutter.doCacheCorrections:
                    brokenShutterKernelCache.clear()

                self.log.info("Correcting for broken red shutter analytically")
                exposure.image.array[:] = ikern@exposure.image.array
                exposure.variance.array[:] = (ikern**2)@exposure.variance.array

        return results
