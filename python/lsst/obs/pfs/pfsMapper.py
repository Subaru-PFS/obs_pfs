import os
import numpy as np

from lsst.obs.base import CameraMapper, MakeRawVisitInfo
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.obs.pfs.pfsCamera import PfsCamera
from lsst.daf.persistence import Policy


class PfsRawVisitInfo(MakeRawVisitInfo):
    def setArgDict(self, md, argDict):
        """Fill an argument dict with arguments for makeVisitInfo and pop associated metadata

        Subclasses are expected to override this method, though the override
        may wish to call this default implementation, which:
        - sets exposureTime from "EXPTIME"
        - sets date by calling getDateAvg

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet;
            items that are used are stripped from the metadata
            (except TIMESYS, because it may apply to more than one other keyword).
        @param[in,out] argdict  a dict of arguments
        """
        super(PfsRawVisitInfo, self).setArgDict(md, argDict)
        argDict["darkTime"] = self.popFloat(md, "DARKTIME")
        #
        # Done setting argDict; check values now that all the header keywords have been consumed
        #
        argDict["darkTime"] = self.getDarkTime(argDict)

    def getDarkTime(self, argDict):
        """Retrieve the dark time from an argDict, waiting to be passed to the VisitInfo ctor"""
        darkTime = argDict.get("darkTime", float("NaN"))
        if np.isfinite(darkTime):
            return darkTime

        self.log.info("darkTime is NaN/Inf; using exposureTime")
        exposureTime = argDict.get("exposureTime", np.nan)
        if not np.isfinite(exposureTime):
            raise RuntimeError("Tried to substitute exposureTime for darkTime but it is also Nan/Inf")

        return exposureTime

    def getDateAvg(self, md, exposureTime):
        """Return date at the middle of the exposure
        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet;
            items that are used are stripped from the metadata
            (except TIMESYS, because it may apply to more than one other keyword).
        @param[in] exposureTime  exposure time (sec)
        """
        
        dateObs = self.popIsoDate(md, "DATE-OBS")
        return self.offsetDate(dateObs, 0.5*exposureTime)


class PfsMapper(CameraMapper):
    """Provides abstract-physical mapping for PFS data"""
    packageName = "obs_pfs"

    MakeRawVisitInfoClass = PfsRawVisitInfo
    
    def __init__(self, **kwargs):
        policyFile = Policy.defaultPolicyFile("obs_pfs", "PfsMapper.yaml", "policy")
        policy = Policy(policyFile)
        if not kwargs.get('root', None):
            try:
                kwargs['root'] = os.path.join(os.environ.get('PFS_DATA_DIR'), 'PFS')
            except:
                raise RuntimeError("Either $PFS_DATA_DIR or root= must be specified")
        if not kwargs.get('calibRoot', None):
            kwargs['calibRoot'] = os.path.join(kwargs['root'], 'CALIB')

        super(PfsMapper, self).__init__(policy, os.path.dirname(policyFile), **kwargs)
        
        # The order of these defineFilter commands matters as their IDs are used to generate at least some
        # object IDs (e.g. on coadds) and changing the order will invalidate old objIDs

        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter(name="UNRECOGNISED", lambdaEff=0,
                                   alias=["NONE", "None", "Unrecognised", "UNRECOGNISED",
                                          "Unrecognized", "UNRECOGNIZED", "NOTSET",])
        afwImageUtils.defineFilter(name='b', lambdaEff=477, alias=['blue', 'PFS-B'])
        afwImageUtils.defineFilter(name='r', lambdaEff=623, alias=['red', 'PFS-R'])
        afwImageUtils.defineFilter(name='n', lambdaEff=623, alias=['nearInfraRed', 'PFS-N'])
        afwImageUtils.defineFilter(name='m', lambdaEff=775, alias=['mediumResolutionRed', 'PFS-M'])
        #
        # self.filters is used elsewhere, and for now we'll set it
        #
        self.filters = {}
        for f in ["b",
                  "r",
                  "n",
                  "m",
                  "NONE",
                  "UNRECOGNISED"
                  ]:
            # Get the canonical name -- see #2113
            self.filters[f] = afwImage.Filter(afwImage.Filter(f).getId()).getName()
        self.defaultFilterName = "UNRECOGNISED"

        #
        # The number of bits allocated for fields in object IDs, appropriate for
        # the default-configured Rings skymap.
        #
        # This shouldn't be the mapper's job at all; see #2797.

        PfsMapper._nbit_id = 64

    def _makeCamera(self, policy, repositoryDir):
        """Make a camera (instance of lsst.afw.cameraGeom.Camera) describing the camera geometry
        """
        return PfsCamera()

    @staticmethod
    def _flipChipsLR(exp, wcs, dataId, dims=None):
        flipLR, flipTB = (True, False)
        exp.setMaskedImage(afwMath.flipImage(exp.getMaskedImage(), flipLR, flipTB))

        return exp

    def _shiftAmpPixels(self, rawExp):
        """Shift pixels in early raw frames.

        Args
        ----
        rawExp : an Exposure for a raw image.

        Early ADC versions insert a spurious pixel at the beginning of
        the readout. This affects all rows, pushing the last overscan
        pixel of a given row into the first leadin pixel on the
        following row.

        We strip out the 0th pixel, and leave the last one untouched.

        The Exposure's Image is modified in place.
        """

        imArr = rawExp.getMaskedImage().getImage().getArray()

        for a in rawExp.getDetector():
            yslice, xslice = a.getRawBBox().getSlices()
            ampIm = imArr[yslice, xslice]
            ampshape = ampIm.shape

            # Don't bother being tricky: make a copy.
            ampPixels = ampIm.flatten()
            ampPixels[:-1] = ampPixels[1:]

            imArr[yslice, xslice] = ampPixels.reshape(ampshape)

    def std_raw(self, item, dataId):
        """Fixup raw image problems.

        1. Fix an early ADC bug.

        Since this is almost certainly an FPGA bug, I'll base the
        decision on the FPGA version number. As of 2016-12-01 the
        keyword is misnamed, so we can fix the format if the keyword
        does not exist.

        See _shiftAmpPixels() for the implementation.
        """

        exp = super(PfsMapper, self).std_raw(item, dataId)

        md = exp.getMetadata()
        try:
            dataVersion = int(md.get('W_VERSIONS_FPGA'), 16)
        except Exception:
            dataVersion = 0

        if dataVersion <= 0x0070:
            self._shiftAmpPixels(exp)

        return exp

    def std_fibertrace(self, item, dataId):  # needed to stop the butler generating a version
        return item

    def map_linearizer(self, dataId, write=False):
        """Map a linearizer."""
        raise RuntimeError("No linearizer available.")

    def _extractAmpId(self, dataId):
        ampId = (self._extractDetectorName(dataId), 0, 0)
        return ampId

    def _extractDetectorName(self, dataId):
        arm = self._getRegistryValue(dataId, "arm")
        spectrograph = self._getRegistryValue(dataId, "spectrograph")
        if arm == 'm':
            arm = 'r'
        return "%s%d" % (arm, spectrograph)

    @staticmethod
    def computeDetectorId(spectrograph, arm):
        """
        Return a DetectorId in the range [0,11] with
           blue                 = {0,3,6,9},
           red/mediumResolution = {1,4,7,10}
           NIR                  = {2,5,8,11}
        corresponding to obs_pfs/pfs/camera/camera.py::config.detectorList[1].id
        """

        return 3*(spectrograph - 1) + dict(b=0, r=1, m=1, n=2)[arm]

    def _extractDetectorId(self, dataId):
        # The returned DetectorId is in the range [0,11] with blue={0,3,6,9},
        # red=mediumResolution={1,4,7,10}, and NIR={2,5,8,11} corresponding to
        # obs_pfs/pfs/camera/camera.py::config.detectorList[1].id
        #
        # If we wanted to get rid of this function we would need to modify
        # genDefectList, genDefectFits, and the defect lookup function
        arm = self._getRegistryValue(dataId, "arm")
        spectrograph = self._getRegistryValue(dataId, "spectrograph")

        return self.computeDetectorId(spectrograph, arm)

    def _getCcdKeyVal(self, dataId):
        """Return CCD key and value used to look a defect in the defect registry

        The default implementation simply returns ("ccd", full detector name)
        """
        return ("ccd", self._extractDetectorName(dataId))

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        @param dataId (dict) Data identifier with visit, ccd
        """
        pathId = self._transformId(dataId)
        visit = pathId['visit']
        ccd = self._extractDetectorId(dataId)
        return visit*200 + ccd

    def bypass_ccdExposureId(self, datasetType, pythonType, location, dataId):
        return self._computeCcdExposureId(dataId)

    def bypass_pfsConfig(self, datasetType, pythonType, location, dataId):
        from pfs.datamodel.pfsConfig import PfsConfig

        pfsConfigId = dataId.get('pfsConfigId', 0x0)

        for pathName in location.locationList:
            pathName = os.path.join(location.storage.root, pathName)
            dirName = os.path.dirname(pathName)

            try:
                pfsConfig = PfsConfig.read(pfsConfigId, dirName=dirName)
            except Exception as e:
                msg = str(e)
                if isinstance(e, IOError):
                    msg = msg.replace(r"[Errno 2] ", "")  # it isn't an error, just a warning
            else:
                return pfsConfig

        if pfsConfigId == 0x0:
            self.log.warn("%s" % (msg,))
            self.log.debug("%s" % (dirName,))
            self.log.info("Creating dummy PfsConfig for pfsConfigId == 0x%x" % (pfsConfigId))

            from pfs.datamodel.pfsConfig import PfsConfig
            return makeDummyPfsConfig(0, pfsConfigId)

        raise RuntimeError("Unable to read pfsConfig for %s: %s" % (dataId.items(), e))

    #
    # This is disabled due to butler changes that mean it is silently ignored
    # when the file requested doesn't exist.  Instead we do disgusting things
    # parsing filenames in PfsArmIO.readFits
    #
    def XXXbypass_pfsArm(self, datasetType, pythonType, location, dataId):
        from pfs.datamodel.pfsArm import PfsArm

        for pathName in location.locationList:
            pathName = os.path.join(location.storage.root, pathName)
            dirName = os.path.dirname(pathName)

            arm = self._getRegistryValue(dataId, "arm")
            spectrograph = self._getRegistryValue(dataId, "spectrograph")
            visit = self._getRegistryValue(dataId, "visit")

            pfsArm = PfsArm(visit, spectrograph, arm)
            try:
                pfsArm.read(dirName=dirName, setPfsConfig=False)
            except Exception as e:
                pass
            else:
                return pfsArm

        raise RuntimeError("Unable to read pfsArm for %s: %s" % (dataId.items(), e))

    def getDetectorId(self, dataId):
        return self._extractDetectorId(dataId)

    def getDetectorName(self, dataId):
        return self._extractDetectorName(dataId)

    @classmethod
    def getCameraName(cls):
        return "pfs"

    def _getRegistryValue(self, dataId, k):
        """Return a value from a dataId, or look it up in the registry if it isn't present"""
        if k in dataId:
            return dataId[k]
        else:
            dataType = "bias" if "taiObs" in dataId else "raw"

            try:
                return self.queryMetadata(dataType, [k], dataId)[0][0]
            except IndexError:
                raise RuntimeError("Unable to lookup %s in \"%s\" registry for dataId %s" %
                                   (k, dataType, dataId))

    def std_detectormap(self, item, dataId):
        return item


def assemble_pfsArm(dataId, componentInfo, cls):
    """Called by the butler to construct the composite type "pfsArm"

    This returns a `pfs.datamodel.PfsArm` instead of a
    `pfs.drp.stella.SpectrumSet`, as that includes the ``pfsConfig``.
    The result can be converted to a ``SpectrumSet`` by:

        spectra = pfs.drp.stella.SpectrumSet.fromPfsArm(pfsArm)

    Parameters
    ----------
    dataId : `lsst.daf.persistence.dataId.DataId`
        the data ID
    componentInfo : `dict`
        dict containing the components, as defined by the composite definition in the mapper policy
    cls : 'object'
        unused

    Returns
    -------
    pfsArm : `pfs.datamodel.PfsArm`
        The assembled pfsArm.
    """
    pfsArm = componentInfo['spectra'].obj.toPfsArm(dataId)
    pfsArm.pfsConfig = componentInfo['pfsConfig'].obj

    return pfsArm


def disassemble_pfsArm(dataId, obj, componentInfo):
    """Called by the butler to deconstruct the composite type "pfsArm"
    """
    componentInfo['spectra'].obj = obj
    # Dummy PfsConfig; this has never worked well, and is in need of a redesign
    componentInfo['pfsConfig'].obj = makeDummyPfsConfig(len(obj), expId=dataId["visit"])


def makeDummyPfsConfig(numFibers, pfiDesignId=0x0, expId=0):
    """Create a dummy PfsConfig

    The current pipeline creates these out of thin air, whereas in reality
    they will be provided for us along with the raw exposures; so this is
    just a placeholder until we can make things more realistic.
    """
    from pfs.drp.stella.pfsConfigIO import PfsConfigIO

    fiberId = np.arange(numFibers, dtype=int)
    floats = np.zeros(numFibers, dtype=float)
    idents = np.ones(numFibers, dtype=int)
    patches = ["x"]*numFibers
    mags = np.zeros((numFibers, 0), dtype=float)
    filters = [[]]*numFibers
    pairs = np.zeros((numFibers, 2), dtype=float)

    # Dummy PfsConfig; this has never worked well, and is in need of a redesign
    return PfsConfigIO(0x0, expId, 0.0, 0.0, fiberId, idents, patches, floats, floats,
                       idents, idents, idents, mags, filters, pairs, pairs)
