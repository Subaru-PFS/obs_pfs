import os
import numpy as np

from astro_metadata_translator import fix_header

from lsst.obs.base import CameraMapper, MakeRawVisitInfo
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.daf.persistence import Policy
import lsst.utils as utils
import lsst.obs.base.yamlCamera as yamlCamera
from .translator import PfsTranslator


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
        argDict["darkTime"] = self.popFloat(md, "DARKTIME") if "DARKTIME" in md else np.nan
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
    _cameraName = "pfs"
    yamlFileList = ("PfsMapper.yaml",)  # list of yaml files to load, keeping the first occurrence
    MakeRawVisitInfoClass = PfsRawVisitInfo
    translatorClass = PfsTranslator

    def __init__(self, **kwargs):
        policyFile = Policy.defaultPolicyFile("obs_pfs", "PfsMapper.yaml", "policy")
        policy = Policy(policyFile)
        if not kwargs.get('root', None):
            pfsDataDir = os.environ.get('PFS_DATA_DIR')
            if pfsDataDir is None:
                raise RuntimeError("Either $PFS_DATA_DIR or root= must be specified")
            kwargs['root'] = os.path.join(pfsDataDir, 'PFS')
        if not kwargs.get('calibRoot', None):
            kwargs['calibRoot'] = os.path.join(kwargs['root'], 'CALIB')

        super(PfsMapper, self).__init__(policy, os.path.dirname(policyFile), **kwargs)

        # Ensure each dataset type of interest knows about the full range of keys available from the registry
        keys = {'site': str,
                'category': str,
                'field': str,
                'visit': int,
                'ccd': int,
                'filter': str,
                'arm': str,
                'spectrograph': int,
                'dateObs': str,
                'expTime': float,
                'dataType': str,
                'taiObs': str,
                'pfsDesignId': int,
                'slitOffset': float,
                }
        for name in ("raw", "pfsArm", "wlFitData", "arcLines"):
            self.mappings[name].keyDict.update(keys)

        # The order of these defineFilter commands matters as their IDs are used to generate at least some
        # object IDs (e.g. on coadds) and changing the order will invalidate old objIDs

        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter(name="UNRECOGNISED", lambdaEff=0,
                                   alias=["NONE", "None", "Unrecognised", "UNRECOGNISED",
                                          "Unrecognized", "UNRECOGNIZED", "NOTSET", ])
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

    @classmethod
    def _makeCamera(cls, policy=None, repositoryDir=None, cameraYamlFile=None):
        """Make a camera  describing the camera geometry.

        policy : ignored
        repositoryDir : ignored
        cameraYamlFile : `str`
           The full path to a yaml file to be passed to `yamlCamera.makeCamera`

        Returns
        -------
        camera : `lsst.afw.cameraGeom.Camera`
            Camera geometry.
        """
        yamlFile = os.path.join(utils.getPackageDir(cls.packageName), "pfs", "camera", "camera.yaml")
        return yamlCamera.makeCamera(yamlFile)

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

        1. Apply header patches. This should be done by the CameraMapper base
        class, but it isn't yet (DM-23959).

        2. Fix an early ADC bug.

        Since this is almost certainly an FPGA bug, I'll base the
        decision on the FPGA version number. As of 2016-12-01 the
        keyword is misnamed, so we can fix the format if the keyword
        does not exist.

        See _shiftAmpPixels() for the implementation.

        Parameters
        ----------
        item : image-like object
            Can be any of lsst.afw.image.Exposure,
            lsst.afw.image.DecoratedImage, lsst.afw.image.Image
            or lsst.afw.image.MaskedImage
            the image-like object whose header information needs to be patched.
        dataId : `dict`
            Dataset identifier

        Returns
        -------
        item : image-like object
            the input object with patched header information.
        """
        exp = super(PfsMapper, self).std_raw(item, dataId)

        filename = self.map("raw_filename", dataId).getLocations()[0]
        md = exp.getMetadata()
        fix_header(md, translator_class=PfsTranslator, filename=filename)
        try:
            dataVersion = int(md.get('W_VERSIONS_FPGA'), 16)
        except Exception:
            dataVersion = None

        if dataVersion is not None and dataVersion <= 0x0070:
            self._shiftAmpPixels(exp)

        return exp

    def std_raw_md(self, item, dataId):
        """Fixup raw header metadata problems.

        This should be done by the CameraMapper base class, but it isn't yet (DM-23959).

        Parameters
        ----------
        item : `lsst.daf.base.PropertyList`
            The raw metadata to be fixed

        dataId : `dict`
            Dataset identifier

        Returns
        -------
        item : `lsst.daf.base.PropertyList`
            The modified raw metadata.
        """
        filename = self.map("raw_filename", dataId).getLocations()[0]
        fix_header(item, translator_class=PfsTranslator, filename=filename)
        return item

    def std_rawAG(self, item, dataId):
        raw = afwImage.makeExposure(afwImage.makeMaskedImage(item.getImage()))
        raw.setMetadata(item.getMetadata())

        raw.setDetector(self._makeCamera()[f"AG{dataId['hdu']}"])

        vi = afwImage.visitInfo.VisitInfo(raw.getMetadata())
        raw.getInfo().setVisitInfo(vi)

        return raw

    def std_fiberProfiles(self, item, dataId):
        """Disable standardization for fiberProfiles

        Because it is not an Exposure.
        """
        return item

    def std_detectorMap(self, item, dataId):
        """Disable standardization for detectorMap

        Because it is not an Exposure.
        """
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

    def _defectLookup(self, dataId, dateKey='taiObs'):
        """Find the defects for a given CCD.

        This is copied from LSST 18.1.0, but using the "raw" table instead of
        "raw_visit", because PFS doesn't use the raw_visit table.

        Parameters
        ----------
        dataId : `dict`
            Dataset identifier.

        Returns
        -------
        path : `str`
            Path to the defects file or None if not available.
        """
        if self.defectRegistry is None:
            return None
        if self.registry is None:
            raise RuntimeError("No registry for defect lookup")

        ccdKey, ccdVal = self._getCcdKeyVal(dataId)

        dataIdForLookup = dict(visit=dataId["visit"], arm=dataId["arm"])
        # .lookup will fail in a posix registry because there is no template to provide.
        rows = self.registry.lookup((dateKey), ('raw'), dataIdForLookup)
        if len(rows) == 0:
            return None
        assert len(rows) == 1
        dayObs = rows[0][0]

        # Lookup the defects for this CCD serial number that are valid at the exposure midpoint.
        rows = self.defectRegistry.executeQuery(("path",), ("defect",),
                                                [(ccdKey, "?")],
                                                ("DATETIME(?)", "DATETIME(validStart)", "DATETIME(validEnd)"),
                                                (ccdVal, dayObs))
        if not rows or len(rows) == 0:
            return None
        if len(rows) == 1:
            return os.path.join(self.defectPath, rows[0][0])
        else:
            raise RuntimeError("Querying for defects (%s, %s) returns %d files: %s" %
                               (ccdVal, dayObs, len(rows), ", ".join([_[0] for _ in rows])))

    def _setFilter(self, mapping, item, dataId):
        """Set the filter object in an Exposure.  If the Exposure had a FILTER
        keyword, this was already processed during load.  But if it didn't,
        use the filter from the registry.

        Parameters
        ----------
        mapping : `lsst.obs.base.Mapping`
            Where to get the filter from.
        item : `lsst.afw.image.Exposure`
            Exposure to set the filter in.
        dataId : `dict`
            Dataset identifier.
        """
        if not (isinstance(item, afwImage.ExposureU) or isinstance(item, afwImage.ExposureI) or
                isinstance(item, afwImage.ExposureF) or isinstance(item, afwImage.ExposureD)):
            return

        if item.getFilter().getId() != afwImage.Filter.UNKNOWN:
            return

        actualId = mapping.need(["arm"], dataId)
        filterName = actualId["arm"]
        if self.filters is not None and filterName in self.filters:
            filterName = self.filters[filterName]
        item.setFilter(afwImage.Filter(filterName))
