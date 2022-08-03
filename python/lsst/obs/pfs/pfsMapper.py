import os
import re
import numpy as np

from astro_metadata_translator import fix_header

from lsst.geom import BoxI, PointI, ExtentI
from lsst.obs.base import CameraMapper, MakeRawVisitInfo
import lsst.afw.image as afwImage
from lsst.afw.fits import FitsError
import lsst.afw.math as afwMath
from lsst.daf.persistence import LogicalLocation, Policy
import lsst.utils as utils
import lsst.obs.base.yamlCamera as yamlCamera
from lsst.obs.base import FilterDefinition, FilterDefinitionCollection
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


pfsFilterDefinitions = FilterDefinitionCollection(
    FilterDefinition(band="b", physical_filter="b", lambdaEff=500),
    FilterDefinition(band="r", physical_filter="r", lambdaEff=800),
    FilterDefinition(band="m", physical_filter="m", lambdaEff=800),
    FilterDefinition(band="n", physical_filter="n", lambdaEff=1100),
)


class PfsMapper(CameraMapper):
    """Provides abstract-physical mapping for PFS data"""
    packageName = "obs_pfs"
    _cameraName = "pfs"
    yamlFileList = ("PfsMapper.yaml",)  # list of yaml files to load, keeping the first occurrence
    MakeRawVisitInfoClass = PfsRawVisitInfo
    translatorClass = PfsTranslator
    filterDefinitions = pfsFilterDefinitions
    _gen3instrument = "lsst.obs.pfs.PrimeFocusSpectrograph"

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

        #
        # The number of bits allocated for fields in object IDs, appropriate for
        # the default-configured Rings skymap.
        #
        # This shouldn't be the mapper's job at all; see #2797.
        PfsMapper._nbit_id = 64

    @classmethod
    def addFilters(cls):
        pfsFilterDefinitions.defineFilters()

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

    @staticmethod
    def lookupHdu(read, imType="IMAGE", md=None):
        """Return the appropriate HDU

        Parameters
        ----------
        read: `int`
            The desired read of the device (1-indexed)
        imType: `str`
            IMAGE or REF
        md: `lsst.daf.base.PropertyList`
            Header from file (optional)
            Used to validate read, and handle extra reset HDUs if present
        """

        if md:
            nread = md.get("W_H4NRED")

            if nread is not None and not (1 <= read <= nread):
                raise RuntimeError(f"The read must be in the range 1..{nread}; saw {read}")

        hdu = 2*read - 1 + dict(IMAGE=0, REF=1)[imType]

        return hdu

    def bypass_datasetExists_raw(self, datasetType, pythonType, butlerLocation, dataId):
        return self.bypass_raw(datasetType, pythonType, butlerLocation, dataId, checkExistence=True)

    def bypass_raw(self, datasetType, pythonType, butlerLocation, dataId, checkExistence=False):
        if dataId['arm'] in "bmr":
            raise IOError("Please ignore this bypass function for [bmr] data")

        def readData(datasetType, butlerLocation, dataId, hdu=1, getMetadata=False, checkExistence=False):
            additionalData = butlerLocation.getAdditionalData()
            locationString = butlerLocation.getLocations()[0]
            locStringWithRoot = os.path.join(butlerLocation.getStorage().root, locationString)
            logLoc = LogicalLocation(locStringWithRoot, additionalData)
            # test for existence of file, ignoring trailing [...]
            # because that can specify the HDU or other information
            filePath = re.sub(r"(\.fits(.[a-zA-Z0-9]+)?)(\[.+\])$", r"\1", logLoc.locString())

            if not os.path.exists(filePath):  # handle reruns
                cfg = butlerLocation.storage.getRepositoryCfg(butlerLocation.storage.root)
                for p in cfg.parents:
                    nfilePath = os.path.join(p.root, locationString)
                    if os.path.exists(nfilePath):
                        filePath = nfilePath
                        break

            if checkExistence:             # we only want to know if the file exists
                return os.path.exists(filePath)

            if not os.path.exists(filePath):
                raise RuntimeError("No such FITS file: " + logLoc.locString())

            if getMetadata:
                from lsst.afw.fits import readMetadata
                return readMetadata(filePath, hdu=0)

            if additionalData.get("llcX"):  # we were asked for a subimage
                bbox = BoxI(PointI(additionalData["llcX"], additionalData["llcY"]),
                            ExtentI(additionalData["width"], additionalData["height"]))
                return afwImage.DecoratedImageF(filePath, hdu=hdu, bbox=bbox)
            else:
                return afwImage.DecoratedImageF(filePath, hdu=hdu)

        if checkExistence:
            return readData(datasetType, butlerLocation, dataId, checkExistence=True)

        md = readData(datasetType, butlerLocation, dataId, getMetadata=True)
        nread = md.get("W_H4NRED")
        hasResetFrame = md["W_4FMTVR"] >= 3 # extra IMG/REF HDUs for the initial reset frame

        hdu_img0 = 2 if hasResetFrame else 0 # hdu for first real read
        #
        # Process the data using CDS
        #
        data = {}
        hdus = 2*np.array([1] if nread == 1 else [1, nread]) - 1
        if "hdu" in dataId:
            hdu = dataId["hdu"]
            read = (hdu + 1)//2         # which read up the ramp do we want?
            if read > nread:
                raise RuntimeError(f"Requested read {read} (HDU {hdu}) is greater than W_H4NRED == {nread}")

            hdus[-1] = hdu

        if hdus[-1]%2 == 0:
            raise RuntimeError(f"HDU {hdus[-1]} is a reference HDU")

        for hdu in hdus:
            read = (hdu + 1)//2         # which read up the ramp do we want?
            try:
                data[hdu] = readData(datasetType, butlerLocation, dataId, hdu=hdu + hdu_img0)
            except FitsError as e:
                self.log.warn("Unable to read IMAGE_%d (hdu %d): %s", read, hdu + hdu_img0, e)
                continue
            assert data[hdu].getMetadata()["EXTNAME"] == f"IMAGE_{read}"

            try:
                ref = readData(datasetType, butlerLocation, dataId, hdu=hdu + hdu_img0 + 1)
            except FitsError as e:
                self.log.warn("Unable to read REF_%d (hdu %d): %s", read, hdu + hdu_img0 + 1, e)
                continue

            assert ref.getMetadata()["EXTNAME"] == f"REF_{read}"

            im = data[hdu].image
            im -= ref.image
            del im

        im = data[hdus[0]].image.array
        #im[im > 30000] = np.NaN

        if hdus[0] == hdus[-1]:         # we can't use CDS
            self.log.warn("You are processing a single frame; not carrying out CDS")
        else:
            im = data[hdus[-1]].image              # data[hdu].image -= data[0].image doesn't work
            im -= data[hdus[0]].image

        return data[hdus[-1]] 

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

        if np.nanmean(exp.variance.array) == 0: # variance isn't set
            gain = exp.getDetector()[0].getGain()

            if exp.getFilterLabel().physicalLabel == 'n':
                gain *= md["W_H4GAIN"]
                var = exp.image.array/gain
                var *= 2                        # CDS
            else:
                var = exp.image.array/gain

            exp.variance.array[:] = var

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

    def std_guider(self, item, dataId):
        """Standardize guider image

        Guider images reside within a PFSD multi-extension FITS file, which
        contains multiple exposures taken with each of six autoguider cameras.
        Each exposure is in its own HDU.
        """
        md = item.getMetadata()
        raw = afwImage.makeExposure(afwImage.makeMaskedImage(item.getImage()))
        raw.setMetadata(md)

        camera = self._makeCamera()
        hdu = dataId["hdu"]  # HDU number; the first data HDU is 1
        guideCam = ((hdu - 1) % 6) + 1

        # The only identifying information we currently have for an HDU is the
        # EXTNAME header; not sure how this will work when we have multiple
        # exposures per guider camera, but hopefully the following will help
        # guard against a mismatch between the data and user expectations.
        extName = md.get("EXTNAME").strip()
        if extName != f"CAM{guideCam}":
            raise RuntimeError(f"Mismatch between assumed camera number ({guideCam}) and header ({extName})")

        detector = camera[f"AG{guideCam}"]
        raw.setDetector(detector)

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

    def _createInitialSkyWcs(self, exposure):
        """Create a SkyWcs from the boresight and camera geometry.

        PFS doesn't have a suitable `SkyWcs`, so this is disabled.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to get data from, and attach the SkyWcs to.
        """
        pass

    def _standardizeExposure(
        self, mapping, item, dataId, filter=True, trimmed=True, setVisitInfo=True, setExposureId=False
    ):
        """Default standardization function for images.

        This sets the Detector from the camera geometry
        and optionally set the Filter. In both cases this saves
        having to persist some data in each exposure (or image).

        This PFS override forces ``filter=False``, since we don't have filters
        and don't use them.

        Parameters
        ----------
        mapping : `lsst.obs.base.Mapping`
            Where to get the values from.
        item : image-like object
            Can be any of lsst.afw.image.Exposure,
            lsst.afw.image.DecoratedImage, lsst.afw.image.Image
            or lsst.afw.image.MaskedImage

        dataId : `dict`
            Dataset identifier
        filter : `bool`
            Set filter? Ignored if item is already an exposure
        trimmed : `bool`
            Should detector be marked as trimmed?
        setVisitInfo : `bool`
            Should Exposure have its VisitInfo filled out from the metadata?
        setExposureId : `bool`
            Should Exposure have its exposure ID filled out from the data ID?

        Returns
        -------
        `lsst.afw.image.Exposure`
            The standardized Exposure.
        """

        if not (isinstance(item, afwImage.ExposureU) or isinstance(item, afwImage.ExposureI) or
                isinstance(item, afwImage.ExposureF) or isinstance(item, afwImage.ExposureD)):
            item = super()._standardizeExposure(
                mapping,
                item,
                dataId,
                filter=False,
                trimmed=trimmed,
                setVisitInfo=setVisitInfo,
                setExposureId=setExposureId
                )

        if item.getFilterLabel() is not None:
            return item

        actualId = mapping.need(["arm"], dataId)
        filterName = actualId["arm"]
        if self.filters is not None and filterName in self.filters:
            filterName = self.filters[filterName]
        item.setFilterLabel(afwImage.FilterLabel(band=filterName, physical=filterName))

        return item
