#!/usr/bin/env python

import os

import lsst.daf.base as dafBase
from lsst.daf.butlerUtils import CameraMapper
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy

class PfsMapper(CameraMapper):
    """Provides abstract-physical mapping for PFS data"""
    packageName = "obs_subaru"

    def __init__(self, **kwargs):
        policyFile = pexPolicy.DefaultPolicyFile("obs_pfs", "PfsMapper.paf", "policy")
        policy = pexPolicy.Policy(policyFile)
        if not kwargs.get('root', None):
            try:
                kwargs['root'] = os.path.join(os.environ.get('SUPRIME_DATA_DIR'), 'PFS')
            except:
                raise RuntimeError("Either $SUPRIME_DATA_DIR or root= must be specified")
        if not kwargs.get('calibRoot', None):
            kwargs['calibRoot'] = os.path.join(kwargs['root'], 'CALIB')

        super(PfsMapper, self).__init__(policy, policyFile.getRepositoryPath(), **kwargs)

        # Ensure each dataset type of interest knows about the full range of keys available from the registry
        keys = {'field': str,
                'visit': int,
                'ccd': int, # for compatibility with HSC: serial no of ccd
                'spectrograph': int, # [0,1,2,3] for each arm in [blue, red, nir, medred]
                'dateObs': str,
                'taiObs': str,
                'filter': str, # 'arm' called filter for compatibility
                'site': str,
                'category': str,
                }
        for name in ("raw",
                     # processCcd outputs
                     "calexp", 
                     "postISRCCD", 
                     ):
            self.mappings[name].keyDict.update(keys)


        # The order of these defineFilter commands matters as their IDs are used to generate at least some
        # object IDs (e.g. on coadds) and changing the order will invalidate old objIDs

        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter(name="UNRECOGNISED", lambdaEff=0,
                                   alias=["NONE", "None", "Unrecognised", "UNRECOGNISED",
                                          "Unrecognized", "UNRECOGNIZED", "NOTSET",])
        afwImageUtils.defineFilter(name='b', lambdaEff=477, alias=['blue','PFS-B'])
        afwImageUtils.defineFilter(name='r', lambdaEff=623, alias=['red','PFS-R'])
        afwImageUtils.defineFilter(name='n', lambdaEff=623, alias=['nearInfraRed','PFS-N'])
        afwImageUtils.defineFilter(name='m', lambdaEff=775, alias=['mediumResolutionRed','PFS-M'])
        #
        # self.filters is used elsewhere, and for now we'll set it
        #
        # It's a bit hard to initialise self.filters properly until #2113 is resolved,
        # including the part that makes it possible to get all aliases
        #
        self.filters = {}
        for f in [
            "b",
            "r",
            "n",
            "m",
            "NONE",
            "UNRECOGNISED"]:
            # Get the canonical name -- see #2113
            self.filters[f] = afwImage.Filter(afwImage.Filter(f).getId()).getName()
        self.defaultFilterName = "UNRECOGNISED"

        #
        # The number of bits allocated for fields in object IDs, appropriate for
        # the default-configured Rings skymap.
        #
        # This shouldn't be the mapper's job at all; see #2797.

        PfsMapper._nbit_id = 64# - (PfsMapper._nbit_tract + 2*PfsMapper._nbit_patch + PfsMapper._nbit_filter)

        #if len(afwImage.Filter.getNames()) >= 2**PfsMapper._nbit_filter:
        #    raise RuntimeError("You have more filters defined than fit into the %d bits allocated" %
        #                       PfsMapper._nbit_filter)

    def map(self, datasetType, dataId, write=False):
        """Need to strip 'flags' argument from map

        We want the 'flags' argument passed to the butler to work (it's
        used to change how the reading/writing is done), but want it
        removed from the mapper (because it doesn't correspond to a
        registry column).
        """
        copyId = dataId.copy()
        copyId.pop("flags", None)
        location = super(PfsMapper, self).map(datasetType, copyId, write=write)

        if 'flags' in dataId:
            location.getAdditionalData().set('flags', dataId['flags'])

        return location

    @staticmethod
    def _flipChipsLR(exp, wcs, dataId, dims=None):
        """Flip the chip left/right or top/bottom. Process either/and the pixels and wcs
           Most chips are flipped L/R, but the rotated ones (100..103) are flipped T/B
        """
        flipLR, flipTB = (False, True) if dataId['ccd'] in (100, 101, 102, 103) else (True, False)
        if exp:
            exp.setMaskedImage(afwMath.flipImage(exp.getMaskedImage(), flipLR, flipTB))
        if wcs:
            wcs.flipImage(flipLR, flipTB, exp.getDimensions() if dims is None else dims)

        return exp

    def std_raw(self, item, dataId):
        exp = super(PfsMapper, self).std_raw(item, dataId)

        md = exp.getMetadata()
        if md.exists("MJD-STR"):
            calib = exp.getCalib()
            expTime = calib.getExptime()
            obsStart = dafBase.DateTime(md.get("MJD-STR"), dafBase.DateTime.MJD, dafBase.DateTime.UTC)
            obsMidpoint = obsStart.nsecs() + long(expTime * 1000000000L / 2)
            calib.setMidTime(dafBase.DateTime(obsMidpoint))

        print "pfsMapper::std_raw: exp.getCalib().getExptime() = ", exp.getCalib().getExptime();

        return exp
    
    def standardizeCalib(self, dataset, item, dataId):
        """Standardize a calibration image read in by the butler

        Some calibrations are stored on disk as Images instead of MaskedImages
        or Exposures.  Here, we convert it to an Exposure.

        @param dataset  Dataset type (e.g., "bias", "dark" or "flat")
        @param item  The item read by the butler
        @param dataId  The data identifier (unused, included for future flexibility)
        @return standardized Exposure
        """
        mapping = self.calibrations[dataset]
        if "MaskedImage" in mapping.python:
            exp = afwImage.makeExposure(item)
        elif "Image" in mapping.python:
            if hasattr(item, "getImage"): # For DecoratedImageX
                item = item.getImage()
            exp = afwImage.makeExposure(afwImage.makeMaskedImage(item))
        elif "Exposure" in mapping.python:
            exp = item
        else:
            raise RuntimeError("Unrecognised python type: %s" % mapping.python)

        parent = super(PfsMapper, self)
        if hasattr(parent, "std_" + dataset):
            return getattr(parent, "std_" + dataset)(exp, dataId)
        return self._standardizeExposure(mapping, exp, dataId)

    def std_bias(self, item, dataId):
        return self.standardizeCalib("bias", item, dataId)

    def std_dark(self, item, dataId):
        exp = self.standardizeCalib("dark", item, dataId)
        exp.getCalib().setExptime(1.0)
        return exp

    def std_flat(self, item, dataId):
        return self.standardizeCalib("flat", item, dataId)

    def _extractAmpId(self, dataId):
        ampId = (self._extractDetectorName(dataId), 0, 0)
#        print 'PfsMapper._extractAmpId: ampId = ',ampId
        return ampId

    def _extractDetectorName(self, dataId):
#        print 'PfsMapper._extractDetectorName: dataId = ',dataId
        detName = "%(filter)s" % dataId
        if detName == 'm':
            detName = 'r'
#        detName = detName + '_' + str("%(spectrograph)s" % dataId)
        detId = int("%(ccd)d" % dataId)
#        print 'PfsMapper._extractDetectorName = <',detName,'>'
        print 'PfsMapper._extractDetectorName = <',detId,'>'
#        return detName
        return detId

    def _extractDetectorId(self, dataId):
#        print 'PfsMapper._extractDetectorId: dataId = ',dataId
        detId = int("%(ccd)d" % dataId)
        return detId
    
    def _getCcdKeyVal(self, dataId):
        """Return CCD key and value used to look a defect in the defect registry

        The default implementation simply returns ("ccd", full detector name)
        """
        return ("ccd", self._extractDetectorId(dataId))
#        return ("ccd", self._extractDetectorName(dataId))

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        @param dataId (dict) Data identifier with visit, ccd
        """
        pathId = self._transformId(dataId)
        visit = pathId['visit']
        ccd = pathId['ccd']
        return visit*200 + ccd

    def bypass_ccdExposureId(self, datasetType, pythonType, location, dataId):
        return self._computeCcdExposureId(dataId)

    def bypass_ccdExposureId_bits(self, datasetType, pythonType, location, dataId):
        """How many bits are required for the maximum exposure ID"""
        return 32 # just a guess, but this leaves plenty of space for sources

    def _computeCoaddExposureId(self, dataId, singleFilter):
        """Compute the 64-bit (long) identifier for a coadd.

        @param dataId (dict)       Data identifier with tract and patch.
        @param singleFilter (bool) True means the desired ID is for a single- 
                                   filter coadd, in which case dataId
                                   must contain filter.
        """

        tract = long(dataId['tract'])
        if tract < 0 or tract >= 2**PfsMapper._nbit_tract:
            raise RuntimeError('tract not in range [0,%d)' % (2**PfsMapper._nbit_tract))
        patchX, patchY = map(int, dataId['patch'].split(','))
        for p in (patchX, patchY):

            if p < 0 or p >= 2**PfsMapper._nbit_patch:
                raise RuntimeError('patch component not in range [0, %d)' % 2**PfsMapper._nbit_patch)
        oid = (((tract << PfsMapper._nbit_patch) + patchX) << PfsMapper._nbit_patch) + patchY
        if singleFilter:
            return (oid << PfsMapper._nbit_filter) + afwImage.Filter(dataId['filter']).getId()
        return oid

    def std_psf(self, calexp, dataId):
        return calexp.getPsf()

    @classmethod
    def getCameraName(cls):
        return "pfs"
