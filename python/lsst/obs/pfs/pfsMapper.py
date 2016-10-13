#!/usr/bin/env python

import os

from lsst.daf.butlerUtils import CameraMapper
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.policy as pexPolicy

class PfsMapper(CameraMapper):
    """Provides abstract-physical mapping for PFS data"""
    packageName = "obs_pfs"

    def __init__(self, **kwargs):
        policyFile = pexPolicy.DefaultPolicyFile("obs_pfs", "PfsMapper.paf", "policy")
        policy = pexPolicy.Policy(policyFile)
        if not kwargs.get('root', None):
            try:
                kwargs['root'] = os.path.join(os.environ.get('PFS_DATA_DIR'), 'PFS')
            except:
                raise RuntimeError("Either $PFS_DATA_DIR or root= must be specified")
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

    @staticmethod
    def _flipChipsLR(exp, wcs, dataId, dims=None):
        flipLR, flipTB = (True, False)
        exp.setMaskedImage(afwMath.flipImage(exp.getMaskedImage(), flipLR, flipTB))

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
    
    def std_fiberTrace(self, item, dataId):
        return item

    def _extractAmpId(self, dataId):
        ampId = (self._extractDetectorName(dataId), 0, 0)
        return ampId

    def _extractDetectorName(self, dataId):
#        detName = "%(filter)s" % dataId
#        if detName == 'm':
#            detName = 'r'
#        detName = detName + '_' + str("%(spectrograph)s" % dataId)
#        return detName
        return self._extractDetectorId(dataId)

    def _extractDetectorId(self, dataId):
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
        
    def bypass_pfsArm(self, datasetType, pythonType, location, dataId):
        from pfs.datamodel.pfsArm import PfsArm

        for pathName in location.locationList:
            dirName = os.path.dirname(pathName)

            pfsArm = PfsArm(dataId["visit"], dataId["spectrograph"], dataId["arm"])
            try:
                pfsArm.read(dirName=dirName, setPfsConfig=False)
            except Exception as e:
                pass
            else:
                return pfsArm

        raise RuntimeError("Unable to read pfsArm for %s: %s" % (dataId.items(), e))

    def bypass_fiberTrace(self, datasetType, pythonType, location, dataId):
        from pfs.datamodel.pfsFiberTrace import PfsFiberTrace

        for pathName in location.locationList:
            dirName = os.path.dirname(pathName)

            pfsFiberTrace = PfsFiberTrace(dataId["dateObs"], dataId["spectrograph"], dataId["arm"])
            try:
                pfsFiberTrace.read(dirName=dirName)
            except Exception as e:
                pass
            else:
                return pfsFiberTrace

        raise RuntimeError("Unable to read pfsFiberTrace for %s: %s" % (dataId.items(), e))

    @classmethod
    def getCameraName(cls):
        return "pfs"
