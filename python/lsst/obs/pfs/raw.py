from typing import Literal, Optional, TYPE_CHECKING, Union

import numpy as np

import lsst.afw.math as afwMath
import lsst.afw.image as afwImage

from astro_metadata_translator import fix_header, ObservationInfo
from lsst.geom import Box2I, Point2I, Extent2I
from lsst.afw.fits import readMetadata
from lsst.afw.image import DecoratedImageF, ImageU, ImageF, VisitInfo
from lsst.obs.base.makeRawVisitInfoViaObsInfo import MakeRawVisitInfoViaObsInfo

from .loadCamera import loadCamera
from .translator import PfsTranslator

if TYPE_CHECKING:
    from lsst.afw.cameraGeom import Detector
    from lsst.daf.base import PropertyList

__all__ = ("PfsRaw",)


class PfsRaw:
    """The raw data from a PFS detector

    PFS raw data comes in two flavors: CCD (for the b,r,m arms) and NIR
    (for the n arm). Dealing with the CCD data is fairly straightforward,
    but the NIR data is more complex (involving multiple reads). This class
    provides a common interface for both, providing an Image to ISR.

    It is essentially a representation of a data cube (the CCD data has a single
    plane; the NIR data has multiple planes, stored as a series of FITS HDUs
    rather than an actual cube, but that's a technicality). It is intended to be
    used as a helper to the Formatter.

    Parameters
    ----------
    path : `str`
        Path to the raw data file.
    """

    def __init__(self, path: str, pfsCategory: Union[None, Literal["F"], Literal["L"]]) -> None:
        """Construct from a filename"""
        self.path = path
        self.pfsCategory = pfsCategory
        self._metadata: Optional[PropertyList] = None
        self._detector: Optional["Detector"] = None
        self._obsInfo: Optional[ObservationInfo] = None
        self._visitInfo: Optional[VisitInfo] = None

        self._h4Config: Optional[Dict] = None

    @property
    def metadata(self) -> "PropertyList":
        """Return the metadata"""
        if self._metadata is None:
            self._metadata = readMetadata(self.path)
            fix_header(self._metadata, translator_class=PfsTranslator, filename=self.path)

        return self._metadata

    @property
    def detector(self) -> "Detector":
        """Return the detector"""
        detId = self.metadata["DET-ID"]
        detector = loadCamera(self.pfsCategory)[detId]
        if self.isNir():
            detector = detector.rebuild()   # returns a Detector Builder
            gain = self.metadata["W_H4GAIN"]  # ASIC pre-amp gain (electrons per ADU)

            # This will need work for non-32-channel ramps.
            for channel in detector.getAmplifiers():  # an H4RG/ROIC/SAM "channel", not really an amplifier
                channel.setGain(channel.getGain() / gain)
            detector = detector.finish()
        return detector

    @property
    def obsInfo(self) -> ObservationInfo:
        """Return the ObservationInfo"""
        if self._obsInfo is None:
            self._obsInfo = ObservationInfo(
                self.metadata, translator_class=PfsTranslator, filename=self.path
            )
        return self._obsInfo

    @property
    def visitInfo(self) -> VisitInfo:
        """Return the ExposureInfo"""
        if self._visitInfo is None:
            self._visitInfo = MakeRawVisitInfoViaObsInfo.observationInfo2visitInfo(self.obsInfo)
        return self._visitInfo

    @property
    def bbox(self) -> Box2I:
        """Return the bounding box"""
        return Box2I(self.xy0, self.dimensions)

    @property
    def dimensions(self) -> Extent2I:
        """Return the image dimensions"""
        if self.metadata.get("ZIMAGE"):
            return Extent2I(self.metadata["ZNAXIS1"], self.metadata["ZNAXIS2"])
        return Extent2I(self.metadata["NAXIS1"], self.metadata["NAXIS2"])

    @property
    def xy0(self) -> Point2I:
        """Return the origin of the image"""
        return Point2I(0, 0)

    def isNir(self) -> bool:
        """Return if this is NIR data"""
        return self.metadata.get("W_ARM") == 3

    def isCcd(self) -> bool:
        """Return if this is CCD data"""
        return not self.isNir()

    def getImage(self, bbox: Optional["Box2I"] = None) -> Union[ImageU, ImageF]:
        """Return the image

        For CCD data, this is the raw image. For NIR data, this is the "CDS"
        image (last read minus the first read).

        Parameters
        ----------
        bbox : `Box2I`
            The bounding box to return.

        Returns
        -------
        image : `ImageU` or `ImageF`
            The image.
        """
        return self.getNirImage(bbox) if self.isNir() else self.getCcdImage(bbox)

    def getCcdImage(self, bbox: Optional["Box2I"] = None) -> ImageU:
        """Return the CCD image

        Parameters
        ----------
        bbox : `Box2I`
            The bounding box to return.

        Returns
        -------
        image : `ImageU`
            The CCD image.
        """
        if not self.isCcd():
            raise RuntimeError("This is not a CCD image")

        args = [bbox] if bbox is not None else []
        image = ImageU(self.path, *args)

        try:
            dataVersion = int(self.metadata.get("W_VERSIONS_FPGA"), 16)
        except Exception:
            dataVersion = None
        if dataVersion is not None and dataVersion <= 0x0070:
            array = image.getArray()
            for amp in self.detector:
                ySlice, xSlice = amp.getRawBBox().getSlices()
                ampIm = array[ySlice, xSlice]
                ampShape = ampIm.shape

                # Don't bother being tricky: make a copy.
                ampPixels = ampIm.flatten()
                ampPixels[:-1] = ampPixels[1:]

                array[ySlice, xSlice] = ampPixels.reshape(ampShape)

        return image

    def getNumReads(self) -> int:
        """Return the number of reads"""
        return self.metadata.get("W_H4NRED") if self.isNir() else 1

    def getNirImage(self, bbox: Optional["Box2I"] = None) -> ImageF:
        """Return a quicklook estimate of the NIR image

        This is the "CDS" image (the last read minus the first read).

        Parameters
        ----------
        bbox : `Box2I`
            The bounding box to return.

        Returns
        -------
        image : `ImageF`
            The NIR image.
        """
        numReads = self.getNumReads()

        first = self.getCorrectedNirRead(1, bbox=bbox)
        if numReads == 1:
            # Special case: can't really do "CDS"
            return first

        last = self.getCorrectedNirRead(numReads, bbox=bbox)
        last -= first
        return last

    def getRawNirRead(self, hdu):
        """The single method which reads from the ramp file.
        
        Only need to call this directly to read reset frames or other oddballs.
        
        Args
        ----
        hdu : `int` or `str`
           If a string, the EXTNAME. See `.getImage()` and `.getRef()`
           If an int, the 1-based HDU to read, where currently the order is:
            1 - reset image
            2 - reset ref image
            3 + 2*readNum - the 0-indexed data image
            4 + 2*readNum - the 0-indexed ref image
            And all that might change, in the face of multiple resets or ramps.
            
        Returns
        -------
        image : `afw.image.maskedImage`
           The image from the given HDU. Seems to be float by default
        """

        # Check whether  is actually efficient -- CPL
        # Hmm, inefficient for reading stacks, by a factor of ~2. open/close overhead?
        img = self.pfsRaw.get(hdu=hdu)
        if False and np.isnan(img.image.array).any():
            self.logger.warning(f'NaNs in {self.dataId} {hdu}: {np.isnan(img.image.array).sum()}')
        return img.maskedImage

    def positiveIndex(self, n: int) -> int:
        """Convert possibly negative index to always positive index."""
        nums = range(0, self.getNumReads())
        return nums[n]
    
    def _readIdxToPFSBIdx(self, n: int) -> int:
        """Convert possibly negative 0-indexed ramp read index into positive 1-indexed PFSB read index"""
        return self.positiveIndex(n) + 1

    def getRawDataImage(
        self,
        readNum: int,
    ) -> ImageF:
        """Get a data image.

        Args
        ----
        readNum : `int`
           the 0-indexed data image to read.

        Returns
        -------
        image : `afw.image.maskedImage`
              'u2' by default.
        """

        readNum = self._readIdxToPFSBIdx(readNum)
        return self.getNirRead(readNum, imageType='IMAGE')

    def getRawIrpImage(
        self, 
        readNum: int,
    ) -> ImageF:
        """Get an IRP ("reference") image. 
        
        Args
        ----
        readNum : `int`
           the 0-indexed IRP image to read.
           
        Returns
        -------
        image : `afw.image.maskedImage`
              'u2' by default.
        """
        readNum = self._readIdxToPFSBIdx(readNum)
        return self.getNirRead(readNum, imageType='REF')

    @staticmethod
    def evilGetArray(image: ImageF) -> np.ndarray:
        """Get the array from a maskedImage, and REPLACE NaNs with 0 (the FITS BLANK value)

        This is horrendous. Subaru requires us to have BLANK cards, and the Rubin reader replaces
        matching pixels with NaNs. BLANKs are supposed to indicate undefined pixels, but 0 and 65535 
        are both valid pixel values from the H4s.

        When we know we do not care and need numpy arrays (e.g. when building stacks)
        call this to replace NaNs with 0. 
        """
        
        arr = image.getArray()
        nanmask = np.isnan(arr)
        if nanmask.any():
            # I feel so dirty -- CPL
            arr[nanmask] = 0.0
        return arr

    def getNirRead(
        self,
        readNum: int,
        imageType: Union[Literal["IMAGE"], Literal["REF"]] = "IMAGE",
        bbox: Optional["Box2I"] = None,
    ) -> ImageF:
        """Return the image for the given read

        Parameters
        ----------
        readNum : `int`
            The read number to return.
        imageType : `str`
            The type of image to return: ``IMAGE`` or ``REF``.
        bbox : `Box2I`
            The bounding box to return.

        Returns
        -------
        image : `ImageF`
            The image for the given read.
        """
        numReads = self.getNumReads()
        if not (1 <= readNum <= self.getNumReads()):
            raise RuntimeError(f"The read must be in the range 1..{numReads}; saw {readNum}")
        hdu = 2*readNum - 1 + dict(IMAGE=0, REF=1)[imageType]

        hasResetFrame = self.metadata["W_4FMTVR"] >= 3  # extra IMG/REF HDUs for the initial reset frame
        if hasResetFrame:
            hdu += 2

        args = [bbox] if bbox is not None else []
        image = DecoratedImageF(self.path, hdu, *args)

        extname = image.getMetadata()["EXTNAME"]
        expectExtname = f"{imageType}_{readNum}"
        if extname != expectExtname:
            raise RuntimeError(f"Expected to see EXTNAME = {expectExtname}; saw {extname}")

        # Rotate immediately and always: completely hide the fact that the detector is physically rotated.
        return afwMath.rotateImageBy90(image.getImage(), 3)

    def _simpleStack(
        self, 
        reader=None, 
        r0: int=0, 
        r1: int=-1, 
        step: int=1, 
        bbox: Optional["Box2I"]=None, 
    ) -> np.ndarray:
        """Return all the raw frames in a single 3d stack.

        Args
        ----
        reader : `callable`
          A function to read the images. If None, uses `self.getRawDataImage`.
        r0 : `int`
          The 0-indexed read to start from.
        r1 : `int`
          The 0-indexed read to end with. Inclusive.
        step : `int`
          A `range`-compliant stepping between r0 and r1.
        bbox

        Returns
        -------
        stack : 3-d float32 numpy array
           the stack, with axis 0 being the reads. Always 'f4'. [ Not sure how to encode to 
           3-d maskedImages. ]
        """

        if reader is None:
            reader = self.getRawDataImage
            
        r0 = self.positiveIndex(r0)
        r1 = self.positiveIndex(r1)
        if r1 <= r0:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}):')
        reads = np.arange(r0, r1, step)
        nreads = len(reads)

        stack = None
        for r_i in reads:
            readImg = reader(r_i)
            if stack is None:
                h, w = readImg.array.shape
                stack = np.empty(shape=(nreads, h, w), dtype='f4')
            stack[r_i,:,:] = self.evilGetArray(readImg)

        return stack

    def rawDataStack(self, r0=0, r1=-1, step=1) -> np.ndarray:
        """Return all the raw image frames in a single 3d stack."""
        return self._simpleStack(reader=self.getRawDataImage, 
                                 r0=r0, r1=r1, step=step)
    def rawIrpStack(self, r0=0, r1=-1, step=1) -> np.ndarray:
        """Return all the raw IRP frames in a single 3d stack."""
        return self._simpleStack(reader=self.getRawIrpImage, 
                                 r0=r0, r1=r1, step=step)

    def interpolateChannelIrp(self, 
                              rawChan: np.ndarray, 
                              doFlip: bool,
    ) -> np.ndarray:
        """Given an IRPn channel from a PFSB file, return IRP1-sized channel image.

        Args
        ----
        rawChan : array
           The raw IRP channel, with the columns possibly subsampled by a factor of self.irpN.
        doFlip : `bool`
           Whether we need to treat this channel as read out R-to-L. Only meaningful if
           we are filtering in time.
           
        Returns
        -------
        im : numpy array
           the full-sized interpolated reference pixel channel.

        We do not yet know how best to interpolate, so simply repeat the pixels refRatio times.
        This needs to be integrated with handling the bad column mask.
        """

        if self.irpN == 0:
            return rawChan

        # Allow the future possibility of using readout order.
        temporalFilter = False
                
        irpHeight, irpWidth = rawChan.shape
        refChan = np.empty(shape=(irpHeight, irpWidth * self.irpN), dtype=rawChan.dtype)

        if doFlip and temporalFilter:
            rawChan = rawChan[:, ::-1]

        # Simply repeat reference pixels
        for i in range(0, self.irpN):
            refChan[:, i::self.irpN] = rawChan

        if doFlip and temporalFilter:
            refChan = refChan[:, ::-1]

        return refChan

    def constructFullIrp(self, refImg):
        """Given an IRPn image, return IRP1 image.

        Args
        ----
        rawImg : ndarray
          A raw reference read from the ASIC.

        Returns
        -------
        img : `afw.ImageF`
          full 4096x4096 image.

        - The detector was read out in nChannel channels, usually 32, but possibly 16, 4, or 1.

        - If oddEven is set (the default), the read order from the
          detector of pairs of channels is --> <--. The ASIC "corrects"
          that order so that the image always "looks" right: the
          column-order of the image from the ASIC is spatially, not temporally, correct.

        - the N:1 ratio of science to reference pixels is deduced from the size of the image.

        - refPix tells us the position of the reference pixel within the N
          science pixels. It must be >= 1 (there must be at least one
          science pixel before the reference pixel). The ASIC default is
          for it to be the last pixel in the group, but we usually try to 
          put the reference pixel in the middle of the block of science pixels.

        """

        rawIrp = self.evilGetArray(refImg)
        height, width = rawIrp.shape

        # If we are a full frame, no interpolation is necessary.
        if width == self.h4Size:
            return refImg

        dataChanWidth = self.h4Size // self.nchan
        refChanWidth = width // self.nchan
        refRatio = self.h4Size // width

        refChans = []
        for c_i in range(self.nchan):
            rawChan = rawIrp[:, c_i*refChanWidth:(c_i+1)*refChanWidth]
            doFlip = self.oddEven and c_i%2 == 1

            # This is where we would intelligently interpolate.
            refChan = self.interpolateChannelIrp(rawChan, doFlip)
            refChans.append(refChan)

        fullRefImg = np.hstack(refChans)

        return afwImage.makeImageFromArray(fullRefImg)

    def getCorrectedNirRead(
        self,
        readNum: int,
        bbox: Optional["Box2I"] = None,
    ) -> ImageF:
        """Return the corrected image for the given read

        The corrected image is the "IMAGE" minus (a statistical representation
        of) the "REF".

        Parameters
        ----------
        readNum : `int`
            The read number to return.
        bbox : `Box2I`
            The bounding box to return. Note that even if the bbox is specified,
            we will read and process the entire image before returning just the
            desired subimage.

        Returns
        -------
        image : `ImageF`
            The corrected image for the given read.
        """
        image = self.getNirRead(readNum, imageType="IMAGE")
        ref = self.getNirRead(readNum, imageType="REF")

        isBadVersion = self.metadata["W_SRH4"] in (18660, "no available value")

        # Do the reference subtraction
        nChannel = 32
        channelWidth = 128
        start = 0
        for ii in range(nChannel):
            stop = start + channelWidth
            imageChannel = image.array[:, start:stop]
            refChannel = ref.array[:, start:stop]
            start = stop

            refPerColumn = np.nanmedian(refChannel, axis=0)
            refChannel -= refPerColumn

            if ii == 26 and isBadVersion:
                refPerColumn = refPerColumn[96:]  # exclude the weird bad column

            refCorrection = np.nanmedian(refChannel, axis=1) + np.mean(refPerColumn)
            imageChannel -= refCorrection.reshape((len(refCorrection), 1))

        if isBadVersion:
            xc, yc = 3280, 3545
            Rx, Ry = 40, 35
            image.array[yc - Ry:yc + Ry + 1, xc - Rx:xc + Rx + 1] = np.NaN

        if bbox is not None:
            image = image[bbox]

        return image
