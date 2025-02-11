from typing import Literal, Optional, TYPE_CHECKING, Union, Tuple

import numpy as np

import lsst.afw.math as afwMath

from astro_metadata_translator import fix_header, ObservationInfo
from lsst.geom import Box2I, Point2I, Extent2I
from lsst.afw.fits import readMetadata
from lsst.afw.image import DecoratedImageF, ImageU, ImageF, VisitInfo
from lsst.afw.image import makeExposure, makeMaskedImage, ExposureF, FilterLabel
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

    def __init__(self, path: str, pfsCategory: Union[None, Literal["F"], Literal["L"], Literal["J"]]) -> None:
        """Construct from a filename"""
        self.path = path
        self.pfsCategory = pfsCategory
        self._metadata: Optional[PropertyList] = None
        self._detector: Optional["Detector"] = None
        self._obsInfo: Optional[ObservationInfo] = None
        self._visitInfo: Optional[VisitInfo] = None

    @property
    def metadata(self) -> "PropertyList":
        """Return the metadata"""
        if self._metadata is None:
            self._metadata = readMetadata(self.path)
            fix_header(self._metadata, translator_class=PfsTranslator, filename=self.path)

        return self._metadata

    def getMetadata(self) -> "PropertyList":
        """Hoist the .exposure.getMetadata() up"""
        return self.metadata

    @property
    def detector(self) -> "Detector":
        """Return the detector"""
        detId = self.metadata["DET-ID"]
        detector = loadCamera(self.pfsCategory)[detId]
        if self.isNir():
            detector = detector.rebuild()  # returns a Detector Builder
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
            self._obsInfo = ObservationInfo(self.metadata, translator_class=PfsTranslator, filename=self.path)
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

    @property
    def nchan(self) -> int:
        """Return the number of H4 readout channels"""
        return self.metadata["W_H4NCHN"]

    @property
    def irpN(self) -> int:
        """Return the ratio of data:IRP pixels."""
        return self.metadata["W_H4IRPN"]

    @property
    def irpOffset(self) -> int:
        """Return the offset of the IRP pixel in the block of irpN data pixels."""
        return self.metadata["W_H4IRPO"]

    @property
    def h4Size(self) -> int:
        return 4096

    def getH4channelReadOrder(self) -> Tuple[bool, bool]:
        """Return the temporal read directions for alternating channels

        Returns
        -------
        evenChannels, oddChannels : `int`
           if 1, channel is read out right-to-left.
        """

        # Mumble. The readout directions are not broken out in the H4 header.
        """
        SPI Register 13 contains:
          HINVDIR field: bit 8 for even channels, bit 9 for odd.
          with 0 = left to right, 1 = right to left
        """
        hinvdir = (self.metadata["W_4SPI13"] & 0x300) >> 8
        return hinvdir & 1, (hinvdir >> 1) & 1

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

    def getExposure(self, bbox: Optional["Box2I"] = None) -> ExposureF:
        """Return an Exposure

        Wraps the image in an Exposure, setting the detector, visitInfo, and
        metadata.

        Parameters
        ----------
        bbox : `Box2I`
            The bounding box to return.

        Returns
        -------
        exposure : `ExposureF`
            The Exposure.
        """
        image = self.getImage(bbox).convertF()
        exposure = makeExposure(makeMaskedImage(image))
        exposure.setDetector(self.detector)
        info = exposure.getInfo()
        info.setVisitInfo(self.visitInfo)
        info.setId(self.visitInfo.id)
        info.setMetadata(self.metadata)
        info.setDetector(self.detector)
        arm = self.obsInfo.ext_arm
        info.setFilter(FilterLabel(arm, arm))
        return exposure

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
           the 0-indexed data image to read. If negative,
           take from the end of the ramp.

        Returns
        -------
        image : `afw.image.maskedImage`
              'u2' by default.
        """

        readNum = self._readIdxToPFSBIdx(readNum)
        return self.getNirRead(readNum, imageType="IMAGE")

    def getRawIrpImage(
        self,
        readNum: int,
    ) -> ImageF:
        """Get an IRP ("reference") image.

        Args
        ----
        readNum : `int`
           the 0-indexed IRP image to read. If negative,
           take from the end of the ramp.

        Returns
        -------
        image : `afw.image.maskedImage`
              'u2' by default.
        """
        readNum = self._readIdxToPFSBIdx(readNum)
        return self.getNirRead(readNum, imageType="REF")

    @staticmethod
    def replaceNansWith0(image: ImageF) -> ImageF:
        """REPLACE NaNs in an ImageF with 0 (the FITS BLANK value)

        This is horrendous. Subaru requires us to have BLANK cards, and the Rubin reader replaces
        matching pixels with NaNs. BLANKs are supposed to indicate undefined pixels, but 0 and 65535
        are both valid pixel values from the H4s: there is no "undefined" value.

        When we know we do not care and need numpy arrays (e.g. when building stacks)
        call this to replace NaNs with 0.
        """

        arr = image.getImage().getArray()
        nanmask = np.isnan(arr)
        if nanmask.any():
            arr[nanmask] = 0.0
        return arr

    def getNirRead(
        self,
        readNum: int,
        imageType: Union[Literal["IMAGE"], Literal["REF"]] = "IMAGE",
        bbox: Optional["Box2I"] = None,
        doRotate: bool = True,
    ) -> ImageF:
        """Return the image for the given read

        Parameters
        ----------
        readNum : `int`
            The read number to return. 1-indexed.
        imageType : `str`
            The type of image to return: ``IMAGE`` or ``REF``.
        bbox : `Box2I`
            The bounding box to return.

        Returns
        -------
        image : `ImageF`
            The image for the given read.
            The image is rotated so that traces run vertically.
        """
        numReads = self.getNumReads()
        if not (1 <= readNum <= self.getNumReads()):
            raise RuntimeError(f"The read must be in the range 1..{numReads}; saw {readNum}")
        hdu = 2 * readNum - 1 + dict(IMAGE=0, REF=1)[imageType]

        hasResetFrame = self.metadata["W_4FMTVR"] >= 3  # extra IMG/REF HDUs for the initial reset frame
        if hasResetFrame:
            hdu += 2

        args = [bbox] if bbox is not None else []
        image = DecoratedImageF(self.path, hdu, *args)

        extname = image.getMetadata()["EXTNAME"]
        expectExtname = f"{imageType}_{readNum}"
        if extname != expectExtname:
            raise RuntimeError(f"Expected to see EXTNAME = {expectExtname}; saw {extname}")

        # The Subaru-required BLANK card does not work for H4s: we have to replace NaNs with 0.
        self.replaceNansWith0(image)

        if doRotate:
            # Rotate immediately and always: completely hide the fact that the detector is physically rotated.
            image = afwMath.rotateImageBy90(image.getImage(), self.detector.getOrientation().getNQuarter())
        else:
            image = image.getImage()
        return image

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
