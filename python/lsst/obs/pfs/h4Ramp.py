import logging
import time

import numpy as np
import scipy

import lsst.afw.image as afwImage

class H4Ramp:
    """A PFSB H4 ramp.
    
    See datamodel/PFSB.md for details.
    """
    h4Size = 4096   # the total size of the full detector
    refPix = 4      # the number of reference border pixels within the full detector. 
                    # According to the H4 docs these can be replaced by live pixels, 
                    # but I don't see us doing that.

    # Get these into some config file(s) -- CPL
    # The irpDefects are generated using the calcIrpDefects method, and are not yet curated. Badly need to be.
    corrections = {18734:dict(leakage=np.ones(32, dtype='f4')*0.012, 
                              irpNoiseClip=5.0),
                   18660:dict(leakage=np.ones(32, dtype='f4')*0.003,
                              irpNoiseClip=5.0,
                              irpDefects0=np.array([ 196,  303,  334,  348,  821,  901,  996, 1282,
                                                   1449, 1469, 1473, 1551, 1582, 1645, 1723, 1831, 
                                                   1945, 2049, 2200, 2448, 2600, 2730, 2865, 3036, 
                                                   3084, 3118, 3169, 3343, 3692, 3753, 3865]),
                              irpDefects=np.array([196,  334,  821, 1449, 1551, 1582, 1645, 1723, 
                                                   2600, 3084, 3118, 3169, 3726, 3794, 3865])),
                   18315:dict(leakage=[0.0048, 0.00486, -0.0, 0.0048, 0.00475, 0.00506, 0.00518, 0.00478, 
                                       0.00479, 0.00477, 0.00493, 0.0047, 0.00456, 0.00467, 0.00453, 0.0044, 
                                       0.00427, 0.00466, 0.00462, 0.00467, 0.00459, 0.00467, 0.0048, 0.00535, 
                                       0.0054, 0.00528, 0.00563, 0.00547, 0.0058, 0.00587, 0.0062, 0.00513], 
                              irpNoiseClip=5.0,
                              irpDefects0=np.array([ 73, 518, 556, 1120, 1137, 1242, 1565, 1607, 1947, 2067, 
                                                   2113, 2194, 2240, 2294, 2460, 2676, 2705, 2858, 
                                                   3238, 3249, 3276, 3327, 3332, 3455, 3479, 3492, 
                                                   3590, 3953, 3961, 3979, 4083]),
                              irpDefects=np.array([471,  518, 1565, 1607, 2676, 2858, 3238, 3455, 
                                                   4060, 4064, 4065, 4066])),
                   18321:dict(leakage=[0.00214, 0.00255, 0.00231, 0.00224, 0.00205, 0.00243, 0.00229, 0.00225, 
                                       0.0022, 0.00221, 0.00224, 0.00187, 0.0022, 0.00193, 0.00177, 0.00111, 
                                       0.00172, 0.00214, 0.00167, 0.00192, 0.00208, 0.00233, 0.00224, 0.00245, 
                                       0.00264, 0.00324, 0.00316, 0.00286, 0.00308, 0.003, 0.00278, 0.00265], 
                              irpNoiseClip=5.0,
                              irpDefects0=np.array([  51,   99,  284,  367,  770, 1047, 1119, 1368, 
                                                   1462, 1483, 1649, 1748, 2061, 2120, 2216, 2608, 
                                                   2686, 2694, 2829, 2837, 2941, 3040, 3046, 3332, 3618]),
                              irpDefects=np.array([284,  330,  367, 1649, 1748, 2216, 2829])),
                   18661:dict(leakage=np.ones(32, dtype='f4')*0.005,
                              irpNoiseClip=5.0,
                              irpDefects=np.array([959, 1219, 1654], dtype='i4')),
    }

    def __init__(self, pfsRaw, disableLeakage=True):
        self.logger = logging.getLogger('H4Ramp')
        self.pfsRaw = pfsRaw
        
        self.md = md = pfsRaw.metadata
        self.nread = md['W_H4NRED']
        self.nchan = md['W_H4NCHN'] # For single-asic instruments like PFS: 32, 16, 4, or 1
        self.irpN = md['W_H4IRPN']  # The number of science pixels per reference pixel.
        self.irpOffset = md['W_H4IRPO'] # The offset of the reference pixel within the N science pixels.

        # Do the channel readout directions flip?
        #   True:  --> <-- 
        #   False: --> -->
        # We *very* rarely use False, and did not encode the other two options. Should have.
        # Mumble. not broken out in the H4 header -- CPL
        readOrder = (md['W_4SPI13'] & 0x300) >> 8
        self.oddEven = readOrder == 0b10
        
        # IRP columns we need to interpolate over because of defects in the IRP row.
        self.irpDefects = self.corrections[self.h4Serial]['irpDefects']

        # How much of the data to subtract from the IRP. 
        # More complicated than just a simple factor, sad to say. So disable for now.
        self.dataToIrpLeakageFactor = self.corrections[self.h4Serial]['leakage']
        if False and disableLeakage:
            self.dataToIrpLeakageFactor = 0.0

        # wire in default routines
        self.borderCorrect = self.refpix4BorderCorrect
        
    def __str__(self):
        return (f'H4Ramp(nread={self.nread}, IRP=({self.irpN},{self.irpOffset})')

    def __repr__(self):
        return str(self)
        
    @classmethod
    def loadFromButler(cls, butler, dataId, disableLeakage=True, **kw):
        dataId = dataId.copy()
        dataId.update(kw)

        raw = butler.get('raw', dataId)
        return cls(raw, disableLeakage=disableLeakage)

    @property
    def haveIrp(self):
        return self.irpN > 0

    @property
    def h4Serial(self):
        return self.md['W_SRH4']
        
    def getRawRead(self, hdu):
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

        # Check whether butler.get('rawb') is actually efficient -- CPL
        # Hmm, inefficient for reading stacks, by a factor of ~2. open/close overhead?
        img = self.pfsRaw.get(hdu=hdu)
        if False and np.isnan(img.image.array).any():
            self.logger.warning(f'NaNs in {self.dataId} {hdu}: {np.isnan(img.image.array).sum()}')
        return img.maskedImage
    
    def positiveIndex(self, n):
        """Convert possibly negative 0-indexed ramp read index into positive 0-indexed read"""
        nums = range(0, self.nread)
        return nums[n]

    def _readIdxToFITSIdx(self, n):
        """Convert possibly negative 0-indexed ramp read index into positive 1-indexed HDU"""
        return self.positiveIndex(n) + 1

    def getRawDataImage(self, readNum):
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

        readNum = self._readIdxToFITSIdx(readNum)
        return self.pfsRaw.getNirRead(readNum, imageType='IMAGE')

    def getRawIrpImage(self, readNum):
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
        readNum = self._readIdxToFITSIdx(readNum)
        return self.pfsRaw.getNirRead(readNum, imageType='REF')

    @staticmethod
    def evilGetArray(image):
        """Get the array from a maskedImage, and REPLACE NaNs with 0 (the FITS BLANK value)

        This is horrendous. We are required to have BLANK cards, and the Rubin reader replaces
        matching pixels with NaNs. 

        When we know we do not care and need numpy arrays (e.g. when building stacks)
        call this to replace NaNs with 0. 
        """
        
        arr = image.getArray()
        nanmask = np.isnan(arr)
        if nanmask.any():
            # I feel so dirty -- CPL
            arr[nanmask] = 0.0
        return arr

    def getRawDataImageArray(self, readNum):
        """Get an data image as an numpy array, with NaNs replaced by 0."""

        img = self.getRawDataImage(readNum)
        return self.evilGetArray(img)
    
    def getRawIrpImageArray(self, readNum):
        """Get an IRP image as an numpy array, with NaNs replaced by 0."""

        img = self.getRawIrpImage(readNum)
        return self.evilGetArray(img)
    
    def refpix4BorderCorrect(self, image, colWindow=4,
                             doRows=True, doCols=True):
        """This is the "standard" Teledyne 'refPixel4' border reference pixel scheme. 

        Step 1:
           For each channel, average all 8 top&bottom rows to one number.
           Subtract that from the channel.

        Step 2:
            Take a 9-row running average of the left&right columns.
            Subtract that from each row.

        The output is not lovely, in particular in the face of per-channel level drift and 1/f. But
        we are just here to duplicate the standard logic. Do *not* modify/improve this particular 
        routine. In general don't use it, either.
        """

        im = image.array
        imHeight, imWidth = im.shape
        if imHeight != self.h4Size or imWidth != self.h4Size:
            raise RuntimeError('not ready to deal with non-full images')

        if doCols:
            top = im[imHeight-self.refPix:,:]
            bottom = im[:self.refPix, :]

            chanWidth = imWidth // self.nchan
            for c_i in range(self.nchan):
                xlow = c_i * chanWidth
                xhigh = xlow + chanWidth
                # Leave ref columns unmolested.
                if c_i == 0:
                    xlow = self.refPix
                elif c_i == self.nchan-1:
                    xhigh = imWidth-self.refPix
                ichan = im[:, xlow:xhigh]

                chanOffset = (top[:, xlow:xhigh].mean()
                              + bottom[:, xlow:xhigh].mean()) / 2
                ichan -= chanOffset

        if doRows:
            sideRefImage = np.ndarray((imHeight, self.refPix*2), dtype=im.dtype)
            sideRefImage[:, :self.refPix] = im[:, :self.refPix]
            sideRefImage[:, self.refPix:] = im[:, imWidth-self.refPix:]
            sideCorr = np.zeros((imHeight,1))

            # Not quite right at the ends -- should not be including the reference rows. 
            for row_i in range(colWindow, imHeight-colWindow+1):
                sideCorr[row_i] = sideRefImage[row_i-colWindow:row_i+colWindow,:].mean()
                im[row_i, :] -= sideCorr[row_i]

        return image

    def _chanFlattenFit(self, img, deg=1):
        """remove per-channel tilt+offset

        Turns out that some ASIC channels have significant drifts. This is just
        a convenience routine to flatten that out. So far linear has been good enough.
        """
        h, w = img.shape
        res = img.copy()
        chan_w = w // self.nchan
        assert(w == chan_w * self.nchan)

        yy = np.arange(h)
        for i_c in range(self.nchan):
            pix_i = np.arange(i_c*chan_w, (i_c+1)*chan_w)
            chanVec = np.median(img[:, pix_i], axis=1)
            fit = np.polynomial.Polynomial.fit(yy, chanVec, deg=deg)
            chan_off = fit(yy)
            res[:, pix_i] -= chan_off[:,None]
        return res

    def calcIrpDefects(self, r0=0, r1=-1, noiseClip=None, 
                      debug=False, force=False):
        """Given two IRP images, return a list of bad columns.

        There are bad pixels in the single IRP row. Stuck high or low, or just noisy. This
        is a crude filter to identify them on the fly, and is useful to seed curated masks.

        Individual IRP frames have significant common structure, so we diff two, just as we do for 
        the science pixels and thus the final CDS images.

        This needs to be done per channel, due to bad ASIC channels on at least n1/18660 or n3/18321. Grr.

        Args
        ----
        r0,r1 : `int`
            0-based indexes of the two IRP images to diff.
        noiseClip : `float`
            The threshold for identifying bad columns. 
            n.b. there are some quite noisy channels from some of the ASICs.
        """

        if noiseClip is None:
            noiseClip = self.corrections[self.h4Serial]['irpNoiseClip']

        i0 = self.getRawIrpImageArray(r0)
        i1 = self.getRawIrpImageArray(r1)
        dIrp = i1-i0

        # Some 18660/18321 ASIC channels have varying slopes: flatten them out a bit.
        dIrp = self._chanFlattenFit(dIrp)
        irpNoise = np.std(dIrp, axis=0) 
        mask = (irpNoise >= noiseClip) | (irpNoise < 4)

        # Local hacks
        # 18321/119 has a bad ASIC channel which looks like a clot of bad columns. For
        # now just ignore real bad columns.
        if self.h4Serial == 18321:
            chanWidth = 128
            mask[7*chanWidth:8*chanWidth] = False

        mask = np.where(mask)[0]

        if debug:
            return mask, irpNoise
        else:
            return mask

    def interpolateChannelIrp(self, rawChan, doFlip):
        """Given an IRPn channel from a PFSB file, return IRP1 channel image.

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

    def calcDataFeedthrough(self, r0=0, r1=-1, debug=False):
        """Fit per-channel feedthrough. Use a flat or something with features: we are trying to flatten the result"""

        r0 = self.positiveIndex(r0)
        r1 = self.positiveIndex(r1)    
        if r1-r0 < 1:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}) + 1:')

        image0 = self.getRawDataImageArray(r0)
        image1 = self.getRawDataImageArray(r1)

        ref0 = self.getRawIrpImageArray(r0)
        ref1 = self.getRawIrpImageArray(r1)

        def fitFunc(x, c, r0, r1, i0, i1):
            xslice = slice(c*128, (c+1)*128)
            d0 = r0[xslice, :] - i0[xslice, :] * x
            d1 = r1[xslice, :] - i1[xslice, :] * x
            dd = d1 - d0

            pLo, pHi = np.nanpercentile(dd, (10, 90))
            if debug:
                print(f'{c}   {x:0.5f} {pHi-pLo:0.5f} {pLo:0.5f} {pHi:0.5f}')
            return pHi-pLo
        
        res = dict()
        for c in range(32):

            ret = scipy.optimize.minimize_scalar(fitFunc, bracket=[0.0, 0.03],
                                                 args=(c, ref0, ref1, image0, image1))
            res[c] = ret.x
        return res
    
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
          for it to be the last pixel in the group.

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

    def _getFinalDiffIrp(self, rawDiffIrp, filterWidth=31, 
                         badCols=None, useFft=True):
        """Given a simple IRP difference image, return the final IRP image.

        We do a few things here:
         - replace the columns from known bad IRP row pixels with the median of their channel
         - per-channel, smooth the IRP image rows over {filterWidth} pixels

        For the moment, we assume that the following has been done before
        we are called:
          - IRPn images have been filled out into IRP1 images.
          - Any crosstalk/leakage from the data images has been removed.
          
        Args
        ----
        rawDiffIrp : `numpy.ndarray`
            The raw difference IRP image. Full-size (i.e. IRP1)
        filterWidth : `int`
            The width of the smoothing filter to use. Should be odd.
        badcols : `numpy.ndarray`, optional
            The IRP columns to mask/interpolate over.
            Here just for development. Will always come from curated data.

        Returns
        -------
        finalIrp : `numpy.ndarray`
            The processed IRP image.            

           """
        
        h, w = rawDiffIrp.shape
        chan_w = w//self.nchan
        ww_half = filterWidth//2
        padWidth = ww_half*2
        
        if badCols is None:
            badCols = self.irpDefects

        # Construct a cos bell filter to run over each channel. 
        # Pad filter with zeros to same width as padded channel.
        if filterWidth > 0:
            filtCore = scipy.signal.hann(filterWidth)
            filtCore /= filtCore.sum()
            filt1 = np.zeros(shape=(chan_w+padWidth), dtype='f4')
            fc = (chan_w+filterWidth)//2
            filt1[fc-ww_half:fc+ww_half+1] = filtCore
            filt = np.tile(filt1, (h, 1))

        # scratch space for the padded channel
        chan1 = np.zeros(shape=(h, chan_w+padWidth), dtype=rawDiffIrp.dtype)
        out = np.zeros_like(rawDiffIrp)
        t0 = time.time()
        for i_c in range(self.nchan):
            colLow = i_c*chan_w
            colHigh = (i_c + 1)*chan_w
            chanPixIdx = np.arange(colLow, colHigh)
            chan0 = rawDiffIrp[:, chanPixIdx].copy()
            if np.any(np.isnan(chan0)):
                print(f'nan in chan0 {i_c} {np.isnan(chan0).sum()}')
            if badCols is not None:
                badcolsChan = badCols[(badCols >= colLow) & (badCols < colHigh)] - colLow
                if len(badcolsChan) > 0:
                    # print(f'{i_c} {badcolsChan} {badcolsChan + colLow}')
                    chan0[:, badcolsChan] = np.nanmedian(chan0, axis=1, keepdims=True)
                    
            if filterWidth > 0:
                # Pad channel with median of channel out to half the filter width
                chan1[:, ww_half:chan_w+ww_half] = chan0
                chan1[:, :ww_half] = np.nanmedian(chan0[:, :ww_half], axis=1, keepdims=True) 
                chan1[:, -ww_half:] = np.nanmedian(chan0[:, -ww_half:], axis=1, keepdims=True)
                if np.any(np.isnan(chan1)):
                    print(f'nan in chan1 {i_c} {np.isnan(chan1).sum()}')

                if useFft:
                    chan_f = scipy.signal.fftconvolve(chan1, filt, mode='same', axes=1)
                else:
                    chan_f = scipy.signal.oaconvolve(chan1, filt, mode='same', axes=1)
                out[:, chanPixIdx] = chan_f[:, ww_half:chan_w+ww_half]
            else:
                out[:, chanPixIdx] = chan0
        t1 = time.time()
        self.logger.debug(f' final IRP: {t1-t0:0.3f}s {filt.shape if filterWidth > 0 else 0} {chan1.shape} raw={rawDiffIrp.sum():0.1f} fac={out.sum()/rawDiffIrp.sum():0.5f}')
        
        return out

    def getSingleRead(self, readNum=0, useIrp=True):
        """Get a single semi-reference-corrected read from a ramp.

        Args
        ----
        readNum : `int`
            The read to return. 0-indexed.
        useIrp: `bool`
            If the ramp has IRP, do use it. Else use the border pixels.

        Returns
        -------
        image : `afw.image.maskedImage`
            reference-corrected image.
        """
        readNum = self.positiveIndex(readNum)
        image = self.getRawDataImage(readNum)
        if self.haveIrp and useIrp:
            ref = self.getRawIrpImage(readNum)
            ref = self.constructFullIrp(ref)
            image -= ref
        else:
            image = self.borderCorrect(image) 

        return image

    def applyLeakageCorrection(self, ref, image, leakageFactor=None):
        if leakageFactor is None:
            leakageFactor = self.dataToIrpLeakageFactor

        for c in range(32):
            xslice = slice(c*128, (c+1)*128)
            ref[xslice, :] -= image[xslice, :] * leakageFactor[c]

        return ref
        
    def _getDiffIrp(self, r0=0, r1=-1, leakageFactor=None,
                    filterWidth=31, badCols=None):
        """Just for development. Get a couple of steps of the difference IRP image."""

        if leakageFactor is None:
            leakageFactor = self.dataToIrpLeakageFactor

        r0 = self.positiveIndex(r0)
        r1 = self.positiveIndex(r1)    
        if r1-r0 < 1:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}) + 1:')

        image0 = self.getRawDataImageArray(r0)
        image1 = self.getRawDataImageArray(r1)

        ref0 = self.getRawIrpImageArray(r0)
        ref1 = self.getRawIrpImageArray(r1)
        dref0 = ref1 - ref0   
        # Naive img - ref

        if leakageFactor is not None:
            ref0 = self.applyLeakageCorrection(ref0, image0, leakageFactor)
            ref1 = self.applyLeakageCorrection(ref1, image1, leakageFactor)
        dref1 = ref1 - ref0
        dref2 = self._getFinalDiffIrp(dref1, filterWidth=0, badCols=badCols)
        # Naive correction plus leakage corrections and masking bad columns

        drefFinal = self._getFinalDiffIrp(dref1, filterWidth=filterWidth, badCols=badCols)
        # Also spacially smooth IRP image

        return dref0, dref2, drefFinal
        
    def getCds(self, r0=0, r1=-1, useIrp=True, filterIRP=True,
               useLeakageFactor=None):
        """Get a single reference-corrected CDS from a ramp.

        Should leave as MaskedImage, both to let errors propagate and to allow indicating bad IRP columns

        Args
        ----
        r0 : `int`
            The read to start from. 0-indexed.
        r1 : `int`
            The read to end at. 0-indexed.
        useIRP : `bool`
            If the ramp has IRP, do use it. Else use the border pixels.
        filterIRP : `bool`
            apply known corrections to the IRP images
        useLeakageFactor : `float`, optional
            If not None, use this value for the data->IRP leakage factor. Only for development
                   
        Returns
        -------
        image : `afw.image.ImageF`
            reference-corrected image.
        """

        r0 = self.positiveIndex(r0)
        r1 = self.positiveIndex(r1)    
        if r1-r0 < 1:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}) + 1:')

        leakageFactor = useLeakageFactor if useLeakageFactor is not None else self.dataToIrpLeakageFactor

        t0 = time.time()
        image0 = self.getRawDataImageArray(r0)
        image1 = self.getRawDataImageArray(r1)

        if self.haveIrp and useIrp:
            ref0 = self.getRawIrpImageArray(r0)
            ref1 = self.getRawIrpImageArray(r1)

            if filterIRP and leakageFactor is not None:
                ref0 = self.applyLeakageCorrection(ref0, image0, leakageFactor)
                ref1 = self.applyLeakageCorrection(ref1, image1, leakageFactor)
            dref = ref1 - ref0
            if filterIRP:
                dref = self._getFinalDiffIrp(dref)
            else:
                dref = self._getFinalDiffIrp(dref, filterWidth=0)

            image1 -= image0
            image1 -= dref
            return afwImage.ImageF(image1)
        else:
            image0 = self.borderCorrect(image0)
            image1 = self.borderCorrect(image1)
            image1 -= image0

        return afwImage.ImageF(image1)

    def _simpleStack(self, reader=None, r0=0, r1=-1, step=1, subtractR0=False):
        """Return all the raw frames in a single 3d stack.

        Args
        ----
        reader : `callable`
          A function to read the images. If None, uses `self.getImage`.
        r0 : `int`
          The 0-indexed read to start from.
        r1 : `int`
          The 0-indexed read to end with. Inclusive.
        step : `int`
          A `range`-compliant stepping between r0 and r1.

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
        nreads = r1 - r0 + 1

        if r1 <= r0:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}):')

        stack = None
        for r_i in range(r0, r1+1, step):
            readImg = reader(r_i)
            if stack is None:
                h, w = readImg.image.array.shape
                stack = np.empty(shape=(nreads, h, w), dtype='f4')
            stack[r_i,:,:] = self.evilGetArray(readImg)

        if subtractR0:
            stack -= stack[0]

        return stack

    def rawDataStack(self, r0=0, r1=-1, step=1, subtractR0=False):
        """Return all the raw image frames in a single 3d stack."""
        return self._simpleStack(reader=self.getRawDataImage, 
                                 r0=r0, r1=r1, step=step, subtractR0=subtractR0)
    def rawIrpStack(self, r0=0, r1=-1, step=1, subtractR0=False):
        """Return all the raw IRP frames in a single 3d stack."""
        return self._simpleStack(reader=self.getRawIrpImage, 
                                 r0=r0, r1=r1, step=step, subtractR0=subtractR0)

    def cdsStack(self, r0=0, r1=-1, step=1, filterIrp=True, bbox=None,
                 useLeakageFactor=None, debug=False):
        """Return all the fully IRP-corrected frames in a single 3d stack.

        Given two raw data images d0 and d1, and two raw IRP images i0 and i1, the net CDS image
        can be either (d1 - i1) - (d0 - i0), or (d1 - d0) - (i1 - i0). The IRP row has various 
        artifacts which make using the latter "nicer", or at least easier to make sense of. So that 
        is the way we do it.
        In particular there are:
          - bad IRP row pixels
          - pixel-to-pixel offsets
          - bad ASIC channels
          
        Note that there is also up to ~1% printthrough from data->IRP, and presumably from IRP->data
        
        Args
        ----
        r0 : `int`
          The 0-indexed read to start from.
        r1 : `int`
          The 0-indexed read to end with. Inclusive.
        step : `int`
          A `range`-compliant stepping between r0 and r1.
        filterIrp : `bool`
            Whether to apply known corrections to the IRP images.
        bbox : `lsst.geom.Box2I`
            UNUSED. The region of the image to return. If None, return the whole image.   
        useLeakageFactor : `float`, optional
            If not None, use this value for the data->IRP leakage factor. Only for development

        Returns
        -------
        stack : 3-d float32 numpy array
           the CDS stack, with axis 0 being the reads. Since the images are CDS,
           one fewer than the number of selected reads. 
        """

        if bbox is not None:
            raise NotImplementedError('bbox not yet implemented')

        leakageFactor = useLeakageFactor if useLeakageFactor is not None else self.dataToIrpLeakageFactor

        r0 = self.positiveIndex(r0)
        r1 = self.positiveIndex(r1)
        if r1-r0 < 1:
            raise ValueError(f'r1 ({r1}) must be greater than r0 ({r0}) + 1:')
        nreads = r1 - r0 + 1

        # Grab the components of read 0, which we will subtract from all the others.
        irp0 = self.getRawIrpImageArray(r0)
        data0 = self.getRawDataImageArray(r0)
        if leakageFactor is not None:
            irp0 = self.applyLeakageCorrection(irp0, data0, leakageFactor)

        stack = np.empty(shape=(nreads-1, *data0.shape), dtype='f4')
        for r_i in range(r0+1, r1+1, step):
            t0 = time.time()
            data1 = self.getRawDataImageArray(r_i)
            irp1 = self.getRawIrpImageArray(r_i)
            t1 = time.time()
            if leakageFactor is not None:
                irp1 = self.applyLeakageCorrection(irp1, data1)
            dirp = irp1 - irp0
            ddata = data1 - data0
            if filterIrp:
                ddata -= self._getFinalDiffIrp(dirp)
            else:
                ddata -= self._getFinalDiffIrp(dirp, filterWidth=0)
            stack[r_i-1,:,:] = ddata
            t2 = time.time()
            if debug:
                print(f'cds {r_i} io1={t1-t0:0.3f} proc={t2-t1:0.3f}')
        return stack

def getCds(butler, dataId, r0=0, r1=-1, useIrp=True):
    ramp = H4Ramp(butler, dataId)
    return ramp.getCds(r0=r0, r1=r1, useIrp=useIrp)
