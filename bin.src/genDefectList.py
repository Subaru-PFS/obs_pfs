#!/usr/bin/env python
import argparse
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.persistence as dafPersist
from lsst.log import Log
import numpy as np
import os
import pfs.drp.stella as drpStella
try:
    import lsst.afw.display as afwDisplay
except ImportError:
    afwDisplay = None

parser = argparse.ArgumentParser()

parser.add_argument("--badCols", help="List of bad columns delimited by ','", type=str, default="1020,1021,1022,1023,1024, 1025, 1026,3068,3069,3070,3071,3072,3073, 3074")
parser.add_argument("--ccd", help="CCD number for which to create the defect list", type=int, default=4)
parser.add_argument("--display", help="Set to display outputs", action="store_true")
parser.add_argument("--debug", help="Set to print debugging outputs", action="store_true")
parser.add_argument("--expTimeHigh", help="Exposure time for stronger flats", type=int, default=30)
parser.add_argument("--expTimeLow", help="Exposure time for weaker flats", type=int, default=2)
parser.add_argument("--gapHigh", help="Number of column where gap between physical devices ends", type=int, default=2052)
parser.add_argument("--gapLow", help="Number of column where gap between physical devices starts", type=int, default=2045)
parser.add_argument("--maxSigma", help="Maximum sigma for good pixels", type=float, default=6.)
parser.add_argument("--medianFlatsOut", help="Median Flats output root (without '.fits'). Leave empty for not writing", type=str, default="")
parser.add_argument("--nCols", help="Number of columns in images", type=int, default=4096)
parser.add_argument("--nPixCut", help="Number of rows/columns around the edge of images to ignore", type=int, default=40)
parser.add_argument("--nRows", help="Number of rows in images", type=int, default=4174)
parser.add_argument("--outFile", help="Output defect list relative to OBS_PFS_DIR", type=str, default="pfs/defects/2015-12-01/defects.dat")
parser.add_argument("--rerun", help="Which rerun directory are the postISRCCD images in if not in root?", type=str, default="")
parser.add_argument("--root", help="Path to PFS data directory", type=str, default="/Volumes/My Passport/Users/azuri/spectra/pfs/PFS")
parser.add_argument("--verbose", help="Set to print info outputs", action="store_true")
parser.add_argument("--visitLow", help="Lowest visit number to search for flats", type=int, default=6301)
parser.add_argument("--visitHigh", help="Highest visit number to search for flats", type=int, default=6758)
args = parser.parse_args()

if afwDisplay:
    if os.environ.get('DISPLAY_DS9_DIR') != "":
        try:
            afwDisplay.setDefaultBackend("ds9" if args.display else "virtualDevice")
        except RuntimeError as e:
            print e

    afwDisplay.setDefaultMaskTransparency(75)

# Quiet down lots of loggers, so we can see genDefectList logs better
Log.getLogger("daf.persistence.butler").setLevel(Log.WARN)
Log.getLogger("afw.image.Mask").setLevel(Log.WARN)
Log.getLogger("afw.ExposureFormatter").setLevel(Log.WARN)
Log.getLogger("daf.persistence.LogicalLocation").setLevel(Log.WARN)
Log.getLogger("CameraMapper").setLevel(Log.WARN)
Log.getLogger("afw.image.ExposureInfo").setLevel(Log.WARN)

logger = Log.getLogger("genDefectList")
logger.setLevel(Log.INFO if args.verbose else Log.WARN)
if args.debug:
    logger.setLevel(Log.DEBUG)

outFileName = os.path.join(os.environ.get('OBS_PFS_DIR'),args.outFile)

# Read existing defects file
textArrRead = False
if args.outFile != '':
    try:
        file = open(outFileName, 'r')
        textArr = file.readlines()
        textArrRead = True
        logger.debug("Existing outFile found")
    except IOError as e:
        logger.debug("No existing outFile found")
    text_file = open(outFileName, "w")
    text_file.write("#\n")
    text_file.write("# Defects file\n")
    text_file.write("#\n")
    text_file.write("# Convert to a fits table using getDefectFits.py;  this is done for you by running scons\n")
    text_file.write("#\n")
    text_file.write("#CCD x0    y0    width height\n")

if textArrRead:
    firstChars = []
    for i in range(len(textArr)):
        firstChars.append(textArr[i][0])
    for i in np.arange(len(firstChars)-1,-1,-1):
        if firstChars[i] == '#':
            textArr.pop(i)
            firstChars.pop(i)
    uniqueFirstChars = np.unique(firstChars)
    uniqueFirstCharsSorted = np.sort(uniqueFirstChars)
    for a in uniqueFirstCharsSorted:
        if int(a) < args.ccd:
            for i in range(len(textArr)):
                if firstChars[i] == a:
                    text_file.write(textArr[i])

badCols = [int(item) for item in args.badCols.split(',')]
dataDir = args.root
if args.rerun != "":
    dataDir = os.path.join(args.root, 'rerun')
    dataDir = os.path.join(dataDir, args.rerun)
butler = dafPersist.Butler( dataDir )

biasesVisit = list()
flatsVisit = list()
flatsExpTimes = list()
flatsLow = list()
flatsHigh = list()
for visit in np.arange( args.visitLow, args.visitHigh + 1 ):
    try:
        dataId = dict( ccd=args.ccd, visit=visit )
        fileName = butler.get( 'raw_filename', dataId )
        md = afwImage.readMetadata( fileName[ 0 ] )

        if md.get("IMAGETYP") == "bias":
            biasesVisit.append(visit)
        if md.get("IMAGETYP") == "flat":
            flatsVisit.append(visit)
            flatsExpTimes.append(md.get('EXPTIME'))
            if md.get('EXPTIME') == args.expTimeLow:
                flatsLow.append(visit)
            elif md.get('EXPTIME') == args.expTimeHigh:
                flatsHigh.append(visit)
    except:
        logger.debug("visit %d not found" % visit)

dataId = dict(visit=flatsLow[0], ccd=args.ccd)
image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
allFlatsLow = np.ndarray(shape=[image.shape[0], image.shape[1], len(flatsLow)], dtype='float32')
allFlatsHigh = np.ndarray(shape=[image.shape[0], image.shape[1], len(flatsHigh)], dtype='float32')
for i in np.arange(len(flatsLow)):
    dataId = dict(visit=flatsLow[i], ccd=args.ccd)
    image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
    allFlatsLow[:,:,i] = image
for i in np.arange(len(flatsHigh)):
    dataId = dict(visit=flatsHigh[i], ccd=args.ccd)
    image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
    allFlatsHigh[:,:,i] = image

nLow = len(flatsLow)
if nLow > 2:
    medianFlatsLow = np.median(allFlatsLow,2)
    logger.info('%d Flats with %ds exposure time found' % (nLow, args.expTimeLow))
elif nLow == 2:
    medianFlatsLow = np.min(allFlatsLow,2)
    logger.warn('Only %d Flats with %ds exposure time found' % (nLow, args.expTimeLow))
elif nLow < 1:
    raise RunTimeError("No Flats with %d s exposure time found" % args.expTimeLow)
else:
    logger.warn('Only 1 Flat with %ds exposure time' % args.expTimeLow)
    medianFlatsLow = allFlatsLow[:,:,0]

nHigh = len(flatsHigh)
if nHigh > 2:
    medianFlatsHigh = np.median(allFlatsHigh,2)
    logger.info('%d Flats with %ds exposure time found' % (nLow, args.expTimeHigh))
elif nHigh == 2:
    medianFlatsHigh = np.min(allFlatsHigh,2)
    logger.warn('Only %d Flats with %ds exposure time found' % (nLow, args.expTimeHigh))
elif nHigh < 1:
    raise RunTimeError("No flats with %d s exposure time found" % args.expTimeHigh)
else:
    logger.warn('Only 1 Flat with %ds exposure time' % args.expTimeHigh)
    medianFlatsHigh = allFlatsHigh[:,:,0]

medianFlatsLow = drpStella.where( medianFlatsLow, '==', 0., 0.000001, medianFlatsLow)
medianFlatsHigh = drpStella.where( medianFlatsHigh, '==', 0., 0.000001, medianFlatsHigh)
divFlat = medianFlatsHigh / medianFlatsLow

sctrl = afwMath.StatisticsControl()
sctrl.setWeighted(False)

divFlatIm = afwImage.ImageF(divFlat)
divFlatMIm = afwImage.makeMaskedImage(divFlatIm)

mean = afwMath.makeStatistics(divFlatMIm, afwMath.MEANCLIP).getValue()
stddev = afwMath.makeStatistics(divFlatMIm,afwMath.STDEVCLIP).getValue()

if args.outFile != '':
    mask = divFlatMIm.getMask()
    maskVal = 1 << divFlatMIm.getMask().getMaskPlane('BAD')
    for badCol in badCols:
        text_file.write("%d 0 %d 1 %d\n" % (args.ccd, badCol, args.nRows))
        mask[badCol,:] |= maskVal
    nBad = 0
    for iRow in np.arange(args.nPixCut,divFlat.shape[0]-args.nPixCut):
        for iCol in np.arange(args.nPixCut,args.gapLow):
            if iCol not in badCols:
                if np.absolute(divFlat[iRow, iCol] - mean) > args.maxSigma * stddev:
                    nBad = nBad + 1
                    text_file.write("%d %d %d 1 1\n" % (args.ccd, iCol, iRow))
                    mask[iCol, iRow] |= maskVal;

    for iRow in np.arange(args.nPixCut,divFlat.shape[0]-args.nPixCut):
        for iCol in np.arange(args.gapHigh,args.nCols-args.nPixCut):
            if iCol not in badCols:
                if np.absolute(divFlat[iRow, iCol] - mean) > args.maxSigma * stddev:
                    nBad = nBad + 1
                    text_file.write("%d %d %d 1 1\n" % (args.ccd, iCol, iRow))
                    mask[iCol, iRow] |= maskVal;

    if textArrRead:
        for a in uniqueFirstCharsSorted:
            if int(a) > args.ccd:
                for i in range(len(textArr)):
                    if firstChars[i] == a:
                        text_file.write(textArr[i])
    text_file.write("# Dead amps")
    text_file.close()

    logger.info('%d bad pixels found' % nBad)

if args.medianFlatsOut != '' or args.display:
    medianFlatsLowExp = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(medianFlatsLow)))
    medianFlatsLowExp.getMaskedImage().getMask()[:,:] = mask[:,:]
    medianFlatsHighExp = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(medianFlatsHigh)))
    medianFlatsHighExp.getMaskedImage().getMask()[:,:] = mask[:,:]

if args.medianFlatsOut != '':
    medianFlatsLowExp.writeFits(args.medianFlatsOut+"Low.fits")
    medianFlatsHighExp.writeFits(args.medianFlatsOut+"High.fits")
    divFlatMIm.writeFits(args.medianFlatsOut+"DivHighByLow.fits")

if afwDisplay:
    if args.display:
        divFlatExp = afwImage.makeExposure(divFlatMIm)

        display0 = afwDisplay.getDisplay()
        display0.mtv(medianFlatsLowExp, title="parent")
        display0.setMaskTransparency(50)
        display0.setMaskPlaneColor("BAD", "orange")
