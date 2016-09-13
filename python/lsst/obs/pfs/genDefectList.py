import lsst.daf.persistence as dafPersist
import lsst.afw.image as afwImage
import pfs.drp.stella as drpStella
import lsst.afw.math as afwMath
import numpy as np
import os
try:
    import lsst.afw.display as afwDisplay
except ImportError:
    afwDisplay = None

if afwDisplay:
    try:
        afwDisplay.setDefaultBackend("ds9" if True else "virtualDevice")
    except RuntimeError as e:
        print e

    afwDisplay.setDefaultMaskTransparency(75)

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-homeDir", help="Path to PFS data directory", default="/Volumes/My Passport/Users/azuri/spectra/pfs/PFS")
parser.add_argument("-outFile", help="Output defect list relative to OBS_PFS_DIR", default="pfs/defects/2015-12-01/defects.dat")
parser.add_argument("-medianFlatsOut", help="Median Flats output root (without '.fits'). Leave empty for not writing", default="/Users/azuri/spectra/pfs/medianFlat")# - leave empty if you don't want to write the median Flats
parser.add_argument("-ccd", help="CCD number for which to create the defect list", type=int, default=5)
#parser.add_argument("display", help="Set to display outputs", action="store_true")
parser.add_argument("-visitLow", help="Lowest visit number to search for flats", type=int, default=6301)
parser.add_argument("-visitHigh", help="Highest visit number to search for flats", type=int, default=6758)
parser.add_argument("-expTimeLow", help="Exposure time for weaker flats", type=int, default=2)
parser.add_argument("-expTimeHigh", help="Exposure time for stronger flats", type=int, default=30)
parser.add_argument("-nPixCut", help="Number of rows/columns around the edge of images to ignore", type=int, default=40)
parser.add_argument("-nCols", help="Number of columns in images", type=int, default=4096)
parser.add_argument("-nRows", help="Number of rows in images", type=int, default=4174)
parser.add_argument("-gapLow", help="Number of column where gap between physical devices starts", type=int, default=2045)
parser.add_argument("-gapHigh", help="Number of column where gap between physical devices ends", type=int, default=2052)
parser.add_argument("-badCols", help="List of bad columns delimited by ','", type=str, default="1020,1021,1022,1023,1024, 1025, 1026,3068,3069,3070,3071,3072,3073, 3074")
args = parser.parse_args()

badCols = [int(item) for item in args.badCols.split(',')]
butler = dafPersist.Butler( args.homeDir )

biasesVisit = list()
flatsVisit = list()
flatsExpTimes = list()
flatsLow = list()
flatsHigh = list()
for visit in np.arange( args.visitLow, args.visitHigh + 1 ):
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

print 'flatsVisit = ',flatsVisit
print 'flatsExpTimes = ',flatsExpTimes
print 'flatsLow = ',len(flatsLow),': ',flatsLow
print 'flatsHigh = ',len(flatsHigh),': ',flatsHigh

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
    
if len(flatsLow) > 2:
    medianFlatsLow = np.median(allFlatsLow,2)
elif len(flatsLow) == 2:
    medianFlatsLow = np.min(allFlatsLow,2)
else:
    print 'only one Flat with ',args.expTimeLow,'s exposure time'
    medianFlatsLow = allFlatsLow[:,:,0]
if len(flatsHigh) > 2:
    medianFlatsHigh = np.median(allFlatsHigh,2)
elif len(flatsHigh) == 2:
    medianFlatsHigh = np.min(allFlatsHigh,2)
else:
    print 'only one Flat with ',args.expTimeHigh,'s exposure time'
    medianFlatsHigh = allFlatsHigh[:,:,0]

medianFlatsLow = drpStella.where( medianFlatsLow, '==', 0., 0.000001, medianFlatsLow)
medianFlatsHigh = drpStella.where( medianFlatsHigh, '==', 0., 0.000001, medianFlatsHigh)
            
if args.medianFlatsOut != '':
    afwImage.makeImageFromArray(medianFlatsLow).writeFits(args.medianFlatsOut+"Low.fits")
    afwImage.makeImageFromArray(medianFlatsHigh).writeFits(args.medianFlatsOut+"High.fits")

divFlat = medianFlatsHigh / medianFlatsLow

if args.medianFlatsOut != '':
    afwImage.makeImageFromArray(divFlat).writeFits(args.medianFlatsOut+"DivHighByLow.fits")
sctrl = afwMath.StatisticsControl()
sctrl.setWeighted(False)

divFlatIm = afwImage.ImageF(divFlat)
divFlatMIm = afwImage.makeMaskedImage(divFlatIm)

mean = afwMath.makeStatistics(divFlatMIm, afwMath.MEANCLIP).getValue()
print 'mean = ',mean

stddev = afwMath.makeStatistics(divFlatMIm,afwMath.STDEVCLIP).getValue()
print 'stddev = ',stddev

if args.outFile != '':
    text_file = open(os.path.join(os.environ.get('OBS_PFS_DIR'),args.outFile), "w")
    text_file.write("#\n")
    text_file.write("# Defects file\n")
    text_file.write("#\n")
    text_file.write("# Convert to a fits table using getDefectFits.py;  this is done for you by running scons\n")
    text_file.write("#\n")
    text_file.write("#CCD x0    y0    width height\n")
    for badCol in badCols:
        text_file.write("%d 0 %d 1 %d\n" % (args.ccd, badCol, args.nRows))
        divFlatMIm.getMask().getArray()[:,badCol] = 1
    nBad = 0
    for iRow in np.arange(args.nPixCut,divFlat.shape[0]-args.nPixCut):
        for iCol in np.arange(args.nPixCut,args.gapLow):
            if iCol not in badCols:
                if np.absolute(divFlat[iRow, iCol] - mean) > 5. * stddev:
                    nBad = nBad + 1
                    text_file.write("%d %d %d 1 1\n" % (args.ccd, iRow, iCol))
                    divFlatMIm.getMask().getArray()[iRow,iCol] = 1

    for iRow in np.arange(args.nPixCut,divFlat.shape[0]-args.nPixCut):
        for iCol in np.arange(args.gapHigh,args.nCols-args.nPixCut):
            if iCol not in badCols:
                if np.absolute(divFlat[iRow, iCol] - mean) > 5. * stddev:
                    nBad = nBad + 1
                    text_file.write("%d %d %d 1 1\n" % (args.ccd, iRow, iCol))
                    divFlatMIm.getMask().getArray()[iRow, iCol] = 1

    text_file.write("# Dead amps")
    text_file.close()

    print nBad,' bad pixels found'

if afwDisplay:
    if True:#args.display:
        divFlatExp = afwImage.makeExposure(divFlatMIm)
        medianFlatsLowExp = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(medianFlatsLow)))
        medianFlatsLowExp.getMaskedImage().getMask().getArray()[:,:] = divFlatExp.getMaskedImage().getMask().getArray()[:,:]

        display0 = afwDisplay.getDisplay()
        print 'display0 = ',display0
        display0.mtv(medianFlatsLowExp, title="parent")
        display0.setMaskTransparency(10)
        display0.setMaskPlaneColor("CROSSTALK", "orange")
