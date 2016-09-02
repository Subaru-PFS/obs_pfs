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

homeDir = '/Volumes/My Passport/Users/azuri/spectra/pfs'
outFile = os.path.join(os.environ.get('OBS_PFS_DIR'),"pfs/defects/2015-12-01/defects.dat")
medianFlatsOut = ''#'/Users/azuri/spectra/pfs/medianFlat'# - leave empty if you don't want to write the median Flats

ccd = 5
display = False

visitLow = 6301
visitHigh = 6758
expTimeLow = 2
expTimeHigh = 30
nPixCut = 40
nCols = 4096
gapLow = 2045
gapHigh = 2052

butler = dafPersist.Butler(os.path.join(homeDir,"PFS"))

biasesVisit = list()
flatsVisit = list()
flatsExpTimes = list()
flatsLow = list()
flatsHigh = list()
for visit in np.arange( visitLow, visitHigh + 1 ):
    dataId = dict( ccd=ccd, visit=visit )
    fileName = butler.get( 'raw_filename', dataId )
    md = afwImage.readMetadata( fileName[ 0 ] )

    if md.get("IMAGETYP") == "bias":
        biasesVisit.append(visit)
    if md.get("IMAGETYP") == "flat":
        flatsVisit.append(visit)
        flatsExpTimes.append(md.get('EXPTIME'))
        if md.get('EXPTIME') == expTimeLow:
            flatsLow.append(visit)
        elif md.get('EXPTIME') == expTimeHigh:
            flatsHigh.append(visit)

print 'flatsVisit = ',flatsVisit
print 'flatsExpTimes = ',flatsExpTimes
print 'flatsLow = ',len(flatsLow),': ',flatsLow
print 'flatsHigh = ',len(flatsHigh),': ',flatsHigh

dataId = dict(visit=flatsLow[0], ccd=ccd)
image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
allFlatsLow = np.ndarray(shape=[image.shape[0], image.shape[1], len(flatsLow)], dtype='float32')
allFlatsHigh = np.ndarray(shape=[image.shape[0], image.shape[1], len(flatsHigh)], dtype='float32')
for i in np.arange(len(flatsLow)):
    dataId = dict(visit=flatsLow[i], ccd=5)
    image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
    allFlatsLow[:,:,i] = image
for i in np.arange(len(flatsHigh)):
    dataId = dict(visit=flatsHigh[i], ccd=5)
    image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
    allFlatsHigh[:,:,i] = image
    
if len(flatsLow) > 2:
    medianFlatsLow = np.median(allFlatsLow,2)
elif len(flatsLow) == 2:
    medianFlatsLow = np.min(allFlatsLow,2)
else:
    print 'only one Flat with ',expTimeLow,'s exposure time'
    medianFlatsLow = allFlatsLow[:,:,0]
if len(flatsHigh) > 2:
    medianFlatsHigh = np.median(allFlatsHigh,2)
elif len(flatsHigh) == 2:
    medianFlatsHigh = np.min(allFlatsHigh,2)
else:
    print 'only one Flat with ',expTimeHigh,'s exposure time'
    medianFlatsHigh = allFlatsHigh[:,:,0]

medianFlatsLow = drpStella.where( medianFlatsLow, '==', 0., 0.000001, medianFlatsLow)
medianFlatsHigh = drpStella.where( medianFlatsHigh, '==', 0., 0.000001, medianFlatsHigh)
            
if medianFlatsOut != '':
    afwImage.makeImageFromArray(medianFlatsLow).writeFits(medianFlatsOut+"Low.fits")
    afwImage.makeImageFromArray(medianFlatsHigh).writeFits(medianFlatsOut+"High.fits")

divFlat = medianFlatsHigh / medianFlatsLow

if medianFlatsOut != '':
    afwImage.makeImageFromArray(divFlat).writeFits(medianFlatsOut+"DivHighByLow.fits")
sctrl = afwMath.StatisticsControl()
sctrl.setWeighted(False)

divFlatIm = afwImage.ImageF(divFlat)
divFlatMIm = afwImage.makeMaskedImage(divFlatIm)

mean = afwMath.makeStatistics(divFlatMIm, afwMath.MEANCLIP).getValue()
print 'mean = ',mean

stddev = afwMath.makeStatistics(divFlatMIm,afwMath.STDEVCLIP).getValue()
print 'stddev = ',stddev

if outFile != '':
    text_file = open(outFile, "w")
    text_file.write("#\n")
    text_file.write("# Defects file\n")
    text_file.write("#\n")
    text_file.write("# Convert to a fits table using getDefectFits.py;  this is done for you by running scons\n")
    text_file.write("#\n")
    text_file.write("#CCD x0    y0    width height\n")
    nBad = 0
    for iRow in np.arange(nPixCut,divFlat.shape[0]-nPixCut):
        for iCol in np.arange(nPixCut,2045):
            if np.absolute(divFlat[iRow, iCol] - mean) > 5. * stddev:
                nBad = nBad + 1
                text_file.write("%d %d %d 1 1\n" % (ccd, iRow, iCol))
                divFlatMIm.getMask().getArray()[iRow,iCol] = 1

    for iRow in np.arange(nPixCut,divFlat.shape[0]-nPixCut):
        for iCol in np.arange(2053,4096-nPixCut):
            if np.absolute(divFlat[iRow, iCol] - mean) > 5. * stddev:
                nBad = nBad + 1
                text_file.write("%d %d %d 1 1\n" % (ccd, iRow, iCol))
                divFlatMIm.getMask().getArray()[iRow, iCol] = 1

    text_file.write("# Dead amps")
    text_file.close()

    print nBad,' bad pixels found'

if display:
    divFlatExp = afwImage.makeExposure(divFlatMIm)
    medianFlatsLowExp = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(medianFlatsLow)))
    medianFlatsLowExp.getMaskedImage().getMask().getArray()[:,:] = divFlatExp.getMaskedImage().getMask().getArray()[:,:]

    display0 = afwDisplay.getDisplay()
    print 'display0 = ',display0
    display0.mtv(medianFlatsLowExp, title="parent")
    display0.setMaskTransparency(10)
    display0.setMaskPlaneColor("CROSSTALK", "orange")
