import lsst.daf.persistence as dafPersist
import lsst.afw.image as afwImage
import numpy as np
import os

homeDir = '/Volumes/My Passport/Users/azuri/spectra/pfs'
outFile = "/Users/azuri/stella-git/obs_pfs/pfs/defects/2016-07-14/defects.dat"

site = 'S'
category = 'A'
filter = 'm'
ccd = 5

visitLow = 6198
visitHigh = 6235
expTimeLow = 2
expTimeHigh = 64
nPixCut = 40
nCols = 4096
gapLow = 2045
gapHigh = 2052

medianFlatsOut = '/Users/azuri/spectra/pfs/medianFlat'#+'.fits' - leave empty if you don't want to write the median Flats

butler = dafPersist.Butler(os.path.join(homeDir,"PFS"))

flatsVisit = list()
flatsExpTimes = list()
flatsLow = list()
flatsHigh = list()
for visit in np.arange(visitLow,visitHigh+1):
    dataId = dict(site=site, category=category, filter=filter, visit=visit)
    fileName = butler.get('raw_filename', dataId)
    print 'fileName = <',fileName[0],'>'
    print ' '
    md = afwImage.readMetadata(fileName[0])

    print md.get('EXPTIME'),md.get('IMAGETYP')
    print ' '
    if md.get("IMAGETYP") == "flat":
        print 'flat found'
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

dataId = dict(site=site, category=category, filter=filter, visit=flatsLow[0], ccd=ccd)
image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
print image.shape
allFlatsLow = np.ndarray(shape=[image.shape[0], image.shape[1], len(flatsLow)], dtype='float32')
allFlatsHigh = np.ndarray(shape=[image.shape[0], image.shape[1], len(flatsHigh)], dtype='float32')
for i in np.arange(len(flatsLow)):
    dataId = dict(site=site, category=category, filter=filter, visit=flatsLow[i], ccd=5)
    image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
    allFlatsLow[:,:,i] = image
for i in np.arange(len(flatsHigh)):
    dataId = dict(site=site, category=category, filter=filter, visit=flatsHigh[i], ccd=5)
    image = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
    allFlatsHigh[:,:,i] = image
    
if len(flatsLow) > 2:
    medianFlatsLow = np.median(allFlatsLow,2)
elif len(flatsLow) == 2:
    medianFlatsLow = np.min(allFlatsLow,2)
else:
    print 'only one Flat with ',expTimeLow,'s exposure time'
if len(flatsHigh) > 2:
    medianFlatsHigh = np.median(allFlatsHigh,2)
elif len(flatsHigh) == 2:
    medianFlatsHigh = np.min(allFlatsHigh,2)
else:
    print 'only one Flat with ',expTimeHigh,'s exposure time'
    
for i in range(medianFlatsLow.shape[0]):
    for j in range(medianFlatsLow.shape[1]):
        if medianFlatsLow[i,j] == 0.:
            print 'medianFlatsLow[',i,',',j,'] = ',medianFlatsLow[i,j]
            medianFlatsLow[i,j] = 0.000001
for i in range(medianFlatsHigh.shape[0]):
    for j in range(medianFlatsHigh.shape[1]):
        if medianFlatsHigh[i,j] == 0.:
            print 'medianFlatsHigh[',i,',',j,'] = ',medianFlatsHigh[i,j]
            medianFlatsHigh[i,j] = 0.000001
            
if medianFlatsOut != '':
    afwImage.makeImageFromArray(medianFlatsLow).writeFits(medianFlatsOut+"Low.fits")
    afwImage.makeImageFromArray(medianFlatsHigh).writeFits(medianFlatsOut+"High.fits")

resFlat = medianFlatsHigh / medianFlatsLow
print 'resFlat = ',resFlat
afwImage.makeImageFromArray(resFlat).writeFits(medianFlatsOut+"Res.fits")

medianResFlat = np.median(resFlat[nPixCut:resFlat.shape[0]-nPixCut,nPixCut:resFlat.shape[1]-nPixCut])
print 'medianResFlat = ',medianResFlat.shape,': ',medianResFlat

stddev = np.std(resFlat[nPixCut:resFlat.shape[0]-nPixCut,nPixCut:resFlat.shape[1]-nPixCut])
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
    
    fiveStdDev = 5. * stddev
    # Break from Column 2045 to Column 2052 to avoid columns where physical CCDs meet
    for iRow in np.arange(nPixCut,resFlat.shape[0]-nPixCut):
        for iCol in np.arange(nPixCut,gapLow):#resFlat.shape[1]):
            if np.absolute(resFlat[iRow, iCol] - medianResFlat) > fiveStdDev:
                print 'bad pixel found at [',iRow,', ',iCol,']: ',resFlat[iRow,iCol]
                nBad = nBad + 1
                text_file.write("5 %d %d 1 1\n" % (iRow, iCol))

    for iRow in np.arange(nPixCut, resFlat.shape[0] - nPixCut):
        for iCol in np.arange(gapHigh + 1,nCols - nPixCut):#resFlat.shape[1]):
            if np.absolute(resFlat[iRow, iCol] - medianResFlat) > fiveStdDev:
                print 'bad pixel found at [',iRow,', ',iCol,']: ',resFlat[iRow,iCol]
                nBad = nBad + 1
                text_file.write("%d %d %d 1 1\n" % (ccd, iRow, iCol))

    text_file.write("# Dead amps")
    text_file.close()

    print nBad,' bad pixels found'