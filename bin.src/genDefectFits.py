#!/usr/bin/env python

import sys
import os.path
import re

import numpy
import pyfits
import collections

import lsst.pex.policy as pexPolicy
from lsst.obs.pfs import PfsMapper
import lsst.afw.cameraGeom.utils as afwCGU

Defect = collections.namedtuple('Defect', ['x0', 'y0', 'width', 'height'])
mapperMap = {'pfs': PfsMapper}

def genDefectFits(cameraName, source, targetDir):
    mapper = mapperMap[cameraName.lower()](root=".")
    camera = mapper.camera

    ccds = dict()
    for ccd in camera:
        print 'dir(ccd) = ',dir(ccd)
        ccdNum = ccd.getId()
        print 'ccdNum = ',ccdNum
        ccds[ccdNum] = ccd.getName()
        print 'ccds[',ccdNum,'] = ',ccds[ccdNum]

    defects = dict()

    f = open(source, "r")
    for line in f:
        line = re.sub("\#.*", "", line).strip()
        if len(line) == 0:
            continue
        ccd, x0, y0, width, height = re.split("\s+", line)
        ccd = int(ccd)
        if not ccds.has_key(ccd):
            raise RuntimeError("Unrecognised ccd: %d" % ccd)
        if not defects.has_key(ccd):
            defects[ccd] = list()
        defects[ccd].append(Defect(x0=int(x0), y0=int(y0), width=int(width), height=int(height)))
    f.close()

    for ccd in ccds:
        # Make empty defect FITS file for CCDs with no defects
        if not defects.has_key(ccd):
            defects[ccd] = list()

        columns = list()
        for colName in Defect._fields:
            colData = numpy.array([d._asdict()[colName] for d in defects[ccd]])
            col = pyfits.Column(name=colName, format="I", array=colData)
            columns.append(col)

        cols = pyfits.ColDefs(columns)
        table = pyfits.new_table(cols)

        table.header['NAME'] = ccd
        print 'ccd = ',ccd
        print 'ccds[ccd] = ',ccds[ccd]
#        table.header['NAME'] = ccds[ccd]

#        print 'defects_%s.fits = <defects_',ccds[ccd],'.fits'
        print 'defects_%s.fits = <defects_',ccd,'.fits'
#        name = os.path.join(targetDir, "defects_%s.fits" % ccds[ccd])
        name = os.path.join(targetDir, "defects_%s.fits" % ccd)
        print "Writing %d defects from CCD %d (%s) to %s" % (table.header['NAXIS2'], ccd, ccds[ccd], name)
        if os.path.exists(name):
            if args.force:
                os.unlink(name)
            else:
                print >> sys.stderr, "File %s already exists; use --force to overwrite" % name
                continue

        table.writeto(name)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("cameraName", type=str, choices=['PFS'],
                        help="Camera name: PFS")
    parser.add_argument("defectsFile", type=str, help="Text file containing list of defects")
    parser.add_argument("targetDir", type=str, nargs="?", help="Directory for generated fits files")
    parser.add_argument("-f", "--force", action="store_true", help="Force operations")
    parser.add_argument("-v", "--verbose", action="store_true", help="Be chattier")
    args = parser.parse_args()

    if not args.targetDir:
        args.targetDir = os.path.split(args.defectsFile)[0]

    genDefectFits(args.cameraName, args.defectsFile, args.targetDir)
