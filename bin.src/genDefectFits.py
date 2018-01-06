#!/usr/bin/env python

import collections
import os.path
import re

import numpy
import pyfits

import lsst.log as log
from lsst.obs.pfs import PfsMapper

Defect = collections.namedtuple('Defect', ['x0', 'y0', 'width', 'height'])

def genDefectFits(source, targetDir, verbose):
    logLevel = log.FATAL
    if verbose:
        logLevel = log.INFO
    for logger in ['daf.persistence.LogicalLocation',
                   'CameraMapper']:
        log.Log.getLogger(logger).setLevel(logLevel)

    mapper = PfsMapper(root=".")
    camera = mapper.camera

    logger = log.Log.getLogger("genDefectFits")
    if not verbose:
        logLevel = log.WARN
    logger.setLevel(logLevel)
    if verbose:
        logger.setLevel(logLevel)

    ccds = dict()
    for ccd in camera:
        ccdName = ccd.getName()
        ccds[ccdName] = ccd.getId()

    defects = dict()

    f = open(source, "r")
    for line in f:
        line = re.sub("\#.*", "", line).strip()
        if len(line) == 0:
            continue
        ccd, x0, y0, width, height = re.split("\s+", line)
        if ccd not in ccds:
            raise RuntimeError("Unrecognised ccd: %s" % ccd)
        if ccd not in defects:
            defects[ccd] = list()
        defects[ccd].append(Defect(x0=int(x0), y0=int(y0), width=int(width), height=int(height)))
    f.close()

    for ccd in ccds:
        # Make empty defect FITS file for CCDs with no defects
        if ccd not in defects:
            defects[ccd] = list()

        columns = list()
        for colName in Defect._fields:
            colData = numpy.array([d._asdict()[colName] for d in defects[ccd]])
            col = pyfits.Column(name=colName, format="I", array=colData)
            columns.append(col)

        cols = pyfits.ColDefs(columns)
        table = pyfits.BinTableHDU.from_columns(cols)

        table.header['NAME'] = ccd
        name = os.path.join(targetDir, "defects_%s.fits" % ccd)
        logger.info("Writing %d defects from CCD %s (%s) to %s" % (table.header['NAXIS2'], ccd, ccds[ccd], name))
        if os.path.exists(name):
            if args.force:
                os.unlink(name)
            else:
                logger.warn("File %s already exists; use --force to overwrite" % name)
                continue

        table.writeto(name)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("defectsFile", type=str, help="Text file containing list of defects")
    parser.add_argument("targetDir", type=str, nargs="?", help="Directory for generated fits files")
    parser.add_argument("-f", "--force", action="store_true", help="Force operations")
    parser.add_argument("-v", "--verbose", action="store_true", help="Be chattier")
    args = parser.parse_args()

    if not args.targetDir:
        args.targetDir = os.path.split(args.defectsFile)[0]

    genDefectFits(args.defectsFile, args.targetDir, args.verbose)
