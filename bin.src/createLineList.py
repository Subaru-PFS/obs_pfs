#!/usr/bin/env python

import sys
import os.path
import re

import lsst.utils
import numpy as np
import pyfits
import collections

Line = collections.namedtuple('Line', ['wavelength', 'element', 'strength'])
distStrength = collections.namedtuple('distStrength', ['distance', 'strengthRatio', 'lowStrength'])

def createLineList(sources, minDist, minStrength, minRatio, minLambda, maxLambda):
    lines = []
    dir = lsst.utils.getPackageDir("obs_pfs")
    subDir = 'pfs/lineLists'
    for el in sources:
        inFile = os.path.join(dir,subDir,"%s_raw.txt" % (el))
        f = open(inFile, "r")
        for line in f:
            if len(line) == 0 or line == "\n":
                continue
            match = re.match("(\w*)\s(I*)\s*\S\s*(\S*)\s*\S\s*\S*\s*\S\s*(\S*)", line)
            if not match:
                print 'error matching'
            element, ion, lam, strength = match.groups()
            element = "%s%s" % (element, ion)
            if '|' not in lam:
                flag = ""
                match = re.match("(\d*)([Echw])", strength)
                if not match:
                    match = re.match("(\d*)", strength)
                    if not match:
                        print "error matching strength"
                    strength = match.groups()[0]
                else:
                    strength, flag = match.groups()
                    print 'flag = <',flag,'>, strength = ',strength
                if flag not in ['c','h','w','E']:
                    print 'flag = <',flag,'> not in [chwE]'
                    lines.append(Line(wavelength=float(lam)*10., element=element, strength=strength))
        f.close()

    lines = sorted(lines, key=lambda x: x.wavelength)
    goodLines = []

    for j in range(len(lines)):
        """For each line calculate distance and strength ratio to all other lines"""
        strengthA = int(lines[j].strength)
        wavelengthA = lines[j].wavelength

        distances = []
        for i in range(len(lines)):
            if i != j:
                match = re.match("(\d*)", lines[i].strength)
                if not match:
                    print "error matching strength"
                strengthB = match.groups()
                strengthB = int(strengthB[0])

                ratio = strengthA / strengthB
                lowStrength = False
                if strengthA < strengthB:
                    ratio = strengthB / strengthA
                    lowStrength = True
                distances.append(distStrength(distance=abs(lines[i].wavelength - wavelengthA), strengthRatio=ratio, lowStrength=lowStrength))
        if strengthA >= minStrength:
            passed = True
            for i in np.arange(len(distances)-1,-1,-1):
                if distances[i].distance < minDist and (distances[i].strengthRatio < minRatio or
                                                        (distances[i].strengthRatio >= minRatio and distances[i].lowStrength)):
                    print 'problem: distances[',i,'].distance = ',distances[i].distance,', distances[',i,'].strengthRatio = ',distances[i].strengthRatio
                    print 'lines[',i,'] = ',lines[i]
                    print 'lines[',j,'] = ',lines[j]
                    passed = False
            if passed:
                if (lines[j].wavelength >= (minLambda*10.)) and (lines[j].wavelength <= (maxLambda * 10.)):
                    goodLines.append(j)
                else:
                    print 'rejecting line ',j,' because its wavelength=', lines[j].wavelength,' is outside [',minLambda*10.,', ',maxLambda*10.,']'
        else:
            print 'rejecting line ',j,' because its predicted strength is too low'
    print 'goodLines = ',goodLines

    # Write output FITS file
    columns = list()
    colData = np.array([lines[d].wavelength for d in goodLines])
    print 'wavelength colData = ',colData
    col = pyfits.Column(name='wavelength', format="E", array=colData)
    columns.append(col)

    colData = np.array([lines[d].element for d in goodLines])
    print 'element colData = ',colData
    col = pyfits.Column(name='element', format="2A", array=colData)
    columns.append(col)

    colData = np.array([lines[d].strength for d in goodLines])
    print 'element colData = ',colData
    col = pyfits.Column(name='strength', format="128A", array=colData)
    columns.append(col)

    cols = pyfits.ColDefs(columns)
    table = pyfits.BinTableHDU.from_columns(cols)

    name = "".join(sources)
    table.header['NAME'] = name
    name = os.path.join(dir,subDir, "%s.fits" % name)
    if os.path.exists(name):
        if args.force:
            os.unlink(name)
        else:
            print >> sys.stderr, "File %s already exists; use --force to overwrite" % name

    table.writeto(name)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("element", nargs='+', type=str, help="Element to create line list for (e.g. Hg)")
    parser.add_argument("--maxLambda", type=float, help="maximum wavelength of lines to be used in nm", default=975.)
    parser.add_argument("--minLambda", type=float, help="minimum wavelength of lines to be used in nm", default=622.)
    parser.add_argument("--minDist", type=float, help="minimum distance between lines", default=1.)
    parser.add_argument("--minStrength", type=int, help="minimum strength of lines", default=200)
    parser.add_argument("--minRatio", type=float, help="minimum strength ratio between adjacent lines", default=100.)
    parser.add_argument("-f", "--force", action="store_true", help="Force operations")
    parser.add_argument("-v", "--verbose", action="store_true", help="Be chattier")
    args = parser.parse_args()

    createLineList(args.element, args.minDist, args.minStrength, args.minRatio, args.minLambda, args.maxLambda)
