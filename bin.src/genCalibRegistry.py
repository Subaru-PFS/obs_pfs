#!/usr/bin/env python
#/Users/azuri/stella-git/obs_subaru/bin.src/genCalibRegistry.py --root='/Volumes/My Passport/data/spectra/pfs/PFS/CALIB' --camera PFS --validity 180
import glob
import os
import re
try:
    import sqlite3 as sqlite
except ImportError:
    import sqlite
import sys
import math
import datetime
import collections

import lsst.daf.base   as dafBase
import lsst.afw.image  as afwImage

import optparse

try:
    import psycopg2 as pgsql
    havePgSql = True
except ImportError:
    try:
        from pg8000 import DBAPI as pgsql
        havePgSql = True
    except ImportError:
        havePgSql = False
if havePgSql:
    from lsst.daf.butlerUtils import PgSqlConfig

# Needs optparse w/ --create, etc. CPL
parser = optparse.OptionParser()
parser.add_option("--create", dest="create", default=False, action="store_true",
                  help="Create new registry (clobber old)?")
parser.add_option("--root", dest="root", default=".", help="Root directory")
parser.add_option("--camera", dest="camera", default="hsc", help="Camera name: HSC|SC|PFS")
parser.add_option("--validity", type="int", dest="validity", default=30, help="Calibration validity (days)")
opts, args = parser.parse_args()

if len(args) > 0 or len(sys.argv) == 1:
    print "Unrecognised arguments:", sys.argv[1:]
    parser.print_help()
    sys.exit(1)

if opts.camera.lower() not in ("suprime-cam", "suprimecam", "sc", "hsc", "hscsim", "pfs"):
    raise RuntimeError("Camera not recognised: %s" % camera)

if os.path.exists(os.path.join(opts.root, 'calibRegistry_pgsql.py')) and havePgSql:
    isSqlite = False
else:
    isSqlite = True

if isSqlite:
    registry = os.path.join(opts.root, "calibRegistry.sqlite3")

    if os.path.exists(registry):
        os.unlink(registry)
    conn = sqlite.connect(registry)
else:
    pgsqlConf = PgSqlConfig()
    pgsqlConf.load(os.path.join(opts.root, 'calibRegistry_pgsql.py'))
    conn = pgsql.connect(host=pgsqlConf.host, port=pgsqlConf.port,
                         user=pgsqlConf.user, password=pgsqlConf.password,
                         database=pgsqlConf.db)
    cur = conn.cursor()

if opts.camera.lower() in ("pfs"):
    Row = collections.namedtuple("Row", ["calibDate", "calibVersion", "spectrograph", "arm", "ccd"])
else:
    Row = collections.namedtuple("Row", ["calibDate", "calibVersion", "ccd"])

for calib in ('bias', 'dark', 'flat', 'fringe', 'traceDef'):
    if isSqlite:
        cmd = "create table " + calib.lower() + " (id integer primary key autoincrement"
        cmd += ", validStart text, validEnd text"
        if opts.camera.lower() in ("pfs"):
            cmd += ", calibDate text, filter text, calibVersion text, spectrograph int, arm text, ccd int"
        else:
            cmd += ", calibDate text, filter text, calibVersion text, ccd int"
        cmd += ")"
        conn.execute(cmd)
    else:
        cmd = "DROP TABLE IF EXISTS " + calib.lower()
        cur.execute(cmd)
        cmd = "create table " + calib.lower() + " (id SERIAL NOT NULL PRIMARY KEY"
        cmd += ", validStart VARCHAR(10), validEnd VARCHAR(10)"
        if opts.camera.lower() in ("pfs"):
            cmd += ", calibDate VARCHAR(10), filter VARCHAR(5), calibVersion VARCHAR(5), spectrogaph INT, arm VARCHAR(1), ccd INT"
        else:
            cmd += ", calibDate VARCHAR(10), filter VARCHAR(16), calibVersion VARCHAR(16), ccd INT"
        cmd += ")"
        cur.execute(cmd)
    conn.commit()

    rowsPerFilter = dict()

    if opts.camera.lower() in ("pfs"):
        checkFits = os.path.join(opts.root, calib.upper(), "*", "*",
                                       "*-20*-*-*-*-*.fits*")
    else:
        checkFits = os.path.join(opts.root, calib.upper(), "20*-*-*", "*", "*",
                             calib.upper() + "-*.fits*")
    print "checkFits = <",checkFits,">"

    for fits in glob.glob(checkFits):
        print "fits = <"+fits+">"
        if opts.camera.lower() in ("suprime-cam", "suprimecam", "sc"):
            m = re.search(r'\w+/(\d{4})-(\d{2})-(\d{2})/([\w+-]+)/([\w-]+)/\w+-(\d).fits', fits)
        elif opts.camera.lower() in ("hsc", "hscsim"):
            m = re.search(r'\w+/(\d{4})-(\d{2})-(\d{2})/([\w+-]+)/([\w-]+)/\w+-(\d{3}).fits', fits)#"BIAS/%(calibDate)s/NONE/%(calibVersion)s/BIAS-%(ccd)03d.fits"
        elif opts.camera.lower() in ("pfs"):
            print 'fits = <',fits,'>'
            m = re.search(r'\w+/([\w-]+)/([\w-]+)/pfs\w+-(\d{4})-(\d{2})-(\d{2})-0-(\d)(\w).fits', fits)#pfsBias-007251-2m.fits
            print 'm = ',m
        if not m:
            if (opts.camera.lower() in ("suprime-cam", "suprimecam", "sc", "pfs") and
                re.search(r'.*/\w+-0000000(\d).fits', fits)):
                print >>sys.stderr, ("Warning: file naming rules have changed. " +
                                     "Try renaming %s to %s" %
                                     (fits, re.sub(r"(.*/\w+)-0000000(\d).fits", r"\1-\2.fits", fits)))
            print >>sys.stderr, "Warning: Unrecognized file name:", fits
            continue

        print "Registering:", fits
        print m.groups()
        if opts.camera.lower() in ("pfs"):
            filterName, version, year, month, day, spectrograph, arm = m.groups()
            print 'year = ',year,', month = ',month,', day = ',day,', filterName = ',filterName,', version = ',version,', spectrograph = ',spectrograph,', arm = ',arm
            ccd = int(spectrograph) - 1
            if arm in ("m"):
                ccd = ccd + 4
            elif arm in ("r"):
                ccd += 4
            elif arm in ("n"):
                ccd += 8
        else:
            year, month, day, filterName, version, ccd = m.groups()

        date = datetime.date(int(year), int(month), int(day))
        if filterName not in rowsPerFilter:
            rowsPerFilter[filterName] = list()
        if opts.camera.lower() in ("pfs"):
            rowsPerFilter[filterName].append(Row(date, version, spectrograph, arm, ccd))
        else:
            rowsPerFilter[filterName].append(Row(date, version, ccd))


    # Fix up the validStart,validEnd so there are no overlaps
    for filterName, rows in rowsPerFilter.items():
        rows.sort(key=lambda row: row.calibDate)
        validity = datetime.timedelta(opts.validity)
        valids = collections.OrderedDict([(row.calibDate, [row.calibDate - validity,
                                                           row.calibDate + validity]) for row in rows])
        dates = valids.keys()
        numDates = len(dates)
        midpoints = [ t1 + (t2 - t1)//2 for t1, t2 in zip(dates[:numDates-1], dates[1:]) ]
        for i, (date, midpoint) in enumerate(zip(dates[:numDates-1], midpoints)):
            if valids[date][1] > midpoint:
                nextDate = dates[i+1]
                #print "Adjusting: %d %s --> %s : %s vs %s" % (i, date, nextDate, valids[date][1], midpoint)
                valids[nextDate][0] = midpoint
                valids[date][1] = midpoint
        del midpoints
        del dates

        for row in rows:
            calibDate = row.calibDate.isoformat()
            validStart = valids[row.calibDate][0].isoformat()
            validEnd = valids[row.calibDate][1].isoformat()

            # print "%f --> %f %f" % (calibDate, validStart, validEnd)

            if isSqlite:
                if opts.camera.lower() in ("pfs"):
                    conn.execute("INSERT INTO " + calib.lower() + " VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?)",
                                (validStart, validEnd, calibDate, filterName, row.calibVersion, row.spectrograph, row.arm, row.ccd))
                else:
                    conn.execute("INSERT INTO " + calib.lower() + " VALUES (NULL, ?, ?, ?, ?, ?, ?)",
                                 (validStart, validEnd, calibDate, filterName, row.calibVersion, row.ccd))
            else:
                if opts.camera.lower() in ("pfs"):
                    cur.execute("INSERT INTO " + calib.lower() + " (validStart, validEnd, calibDate, filter, calibVersion, spectrograph, arm, ccd) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)",
                                (validStart, validEnd, calibDate, filterName, row.calibVersion, row.spectrograph, row.arm, row.ccd))
                else:
                    cur.execute("INSERT INTO " + calib.lower() + " (validStart, validEnd, calibDate, filter, calibVersion, ccd) VALUES (%s, %s, %s, %s, %s, %s)",
                                (validStart, validEnd, calibDate, filterName, row.calibVersion, row.ccd))


conn.commit()
if not isSqlite:
    cur.close()
conn.close()
