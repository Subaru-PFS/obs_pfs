#!/usr/bin/env python
import glob
import os
import re
import sqlite3 as sqlite
import sys
import datetime
import collections
import argparse

def main(root, validityDays):

    registry = os.path.join(root, "calibRegistry.sqlite3")

    if os.path.exists(registry):
        os.unlink(registry)
    conn = sqlite.connect(registry)

    Row = collections.namedtuple("Row", ["calibDate", "calibVersion", "spectrograph", "arm", "ccd"])

    for calib in ('bias', 'dark', 'flat', 'fiberFlat', 'imageFlat', 'fiberTrace'):
        cmd = "create table " + calib.lower() + " (id integer primary key autoincrement"
        cmd += ", validStart text, validEnd text"
        cmd += ", calibDate text, filter text, calibVersion text, spectrograph int, arm text, ccd int"
        cmd += ")"
        conn.execute(cmd)
        conn.commit()

        rowsPerFilter = dict()

        checkFits = os.path.join(root, calib.upper(), "*", "*", "*.fits")

        for fits in glob.glob(checkFits):
            m = re.search(r'\w+/([\w-]+)/([\w-]+)/pfs\w+-(\d{4})-(\d{2})-(\d{2})-0-(\d)(\w).fits', fits)
            if not m:
                m = re.search(r'\w+/([\w-]+)/(\d{4})-(\d{2})-(\d{2})/pfs(\w+)-(\d{6})-(\d)(\w).fits', fits)
                if not m:
                    print >>sys.stderr, "Warning: Unrecognized file name:", fits
                    raise Exception("Unrecognized file name")
                else:
                    filterName, year, month, day, version, visit, spectrograph, arm = m.groups()
                    version = version[0].lower() + version[1:]
            else:
                filterName, version, year, month, day, spectrograph, arm = m.groups()

            print "Registering:", fits
            date = datetime.date(int(year), int(month), int(day))
            ccd = int(spectrograph) - 1
            if arm in ("m"):
                ccd = ccd + 4
            elif arm in ("r"):
                ccd += 4
            elif arm in ("n"):
                ccd += 8

            if filterName not in rowsPerFilter:
                rowsPerFilter[filterName] = list()
            rowsPerFilter[filterName].append(Row(date, version, spectrograph, arm, ccd))

        # Fix up the validStart,validEnd so there are no overlaps
        for filterName, rows in rowsPerFilter.items():
            rows.sort(key=lambda row: row.calibDate)
            validity = datetime.timedelta(validityDays)
            valids = collections.OrderedDict([(row.calibDate, [row.calibDate - validity,
                                                               row.calibDate + validity]) for row in rows])
            dates = valids.keys()
            numDates = len(dates)
            midpoints = [ t1 + (t2 - t1)//2 for t1, t2 in zip(dates[:numDates-1], dates[1:]) ]
            for i, (date, midpoint) in enumerate(zip(dates[:numDates-1], midpoints)):
                if valids[date][1] > midpoint:
                    nextDate = dates[i+1]
                    valids[nextDate][0] = midpoint
                    valids[date][1] = midpoint
            del midpoints
            del dates

            for row in rows:
                calibDate = row.calibDate.isoformat()
                validStart = valids[row.calibDate][0].isoformat()
                validEnd = valids[row.calibDate][1].isoformat()

                conn.execute("INSERT INTO " + calib.lower() + " VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?)",
                             (validStart, validEnd, calibDate, filterName, row.calibVersion, row.spectrograph,
                              row.arm, row.ccd))

    conn.commit()
    conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", required=True, help="Root directory")
    parser.add_argument("--validity", type=int, dest="validity", default=30,
                        help="Calibration validity (days)")
    args = parser.parse_args()

    main(args.root, args.validity)

