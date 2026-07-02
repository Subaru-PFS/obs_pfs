#!/usr/bin/env python
"""Install PFS calibration products into a Gen3 butler.

The ``nirLinearity``, ``badRefPixels`` and ``defects`` calibrations live as data
files in the ``drp_pfs_data`` package. This script reads them, ``butler.put``s
them into a RUN collection, and ``registry.certify``s them into a
CALIBRATION collection with a validity timespan, so that the pipeline selects
the appropriate version by observation date (rather than by EUPS setup).

Collections follow the DMTN-222 convention
``{root}/{ticket}/{tag}/{product}Gen.{iteration}``, where ``{iteration}`` is a
``YYYYMMDDv`` date string (creation date plus a trailing letter, incremented if
a generation is retried). The root defaults to ``PFS/calib`` but may be
overridden (e.g. ``u/<user>/<ticket>``) for testing.
"""

from __future__ import annotations

import datetime
import os
from argparse import ArgumentParser
from glob import glob

import astropy.time
import ruamel.yaml as yaml

from lsst.daf.butler import Butler, CollectionType, Timespan
from lsst.ip.isr import Defects
from lsst.utils import getPackageDir

from lsst.obs.pfs.nirBadRefPixels import NirBadRefPixels
from lsst.obs.pfs.nirLinearity import NirLinearity

NIR_CAMERAS = ("n1", "n2", "n3", "n4")


def isoTai(value: str | None) -> astropy.time.Time | None:
    """Parse an ISO-8601 TAI string into a `~astropy.time.Time` (or `None`)."""
    if value is None:
        return None
    return astropy.time.Time(value, format="isot", scale="tai")


def dateFromFilename(path: str) -> astropy.time.Time:
    """Parse the validity-start date encoded in a dated calibration filename.

    The stem is an ISO-8601 basic-format datetime, e.g. ``20240815T000000``.
    Interpreted as UTC, matching the historical curated-calibration behaviour.
    """
    stem = os.path.splitext(os.path.basename(path))[0]
    when = datetime.datetime.fromisoformat(stem)
    return astropy.time.Time(when, format="datetime", scale="utc")


def collectionName(root: str, ticket: str, tag: str, product: str, iteration: str) -> str:
    """Compose the DMTN-222 CALIBRATION collection name for a product."""
    return f"{root}/{ticket}/{tag}/{product}Gen.{iteration}"


def certifyGroup(butler, calibCollection, run, datasetTypeName, entries):
    """Put and certify a group of calibrations.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Writeable butler.
    calibCollection : `str`
        CALIBRATION collection to certify into.
    run : `str`
        RUN collection to put the datasets into.
    datasetTypeName : `str`
        Dataset type name.
    entries : iterable of (calib, dataId, `~lsst.daf.butler.Timespan`)
        The calibration objects, their data IDs, and validity timespans.
    """
    butler.registry.registerRun(run)
    with butler.transaction():
        for calib, dataId, timespan in entries:
            ref = butler.put(calib, datasetTypeName, dataId, run=run)
            butler.registry.certify(calibCollection, [ref], timespan)


def installNirLinearity(butler, drpData, calibCollection, instrument, begin, end):
    """Install the per-detector NIR linearity calibrations."""
    run = f"{calibCollection}/put"
    entries = []
    for cam in NIR_CAMERAS:
        path = os.path.join(drpData, "nirLinearity", f"nirLinearity-{cam}.fits")
        if not os.path.exists(path):
            print(f"nirLinearity: {path} not found; skipping {cam}")
            continue
        calib = NirLinearity.readFits(path)
        dataId = dict(instrument=instrument, arm="n", spectrograph=int(cam[1]))
        entries.append((calib, dataId, Timespan(begin, end)))
    certifyGroup(butler, calibCollection, run, "nirLinearity", entries)
    print(f"nirLinearity: certified {len(entries)} detector(s) into {calibCollection}")


def installBadRefPixels(butler, drpData, calibCollection, instrument, begin, end):
    """Install the per-detector bad IRP reference-pixel calibrations.

    The multi-camera ``h4/badRefPixels.yaml`` source is split into one
    per-detector calibration.
    """
    path = os.path.join(drpData, "h4", "badRefPixels.yaml")
    with open(path) as fd:
        config = yaml.YAML(typ="safe", pure=True).load(fd)

    run = f"{calibCollection}/put"
    entries = []
    for cam in NIR_CAMERAS:
        if cam not in config:
            print(f"badRefPixels: {cam} not defined in {path}; skipping")
            continue
        calib = NirBadRefPixels.fromList(config[cam], cam)
        dataId = dict(instrument=instrument, arm="n", spectrograph=int(cam[1]))
        entries.append((calib, dataId, Timespan(begin, end)))
    certifyGroup(butler, calibCollection, run, "badRefPixels", entries)
    print(f"badRefPixels: certified {len(entries)} detector(s) into {calibCollection}")


def detectorIds(butler, instrument):
    """Return a mapping of detector full_name -> detector id."""
    records = butler.registry.queryDimensionRecords("detector", instrument=instrument)
    return {record.full_name: record.id for record in records}


def installDefects(butler, drpData, calibCollection, instrument):
    """Install the per-detector defects, certifying the existing dated chain.

    Each detector directory holds one or more ISO-dated ``.ecsv`` files. The
    versions are certified with consecutive validity timespans: version ``i``
    is valid from its own date until the next version's date, and the final
    version is valid open-ended.
    """
    nameToId = detectorIds(butler, instrument)
    defectsDir = os.path.join(drpData, "curated", "pfs", "defects")

    numCertified = 0
    for detectorName in sorted(nameToId):
        detectorDir = os.path.join(defectsDir, detectorName)
        paths = sorted(glob(os.path.join(detectorDir, "*.ecsv")))
        if not paths:
            continue
        times = [dateFromFilename(path) for path in paths]
        ends = times[1:] + [None]

        arm = detectorName[0]
        spectrograph = int(detectorName[1])
        detector = nameToId[detectorName]

        for path, begin, end in zip(paths, times, ends):
            calib = Defects.readText(path)
            dataId = dict(
                instrument=instrument, detector=detector, arm=arm, spectrograph=spectrograph
            )
            # Each version shares a dataId, so it needs its own RUN.
            run = f"{calibCollection}/{os.path.splitext(os.path.basename(path))[0]}"
            certifyGroup(
                butler, calibCollection, run, "defects", [(calib, dataId, Timespan(begin, end))]
            )
            numCertified += 1
    print(f"defects: certified {numCertified} version(s) into {calibCollection}")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("repo", help="Path to the butler repository")
    parser.add_argument("--instrument", default="PFS", help="Instrument name (default: PFS)")
    parser.add_argument(
        "--products",
        nargs="+",
        default=["nirLinearity", "badRefPixels", "defects"],
        choices=["nirLinearity", "badRefPixels", "defects"],
        help="Products to install (default: all)",
    )
    parser.add_argument("--root", default="PFS/calib", help="Collection root (default: PFS/calib)")
    parser.add_argument("--ticket", required=True, help="Ticket name, e.g. PIPE2D-1856")
    parser.add_argument("--tag", required=True, help="Tag/label for this calibration set")
    parser.add_argument(
        "--iteration",
        default=None,
        help="DMTN-222 rerun iteration (YYYYMMDDv); default: today + 'a'",
    )
    parser.add_argument(
        "--begin-date",
        default=None,
        help="Validity start (ISO-8601 TAI) for nirLinearity/badRefPixels; default: unbounded",
    )
    parser.add_argument(
        "--end-date",
        default=None,
        help="Validity end (ISO-8601 TAI) for nirLinearity/badRefPixels; default: unbounded",
    )
    args = parser.parse_args()

    iteration = args.iteration
    if iteration is None:
        iteration = datetime.date.today().strftime("%Y%m%d") + "a"

    begin = isoTai(args.begin_date)
    end = isoTai(args.end_date)

    drpData = getPackageDir("drp_pfs_data")
    butler = Butler(args.repo, writeable=True)

    for product in args.products:
        calibCollection = collectionName(args.root, args.ticket, args.tag, product, iteration)
        butler.registry.registerCollection(calibCollection, CollectionType.CALIBRATION)
        if product == "nirLinearity":
            installNirLinearity(butler, drpData, calibCollection, args.instrument, begin, end)
        elif product == "badRefPixels":
            installBadRefPixels(butler, drpData, calibCollection, args.instrument, begin, end)
        elif product == "defects":
            installDefects(butler, drpData, calibCollection, args.instrument)


if __name__ == "__main__":
    main()
