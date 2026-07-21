#!/usr/bin/env python
"""Install PFS calibration products into a Gen3 butler.

The ``badRefPixels``, ``h4Linearity`` and ``defects`` calibrations live as data
files in the ``drp_pfs_data`` package (or a directory given with ``--data-root``;
see that option for the expected layout). This script reads them, ``butler.put``s
them into a RUN collection, and ``registry.certify``s them into a CALIBRATION
collection with a validity timespan, so that the pipeline selects the
appropriate version by observation date (rather than by EUPS setup).

Collections follow the DMTN-222 convention (https://dmtn-222.lsst.io), as the
pipeline-generated ``dark`` and ``bias`` do. Under ``{root}/{ticket}/{tag}/``,
with ``{iteration}`` a ``YYYYMMDDv`` date string:

- ``{product}.{iteration}`` -- CALIBRATION, certified into, and what the
  pipeline selects the calibration from;
- ``{product}Gen.{iteration}`` -- CHAINED, gathering this generation's RUNs;
- ``{product}Gen.{iteration}/{run}`` -- the RUN(s) holding the datasets. For
  ``defects`` each dated version gets its own run, named for its validity date;
  for ``badRefPixels`` and ``h4Linearity`` a single run named for the
  generation time.

Note the certified collection drops the ``Gen`` suffix, which marks the
generation chain. The root defaults to ``PFS/calib`` but may be overridden (e.g.
``u/<user>/<ticket>``) for testing.
"""

from __future__ import annotations

import datetime
import os
from argparse import ArgumentParser, RawTextHelpFormatter
from glob import glob

import astropy.time
import ruamel.yaml as yaml

from lsst.daf.butler import Butler, CollectionType, DatasetType, Timespan
from lsst.ip.isr import Defects
from lsst.utils import getPackageDir

from lsst.obs.pfs.nirBadRefPixels import NirBadRefPixels

NIR_CAMERAS = ("n1", "n2", "n3", "n4")

# Calibration dataset type definitions (storage class, dimensions) per product.
# These match PrimeFocusSpectrograph.register(); used by --register-dataset-types.
DATASET_TYPES = {
    "badRefPixels": ("NirBadRefPixels", ("instrument", "arm", "spectrograph")),
    "h4Linearity": ("H4Linearity", ("instrument", "arm", "spectrograph")),
    "defects": ("Defects", ("instrument", "detector", "arm", "spectrograph")),
}


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


def collectionRoot(root: str) -> str:
    """Normalize a collection base, stripping a trailing slash.

    Butler collection names are opaque strings, so ``PFS/calib/`` would otherwise
    yield ``PFS/calib//PIPE2D-...`` and be registered under that name.
    """
    root = root.rstrip("/")
    if not root:
        raise ValueError("collection root must not be empty")
    return root


def calibCollectionName(root: str, ticket: str, tag: str, product: str, iteration: str) -> str:
    """The CALIBRATION collection the pipeline selects the calibration from.

    ``{root}/{ticket}/{tag}/{product}.{iteration}`` -- no ``Gen`` suffix, which
    belongs to the generation chain.
    """
    return f"{collectionRoot(root)}/{ticket}/{tag}/{product}.{iteration}"


def genCollectionName(root: str, ticket: str, tag: str, product: str, iteration: str) -> str:
    """The CHAINED collection gathering everything one generation produced."""
    return f"{collectionRoot(root)}/{ticket}/{tag}/{product}Gen.{iteration}"


def runName(genCollection: str, suffix: str) -> str:
    """The RUN the datasets are written to, within the generation chain."""
    return f"{genCollection}/{suffix}"


def defaultTimestamp() -> str:
    """The default RUN suffix: the current UTC time, as ``YYYYMMDDTHHMMSSZ``."""
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def chainRun(butler, genCollection: str, run: str) -> None:
    """Register the generation CHAINED collection if needed and prepend ``run``.

    Members already in the chain keep their position; the new run is prepended so
    the most recent resolves first. Only this generation's own chain is touched;
    shared chains such as ``PFS/calib`` are maintained by hand.
    """
    registry = butler.registry
    registry.registerCollection(genCollection, CollectionType.CHAINED)
    chain = list(registry.getCollectionChain(genCollection))
    if run not in chain:
        registry.setCollectionChain(genCollection, [run, *chain])


def certifyGroup(butler, calibCollection, genCollection, run, datasetTypeName, entries):
    """Put a group of calibrations into ``run``, certify them, and chain the run.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Writeable butler.
    calibCollection : `str`
        CALIBRATION collection to certify into.
    genCollection : `str`
        Generation CHAINED collection the run is added to.
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
    chainRun(butler, genCollection, run)


def registerDatasetType(butler, product):
    """Register the calibration dataset type for a product (idempotent)."""
    storageClass, dimensions = DATASET_TYPES[product]
    datasetType = DatasetType(
        product, dimensions, storageClass, universe=butler.dimensions, isCalibration=True
    )
    inserted = butler.registry.registerDatasetType(datasetType)
    print(f"{product}: {'registered' if inserted else 'already registered'} dataset type")


def installBadRefPixels(butler, drpData, calibCollection, genCollection, timestamp,
                        instrument, begin, end):
    """Install the per-detector bad IRP reference-pixel calibrations.

    The multi-camera ``h4/badRefPixels.yaml`` source is split into one
    per-detector calibration. All detectors share one run, named for the
    generation time.
    """
    path = os.path.join(drpData, "h4", "badRefPixels.yaml")
    with open(path) as fd:
        config = yaml.YAML(typ="safe", pure=True).load(fd)

    run = runName(genCollection, timestamp)
    entries = []
    for cam in NIR_CAMERAS:
        if cam not in config:
            print(f"badRefPixels: {cam} not defined in {path}; skipping")
            continue
        calib = NirBadRefPixels.fromList(config[cam], cam)
        dataId = dict(instrument=instrument, arm="n", spectrograph=int(cam[1]))
        entries.append((calib, dataId, Timespan(begin, end)))
    certifyGroup(butler, calibCollection, genCollection, run, "badRefPixels", entries)
    print(f"badRefPixels: certified {len(entries)} detector(s) into {calibCollection}")


def installH4Linearity(butler, drpData, calibCollection, genCollection, timestamp,
                       instrument, begin, end):
    """Install the per-detector H4 nonlinearity corrections.

    All detectors share one run, named for the generation time.
    """
    from lsst.obs.pfs import h4Linearity

    run = runName(genCollection, timestamp)
    entries = []
    for cam in NIR_CAMERAS:
        path = os.path.join(drpData, "nirLinearity", f"nirLinearity-{cam}.fits")
        if not os.path.exists(path):
            print(f"h4Linearity: {path} not found; skipping {cam}")
            continue
        if not h4Linearity.isH4LinearityFile(path):
            print(f"h4Linearity: {path} is not an h4Linearity-format file; skipping {cam}")
            continue
        calib = h4Linearity.loadFits(path)
        dataId = dict(instrument=instrument, arm="n", spectrograph=int(cam[1]))
        entries.append((calib, dataId, Timespan(begin, end)))
    certifyGroup(butler, calibCollection, genCollection, run, "h4Linearity", entries)
    print(f"h4Linearity: certified {len(entries)} detector(s) into {calibCollection}")


def detectorIds(butler, instrument):
    """Return a mapping of detector full_name -> detector id."""
    records = butler.registry.queryDimensionRecords("detector", instrument=instrument)
    return {record.full_name: record.id for record in records}


def installDefects(butler, drpData, calibCollection, genCollection, instrument):
    """Install the per-detector defects, certifying the existing dated chain.

    Each detector directory holds one or more ISO-dated ``.ecsv`` files. The
    versions are certified with consecutive validity timespans: version ``i``
    is valid from its own date until the next version's date, and the final
    version is valid open-ended. Each version shares a dataId, so it gets its
    own run, named for its validity date.
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
            run = runName(genCollection, os.path.splitext(os.path.basename(path))[0])
            certifyGroup(
                butler, calibCollection, genCollection, run, "defects",
                [(calib, dataId, Timespan(begin, end))]
            )
            numCertified += 1
    print(f"defects: certified {numCertified} version(s) into {calibCollection}")


def makeParser():
    """Build the command-line parser.

    Uses ``RawTextHelpFormatter`` so the description and the multi-line
    ``--data-root`` help keep their layout instead of being reflowed into a
    single block.
    """
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "repo",
        help="Path to the butler repository.",
    )
    parser.add_argument(
        "--data-root",
        default=None,
        metavar="DIR",
        help=(
            "Directory to read the calibration data files from, instead of the\n"
            "installed drp_pfs_data package. Files are read from, under DIR:\n"
            "    badRefPixels : h4/badRefPixels.yaml\n"
            "    h4Linearity  : nirLinearity/nirLinearity-<cam>.fits   (cam = n1..n4)\n"
            "    defects      : curated/pfs/defects/<detector>/<date>.ecsv\n"
            "(default: the installed drp_pfs_data package)"
        ),
    )
    parser.add_argument(
        "--instrument",
        default="PFS",
        help="Instrument name (default: PFS).",
    )
    parser.add_argument(
        "--products",
        nargs="+",
        default=["badRefPixels", "h4Linearity", "defects"],
        choices=["badRefPixels", "h4Linearity", "defects"],
        help="Products to install (default: all).",
    )
    parser.add_argument(
        "--root",
        default="PFS/calib",
        help="Collection root (default: PFS/calib). Override (e.g. u/<user>/<ticket>) for testing.",
    )
    parser.add_argument(
        "--ticket",
        required=True,
        help="Ticket name, e.g. PIPE2D-1856.",
    )
    parser.add_argument(
        "--tag",
        required=True,
        help="Tag/label for this calibration set.",
    )
    parser.add_argument(
        "--register-dataset-types",
        action="store_true",
        help="Register the calibration dataset types before certifying (idempotent).",
    )
    parser.add_argument(
        "--iteration",
        default=None,
        help="DMTN-222 rerun iteration (YYYYMMDDv); default: today + 'a'.",
    )
    parser.add_argument(
        "--begin-date",
        default=None,
        help="Validity start (ISO-8601 TAI) for badRefPixels/h4Linearity; default: unbounded.",
    )
    parser.add_argument(
        "--end-date",
        default=None,
        help="Validity end (ISO-8601 TAI) for badRefPixels/h4Linearity; default: unbounded.",
    )
    return parser


def resolveDataRoot(args):
    """Return the directory holding the calibration data files: ``--data-root``
    if given, else the installed ``drp_pfs_data`` package directory."""
    if args.data_root is not None:
        return args.data_root
    return getPackageDir("drp_pfs_data")


def main():
    args = makeParser().parse_args()

    iteration = args.iteration
    if iteration is None:
        iteration = datetime.date.today().strftime("%Y%m%d") + "a"
    timestamp = defaultTimestamp()

    begin = isoTai(args.begin_date)
    end = isoTai(args.end_date)

    drpData = resolveDataRoot(args)
    butler = Butler(args.repo, writeable=True)

    for product in args.products:
        if args.register_dataset_types:
            registerDatasetType(butler, product)
        calibCollection = calibCollectionName(args.root, args.ticket, args.tag, product, iteration)
        genCollection = genCollectionName(args.root, args.ticket, args.tag, product, iteration)
        butler.registry.registerCollection(calibCollection, CollectionType.CALIBRATION)
        if product == "badRefPixels":
            installBadRefPixels(butler, drpData, calibCollection, genCollection, timestamp,
                                args.instrument, begin, end)
        elif product == "h4Linearity":
            installH4Linearity(butler, drpData, calibCollection, genCollection, timestamp,
                               args.instrument, begin, end)
        elif product == "defects":
            installDefects(butler, drpData, calibCollection, genCollection, args.instrument)


if __name__ == "__main__":
    main()
