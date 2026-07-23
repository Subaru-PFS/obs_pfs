import datetime
import logging
import os
import time
from typing import NamedTuple

import astropy.time
import numpy as np

import fitsio

import lsst.daf.butler as dafButler
from lsst.daf.butler import CollectionType, DatasetType, MissingCollectionError, Timespan
from lsst.obs.pfs import imageCube

from pfs.drp.stella.calibs import setCalibHeader

logger = logging.getLogger(__name__)

# Dimensions and storage class of the nirDark family of dataset types, matching
# the registrations in `PrimeFocusSpectrograph.registerDatasetTypes`.
DARK_DIMENSIONS = ("instrument", "arm", "spectrograph")
DARK_STORAGE_CLASS = "ImageCube"


def getExpInfo(exp):
    """Get exposure start, read noise, and gain for an exposure.

    Returns
    -------
    startTime: `datetime.datetime`
    gain: `float`
    readNoise: `float`
    """
    startTime: datetime.datetime = exp.getInfo().getVisitInfo().getDate().toPython()
    det = exp.getDetector()
    amp = det.getAmplifiers()[0]
    gain = amp.getGain()
    readNoise = amp.getReadNoise()
    return startTime, gain, readNoise


def getStartDate(butler: dafButler.Butler, dataId: dict, visit: int) -> datetime.datetime:
    """The observation date of a visit, for use as a calibration's validity start."""
    return getExpInfo(butler.get('postISRCCD', dataId, visit=visit))[0]


def datasetTypeForIrp(irpN: int) -> str:
    """Return the nirDark dataset type for a given IRP ratio.

    The dark is subtracted from an exposure ramp read-by-read, so it is only
    valid for exposures with the same per-read cadence. The IRP ratio changes
    that cadence, so each ratio has its own product: the default (IRP1) writes
    ``nirDark``, and other ratios write ``nirDark_irp<N>`` (e.g. ``nirDark_irp4``
    for IRP4).
    """
    return "nirDark" if irpN == 1 else f"nirDark_irp{irpN}"


def rampInfoFromCubes(cubes: dict, vlist: list[int]) -> tuple[dict, int]:
    """Validate a set of open input ramps and extract their common IRP/timing cards.

    The primary header of a rawISRCube carries the raw H4 metadata (see
    `PfsIsrTask._makeExposure`).

    Parameters
    ----------
    cubes : `dict` mapping visit to `fitsio.FITS`
        The open input ramps.
    vlist : list[int]
        The visits to check, in order.

    Returns
    -------
    rampInfo : `dict`
        ``W_H4IRPN``, ``W_H4IRPO``, ``W_H4NRED``, ``W_H4FRMT``, common to all inputs.
    nreads : `int`
        The number of reads per ramp.

    Raises
    ------
    ValueError
        If the input ramps do not all share the same number of reads or the same
        IRP ratio: they cannot be meaningfully combined.
    """
    rampInfo = None
    nreads = None
    for v in vlist:
        hdr = cubes[v][0].read_header()
        vReads = len(cubes[v]) - 1
        # W_H4IRPO is the IRP pixel's position within each block of IRPN
        # data pixels: it is 1-based and lies in 1..W_H4IRPN, so 1 (not 0)
        # is the fallback when the card is absent.
        vInfo = dict(W_H4IRPN=int(hdr.get("W_H4IRPN", 1)),
                     W_H4IRPO=int(hdr.get("W_H4IRPO", 1)),
                     W_H4NRED=int(hdr.get("W_H4NRED", vReads + 1)),
                     W_H4FRMT=float(hdr.get("W_H4FRMT", 0.0)),
                     )
        # A rawISRCube holds the flux accumulated between consecutive reads, so
        # it has one plane fewer than the ramp has reads: the first read is the
        # reference every plane is measured against, not a plane of its own (see
        # PfsIsrTask.makeUTRdeltas). A ramp that breaks this would give the dark a
        # different plane count from the exposures it is subtracted from, which
        # the read-by-read subtraction cannot detect.
        if vReads != vInfo["W_H4NRED"] - 1:
            raise ValueError(
                f"input ramp {v} has {vReads} planes but W_H4NRED={vInfo['W_H4NRED']} "
                f"reads; expected {vInfo['W_H4NRED'] - 1} planes."
            )
        if nreads is None:
            nreads = vReads
            rampInfo = vInfo
        else:
            if vReads != nreads:
                raise ValueError(
                    f"input ramp {v} has {vReads} reads; expected {nreads}. "
                    "Dark ramps must have the same number of reads to combine."
                )
            if vInfo["W_H4IRPN"] != rampInfo["W_H4IRPN"]:
                raise ValueError(
                    f"input ramp {v} has IRP ratio {vInfo['W_H4IRPN']}; "
                    f"expected {rampInfo['W_H4IRPN']}. Dark ramps must share an "
                    "IRP ratio to combine (the ratio sets the read cadence)."
                )
    return rampInfo, nreads


def getRampInfo(butler: dafButler.Butler,
                dataId: dict[str, str],
                vlist: list[int]) -> dict:
    """Read the common IRP/timing cards of the input ramps, without their pixels.

    Only the primary headers are read, so this is cheap enough to call before
    committing to the combine.
    """
    cubes = dict()
    try:
        for v in vlist:
            cubes[v] = fitsio.FITS(butler.getURI('rawISRCube', dataId, visit=v))
        rampInfo, _ = rampInfoFromCubes(cubes, vlist)
    finally:
        for fits in cubes.values():
            fits.close()
    return rampInfo


def collectionRoot(root: str) -> str:
    """Normalize a collection base.

    A trailing slash is stripped: butler collection names are opaque strings, so
    ``u/me/calib/`` would otherwise yield ``u/me/calib//PIPE2D-1664/...`` and be
    registered, silently, under that name.
    """
    root = root.rstrip("/")
    if not root:
        raise ValueError("collection base must not be empty")
    return root


def calibCollectionName(root: str, ticket: str, tag: str, product: str, iteration: str) -> str:
    """The CALIBRATION collection the pipeline selects the calibration from.

    ``{root}/{ticket}/{tag}/{product}.{iteration}``, following DMTN-222
    (https://dmtn-222.lsst.io) as the pipeline-generated ``dark`` and ``bias`` do
    (e.g. ``PFS/calib/pipe2d-1842/run28/dark.20260503a``). Note the absence of the
    ``Gen`` suffix, which belongs to the generation chain.
    """
    return f"{collectionRoot(root)}/{ticket}/{tag}/{product}.{iteration}"


def genCollectionName(root: str, ticket: str, tag: str, product: str, iteration: str) -> str:
    """The CHAINED collection gathering everything one generation produced and used.

    ``{root}/{ticket}/{tag}/{product}Gen.{iteration}``.
    """
    return f"{collectionRoot(root)}/{ticket}/{tag}/{product}Gen.{iteration}"


def runName(genCollection: str, timestamp: str) -> str:
    """The RUN the datasets are written to, within the generation chain.

    Each attempt gets its own timestamped RUN, so a retry after a failed write
    neither collides with nor overwrites the previous one.
    """
    return f"{genCollection}/{timestamp}"


def defaultIteration() -> str:
    """The default DMTN-222 rerun iteration: today's date plus ``a``."""
    return datetime.date.today().strftime("%Y%m%d") + "a"


def defaultTimestamp() -> str:
    """The default RUN timestamp: the current UTC time, as ``YYYYMMDDTHHMMSSZ``."""
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def asTime(when) -> astropy.time.Time | None:
    """Coerce a `datetime.datetime` (or `None`) to an `~astropy.time.Time`.

    Naive datetimes are interpreted as UTC, matching the curated-calibration
    behaviour in `makePfsCalibs`.
    """
    if when is None or isinstance(when, astropy.time.Time):
        return when
    return astropy.time.Time(when, format="datetime", scale="utc")


def ensureCalibCollection(butler: dafButler.Butler, calibCollection: str) -> None:
    """Make sure the CALIBRATION collection to certify into exists."""
    if butler.registry.registerCollection(calibCollection, CollectionType.CALIBRATION):
        logger.info("registered CALIBRATION collection %s", calibCollection)


def ensureGenCollection(butler: dafButler.Butler, genCollection: str,
                        members: list[str]) -> None:
    """Make sure the generation CHAINED collection exists and holds ``members``.

    Members already in the chain keep their position; new ones are prepended, so
    the most recent RUN resolves first. This mirrors what ``pipetask -o`` builds
    for the pipeline-generated calibrations.

    Only this one chain is ever written. It is named for a single product,
    ticket, tag and iteration, so it belongs to this generation alone; the shared
    chains such as ``PFS/calib`` are maintained by hand and are never touched
    here. Promoting a dark into one of those is a separate, deliberate step (see
    `lsst.obs.pfs.datasetCopy`).
    """
    registry = butler.registry
    if registry.registerCollection(genCollection, CollectionType.CHAINED):
        logger.info("registered CHAINED collection %s", genCollection)
    existing = list(registry.getCollectionChain(genCollection))
    added = [member for member in members if member not in existing]
    if added:
        registry.setCollectionChain(genCollection, added + existing)
        logger.info("chained %s into %s", ", ".join(added), genCollection)


def ensureRunCollection(butler: dafButler.Butler, outputRun: str) -> None:
    """Make sure the output collection exists and is a RUN, creating it if absent.

    `Butler.put` cannot write to a CHAINED collection, so a CHAINED ``outputRun``
    is rejected here rather than after the combine.

    Raises
    ------
    ValueError
        If ``outputRun`` exists but is not a RUN collection.
    """
    registry = butler.registry
    try:
        collectionType = registry.getCollectionType(outputRun)
    except MissingCollectionError:
        registry.registerCollection(outputRun, CollectionType.RUN)
        logger.info("registered RUN collection %s", outputRun)
    else:
        if collectionType is not CollectionType.RUN:
            raise ValueError(
                f"output collection {outputRun!r} is of type {collectionType.name}; "
                "a RUN collection is required to write the dark. Pass a RUN collection "
                f"(e.g. {outputRun}/nirDark) as --output; it may then be chained into "
                f"{outputRun!r} afterwards."
            )


def ensureDatasetType(butler: dafButler.Butler, datasetType: str) -> None:
    """Register the given nirDark dataset type if the repository lacks it.

    A newly-introduced ``nirDark_irp<N>`` will not exist in a repository created
    before it was defined, so it is registered on use, as `makePfsCalibs` does
    for its collections. Registration is idempotent.
    """
    registry = butler.registry
    if registry.registerDatasetType(DatasetType(datasetType, DARK_DIMENSIONS,
                                                DARK_STORAGE_CLASS,
                                                universe=registry.dimensions,
                                                isCalibration=True)):
        logger.info("registered dataset type %s", datasetType)


def ensureOutputs(butler: dafButler.Butler, outputRun: str, calibCollection: str,
                  datasetType: str) -> None:
    """Make sure the output collections and dataset type exist and are usable.

    The combine takes hours, so everything that can make the final `Butler.put`
    or `Registry.certify` fail is checked -- and fixed where we legitimately can
    -- up front. Registration is idempotent, so it is harmless for the workers to
    repeat what the wrapper already did.
    """
    ensureCalibCollection(butler, calibCollection)
    ensureRunCollection(butler, outputRun)
    ensureDatasetType(butler, datasetType)


class Plan(NamedTuple):
    """Everything the per-spectrograph combines must agree on, resolved once."""

    datasetType: str
    """``nirDark`` or ``nirDark_irp<N>``, from the inputs' IRP ratio."""
    calibCollection: str
    """The CALIBRATION collection to certify into."""
    genCollection: str
    """The CHAINED collection gathering this generation's run and its input."""
    outputRun: str
    """The timestamped RUN to put the dark into."""
    startDate: datetime.datetime
    """Start of the dark's validity period."""


def planOutputs(butler: dafButler.Butler,
                dataId: dict,
                visits: list[int],
                inputRun: str,
                outputBase: str,
                ticket: str,
                tag: str,
                iteration: str,
                timestamp: str) -> Plan:
    """Name and register the collections this combine will write, before it runs.

    The dataset type follows from the inputs' IRP ratio, and names the product, so
    the input headers -- but not their pixels -- have to be read first.
    """
    datasetType = datasetTypeForIrp(getRampInfo(butler, dataId, visits)["W_H4IRPN"])
    calibCollection = calibCollectionName(outputBase, ticket, tag, datasetType, iteration)
    genCollection = genCollectionName(outputBase, ticket, tag, datasetType, iteration)
    outputRun = runName(genCollection, timestamp)
    ensureOutputs(butler, outputRun, calibCollection, datasetType)
    ensureGenCollection(butler, genCollection, [outputRun, inputRun])
    return Plan(datasetType, calibCollection, genCollection, outputRun,
                getStartDate(butler, dataId, visits[0]))


def preflight(inputRun: str, outputBase: str,
              dataId: dict,
              visits: list[int],
              ticket: str,
              tag: str,
              iteration: str,
              timestamp: str,
              repo_path: str = '/work/datastore') -> Plan:
    """Resolve, once, everything the per-spectrograph combines must agree on.

    Besides naming and registering the output collections, this reads the validity
    start date from the first visit. Each detector latches its own timestamp a
    second or two apart, so letting every spectrograph derive its own start date
    would certify one logical dark set with several staggered validity periods --
    and letting each derive its own RUN timestamp would scatter them across four
    RUNs.
    """
    butler = makeButler(inputRun, None, repo_path)
    return planOutputs(butler, dataId, visits, inputRun, outputBase, ticket, tag,
                       iteration, timestamp)


def saveNirDark(butler: dafButler.Butler,
                dataId: dict[str, str],
                run: str,
                visits: list[int],
                nirCube: np.ndarray,
                start: datetime.datetime,
                readNoise: float,
                gain: float,
                rampInfo: dict,
                saveDir: str | None = None,
                calibCollection: str | None = None):
    """Write a combined NIR dark cube to the butler, and certify it as a calib.

    The output dataset type (``nirDark`` or ``nirDark_irp<N>``) is selected from
    the IRP ratio of the input ramps, and the IRP ratio and read-timing cards are
    recorded in the header so the dark can be matched to exposures with the same
    cadence.

    Parameters
    ----------
    start : `datetime.datetime`
        Start of the dark's validity period; normally the date of the first
        input visit.
    rampInfo : `dict`
        IRP/timing cards common to the input ramps, as returned by
        `makeMasterDark`: ``W_H4IRPN``, ``W_H4IRPO``, ``W_H4NRED``,
        ``W_H4FRMT`` (per-read frame time).
    saveDir : `str`, optional
        If given, the cube is written here as a plain FITS file before it is put
        to the butler, so that a butler failure does not throw away a combine
        that takes hours. The file can be ingested afterwards.
    calibCollection : `str`, optional
        CALIBRATION collection to certify the dark into, so that ISR selects it
        by observation date. If `None` the dark is only put into ``run``. The
        validity period is open-ended: a later dark supersedes this one through
        the calibration chain, not through an end date set here.
    """

    cam = f'n{dataId["spectrograph"]}'
    metadata = dict(GAIN=gain,
                    READNOISE=readNoise,
                    CAM=cam,
                    CALIBDATE=start.strftime('%Y-%m-%dT%H:%M:%S'),
                    )
    # Record the IRP ratio and read cadence so the dark can be matched to
    # exposures read out the same way (the dark is subtracted read-by-read).
    metadata.update(rampInfo)
    setCalibHeader(metadata, 'dark', visitList=visits, outputId=dataId)
    ic = imageCube.ImageCube.fromCube(nirCube, metadata)

    datasetType = datasetTypeForIrp(rampInfo["W_H4IRPN"])

    if saveDir is not None:
        path = os.path.join(saveDir, f"{datasetType}-{cam}.fits")
        t0 = time.time()
        ic.writeFits(path)
        logger.info("%s: wrote %s (%.1fs)", cam, path, time.time() - t0)

    try:
        with butler.transaction():
            ref = butler.put(ic, datasetType, dataId, run=run)
            if calibCollection is not None:
                timespan = Timespan(asTime(start), None)
                butler.registry.certify(calibCollection, [ref], timespan)
                logger.info("%s: certified %s into %s for %s",
                            cam, datasetType, calibCollection, timespan)
    except Exception:
        if saveDir is not None:
            logger.error("%s: butler.put of %s failed; the combined dark is safe in %s",
                         cam, datasetType, path)
        raise


def makeMasterDark(butler: dafButler.Butler,
                   dataId: dict[str, str],
                   vlist: list[int]):
    """Construct a single dark cube from a set of already processed rawISRCubes.

    Parameters
    ----------
    butler : `dafButler.Butler`
       Configured to .get our rawISRCubes.
    dataId : `dict`
       Everything but the visit
    vlist : list[int]
       The visits to merge

    Returns
    -------
    nirDark : `np.ndarray`
       The merged master dark, shaped ``(nreads, height, width)``.
    offsets : `list`
       The offsets, per-read, between the individual input darks.
    rampInfo : `dict`
       IRP/timing cards common to all inputs (``W_H4IRPN``, ``W_H4IRPO``,
       ``W_H4NRED``, ``W_H4FRMT``), for the output header.

    Raises
    ------
    ValueError
        If the input ramps do not all share the same number of reads or the same
        IRP ratio: they cannot be meaningfully combined.
    """

    # Open all input cubes. We then step through them one read at a time,
    # combining them into a single new super dark.
    #
    # We use fitsio directly rather than astropy.io.fits or the ImageCube
    # wrapper: astropy.io.fits internally caches the images, which is
    # catastrophic here. You can flush the cache or bypass it, but it is cleaner
    # just to access the data directly.
    cubes = dict()
    try:
        for v in vlist:
            cubeLoc = butler.getURI('rawISRCube', dataId, visit=v)
            cubes[v] = fitsio.FITS(cubeLoc)

        # Validate that the inputs are compatible, and gather the IRP/timing
        # cards to stamp on the output.
        rampInfo, nreads = rampInfoFromCubes(cubes, vlist)

        ndarks = len(cubes)
        height, width = cubes[vlist[0]]['IMAGE_1'].read().shape

        oneRead = np.zeros(shape=(ndarks, height, width), dtype='f4')
        newCube = np.zeros(shape=(nreads, height, width), dtype='f4')
        offsets = []
        cam = f'n{dataId["spectrograph"]}'
        visits = list(cubes.keys())
        logger.info("%s: combining %d dark ramps of %d reads (%dx%d)",
                    cam, ndarks, nreads, height, width)
        ioTime = procTime = 0.0
        for r_i in range(nreads):
            t0 = time.time()
            for v_i, v in enumerate(visits):
                oneRead[v_i, :, :] = cubes[v][f'IMAGE_{r_i+1}'].read()
            t1 = time.time()

            # Correct any offsets between ramps for this single read. Make sure
            # not to take out the common increase in flux.
            readMeds = np.median(oneRead, axis=(1, 2))
            readMeds -= np.median(readMeds)
            offsets.append(readMeds)
            oneRead -= readMeds[:, None, None]

            # Finish this read
            newCube[r_i] = np.median(oneRead, axis=0)
            t2 = time.time()

            ioTime += t1 - t0
            procTime += t2 - t1
            # Report progress and the single worst ramp offset, so an outlier
            # dark is visible without dumping the whole per-ramp offset list.
            worst = int(np.argmax(np.abs(readMeds)))
            logger.info("%s read %d/%d: max|offset|=%.1f (visit %d)  io=%.1fs proc=%.1fs",
                        cam, r_i + 1, nreads, abs(readMeds[worst]), visits[worst],
                        t1 - t0, t2 - t1)

        logger.info("%s: combined %d reads in %.1fs (read=%.1fs combine=%.1fs)",
                    cam, nreads, ioTime + procTime, ioTime, procTime)
    finally:
        for fits in cubes.values():
            fits.close()

    return newCube, offsets, rampInfo


def makeButler(inputRun: str, outputRun: str | None,
               repo_path: str = '/work/datastore') -> dafButler.Butler:
    """Open a writeable butler reading from ``inputRun`` and writing to ``outputRun``.

    ``outputRun`` may be `None` when the output collection is not yet named, as
    when preflighting.
    """
    collections = [inputRun, 'PFS/defaults']
    if outputRun is not None:
        collections.insert(0, outputRun)
    return dafButler.Butler(repo_path,
                            collections=collections,
                            run=outputRun,
                            writeable=True)


def configureLogging() -> None:
    """Set up console logging for the dark combine and its worker processes.

    Called once in the parent and again, as the worker-pool initializer, in each
    worker. The pool uses a ``forkserver`` start method, so workers begin from a
    clean interpreter and do not inherit the parent's logging configuration;
    reapplying it here keeps their per-read progress lines reaching the terminal.
    """
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")


def processMasterDark(inputRun: str, outputRun: str,
                      dataId: dict,
                      visits: list[int],
                      startDate=None,
                      repo_path='/work/datastore',
                      saveDir=None,
                      calibCollection=None):
    """Combine a spectrograph's dark ramps, then put and certify the result.

    ``startDate`` defaults to the observation date of the first input visit, and
    is both the ``CALIBDATE`` header card and the start of the certified,
    open-ended validity period.
    """

    mdButler = makeButler(inputRun, outputRun, repo_path)

    # Register the outputs before the combine, so that hours of work are not
    # thrown away by a butler problem at the put. The wrapper has already done
    # this for every spectrograph, making these no-ops; repeating it here keeps a
    # direct caller honest.
    if calibCollection is not None:
        ensureCalibCollection(mdButler, calibCollection)
    ensureRunCollection(mdButler, outputRun)
    ensureDatasetType(mdButler, datasetTypeForIrp(getRampInfo(mdButler, dataId,
                                                              visits)["W_H4IRPN"]))

    exp = mdButler.get('postISRCCD', dataId, visit=visits[0])
    expStartTime, gain, readNoise = getExpInfo(exp)
    if startDate is None:
        startDate = expStartTime
    t0 = time.time()
    darkCube, offsets, rampInfo = makeMasterDark(mdButler, dataId, visits)
    t1 = time.time()
    saveNirDark(mdButler, dataId, outputRun, visits,
                darkCube, start=startDate,
                readNoise=readNoise, gain=gain, rampInfo=rampInfo,
                saveDir=saveDir, calibCollection=calibCollection)
    t2 = time.time()
    del darkCube
    logger.info("n%s: make=%0.1fs save=%0.1fs", dataId["spectrograph"], t1 - t0, t2 - t1)
