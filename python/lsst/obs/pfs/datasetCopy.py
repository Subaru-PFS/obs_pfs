"""Copy or link datasets between collections within one butler repository.

A dataset's RUN collection is immutable: the butler offers no move or rename, and
``Butler.transfer_from`` is a no-op when its source and destination are the same
repository. The only way to re-home a dataset is to give its file a new dataset
ID under a new RUN, which is what `copyDatasets` does.

Where a genuine copy is not wanted, `linkDatasets` leaves the datasets in the RUN
they were written to. Either way, calibrations may then be certified into a
CALIBRATION collection with `certifyDatasets`, which is what makes the pipeline
select them by observation date.
"""

from __future__ import annotations

import logging

import astropy.time

from lsst.daf.butler import Butler, CollectionType, DatasetRef, FileDataset, Timespan
from lsst.daf.butler.registry import ConflictingDefinitionError, MissingCollectionError

logger = logging.getLogger(__name__)

# Bulky intermediates that exist only to be consumed by a later step. Promoting
# them into a shared collection is almost always a mistake, so copying one has to
# be asked for explicitly.
INTERMEDIATE_DATASET_TYPES = frozenset({"rawISRCube"})


def findDatasets(butler: Butler,
                 datasetTypes: list[str],
                 inputCollections: list[str],
                 where: str = "") -> list[DatasetRef]:
    """Find the datasets to copy or link.

    Only the first match of each dataset type and dataId is returned, so a
    CHAINED input collection resolves the way the pipeline would resolve it.
    """
    refs: list[DatasetRef] = []
    for datasetType in datasetTypes:
        found = butler.query_datasets(datasetType, collections=inputCollections,
                                      find_first=True, where=where, explain=False)
        logger.info("found %d %s dataset(s) in %s",
                    len(found), datasetType, ", ".join(inputCollections))
        refs.extend(found)
    if not refs:
        raise RuntimeError(
            f"no datasets of type(s) {', '.join(datasetTypes)} found in "
            f"{', '.join(inputCollections)}"
            + (f" matching {where!r}" if where else "")
        )
    return refs


def ensureRun(butler: Butler, outputRun: str) -> None:
    """Make sure ``outputRun`` exists and is a RUN collection."""
    registry = butler.registry
    try:
        collectionType = registry.getCollectionType(outputRun)
    except MissingCollectionError:
        registry.registerCollection(outputRun, CollectionType.RUN)
        logger.info("registered RUN collection %s", outputRun)
    else:
        if collectionType is not CollectionType.RUN:
            raise ValueError(f"output collection {outputRun!r} is of type "
                             f"{collectionType.name}; a RUN collection is required")


def artifactUri(butler: Butler, ref: DatasetRef):
    """The single file backing a dataset.

    Raises
    ------
    NotImplementedError
        If the dataset is stored as a disassembled composite, which has no single
        artifact to copy.
    """
    uris = butler.getURIs(ref)
    if uris.componentURIs:
        raise NotImplementedError(
            f"{ref.datasetType.name} {ref.dataId} is stored as a disassembled composite "
            f"({len(uris.componentURIs)} components); copying it is not supported")
    return uris.primaryURI


def copyDatasets(butler: Butler,
                 refs: list[DatasetRef],
                 outputRun: str,
                 skipExisting: bool = False,
                 transferMode: str = "copy") -> dict[DatasetRef, DatasetRef]:
    """Copy datasets into a new RUN collection, returning a source -> copy mapping.

    The file backing each dataset is copied and ingested under a fresh dataset ID,
    so the result does not depend on the source collection surviving.

    The obvious ``butler.put(butler.get(ref))`` is wrong here. It round-trips the
    data through the storage class, which for `ImageCube` silently drops every
    image -- `ImageCube.writeFits` only writes images that were explicitly read or
    set, and a freshly-fetched cube has none. It would also deserialize a
    multi-gigabyte NIR dark ramp merely to serialize it again. Copying the
    artifact is byte-exact and cheap.

    Parameters
    ----------
    skipExisting : `bool`
        Skip datasets already present in ``outputRun`` rather than failing, so
        that an interrupted copy can be resumed.
    transferMode : `str`
        How to transfer the file: ``copy`` (the default, and the only mode that
        leaves the destination independent of the source), or any other mode
        `Butler.ingest` accepts, such as ``hardlink`` or ``direct``.
    """
    ensureRun(butler, outputRun)
    copied: dict[DatasetRef, DatasetRef] = {}
    for i, ref in enumerate(refs, 1):
        existing = butler.query_datasets(ref.datasetType, collections=[outputRun],
                                         data_id=ref.dataId, explain=False, limit=1)
        if existing:
            if not skipExisting:
                raise ConflictingDefinitionError(
                    f"{ref.datasetType.name} {ref.dataId} already exists in {outputRun}; "
                    "pass skipExisting to resume an interrupted copy")
            logger.info("[%d/%d] %s %s already in %s; skipping",
                        i, len(refs), ref.datasetType.name, ref.dataId, outputRun)
            copied[ref] = existing[0]
            continue
        newRef = DatasetRef(ref.datasetType, ref.dataId, run=outputRun)
        butler.ingest(FileDataset(path=artifactUri(butler, ref), refs=[newRef]),
                      transfer=transferMode)
        logger.info("[%d/%d] copied %s %s -> %s",
                    i, len(refs), ref.datasetType.name, ref.dataId, outputRun)
        copied[ref] = newRef
    return copied


def linkDatasets(butler: Butler, refs: list[DatasetRef]) -> dict[DatasetRef, DatasetRef]:
    """Use the datasets where they already are, without copying them.

    The datasets keep their existing RUN, so anything certifying them stays
    dependent on that RUN continuing to exist.
    """
    for ref in refs:
        logger.info("linking %s %s in %s", ref.datasetType.name, ref.dataId, ref.run)
    return {ref: ref for ref in refs}


def sourceTimespans(butler: Butler,
                    refs: list[DatasetRef],
                    inputCollections: list[str]) -> dict[DatasetRef, Timespan]:
    """Read back the validity timespans the datasets were certified with.

    Used when promoting an already-certified calibration, so that its validity
    period need not be restated. Only CALIBRATION collections carry timespans.

    Raises
    ------
    RuntimeError
        If any of ``refs`` is not certified in ``inputCollections``.
    """
    # queryDatasetAssociations is the only route to a dataset's timespan; its own
    # docstring calls it a placeholder, so it is confined to this function.
    wanted = {ref.id: ref for ref in refs}
    found: dict[DatasetRef, Timespan] = {}
    for datasetType in {ref.datasetType.name for ref in refs}:
        for assoc in butler.registry.queryDatasetAssociations(
                datasetType, collections=inputCollections,
                collectionTypes=[CollectionType.CALIBRATION], flattenChains=True):
            if assoc.ref.id in wanted and assoc.timespan is not None:
                found[wanted[assoc.ref.id]] = assoc.timespan
    missing = [ref for ref in refs if ref not in found]
    if missing:
        raise RuntimeError(
            f"{len(missing)} dataset(s) are not certified in "
            f"{', '.join(inputCollections)}, so their validity period cannot be "
            "inherited; pass an explicit begin date. First: "
            f"{missing[0].datasetType.name} {missing[0].dataId}")
    return found


def certifyDatasets(butler: Butler,
                    calibCollection: str,
                    timespans: dict[DatasetRef, Timespan]) -> None:
    """Certify datasets into a CALIBRATION collection with their timespans.

    Datasets sharing a timespan are certified together.
    """
    if butler.registry.registerCollection(calibCollection, CollectionType.CALIBRATION):
        logger.info("registered CALIBRATION collection %s", calibCollection)

    grouped: dict[Timespan, list[DatasetRef]] = {}
    for ref, timespan in timespans.items():
        grouped.setdefault(timespan, []).append(ref)
    with butler.transaction():
        for timespan, group in grouped.items():
            butler.registry.certify(calibCollection, group, timespan)
            logger.info("certified %d dataset(s) into %s for %s",
                        len(group), calibCollection, timespan)


def isoTai(value: str | None) -> astropy.time.Time | None:
    """Parse an ISO-8601 TAI string into a `~astropy.time.Time` (or `None`)."""
    if value is None:
        return None
    return astropy.time.Time(value, format="isot", scale="tai")


def transfer(butler: Butler,
             datasetTypes: list[str],
             inputCollections: list[str],
             outputRun: str | None = None,
             where: str = "",
             calibCollection: str | None = None,
             beginDate: str | None = None,
             endDate: str | None = None,
             skipExisting: bool = False,
             transferMode: str = "copy",
             allowIntermediates: bool = False,
             dryRun: bool = False) -> list[DatasetRef]:
    """Copy (or link) datasets into ``outputRun``, optionally certifying them.

    Parameters
    ----------
    outputRun : `str`, optional
        RUN collection to copy the datasets into. If `None` the datasets are
        left in the RUN they are already in.
    calibCollection : `str`, optional
        CALIBRATION collection to certify the results into.
    beginDate, endDate : `str`, optional
        ISO-8601 TAI bounds of the validity period. If ``beginDate`` is omitted,
        the timespans are inherited from the source certification.
    allowIntermediates : `bool`
        Permit copying an `INTERMEDIATE_DATASET_TYPES` member, such as the bulky
        ``rawISRCube`` ramps that a nirDark is combined from. Refused by default.

    Raises
    ------
    ValueError
        If an intermediate dataset type is requested without
        ``allowIntermediates``.
    """
    if not allowIntermediates:
        intermediates = sorted(set(datasetTypes) & INTERMEDIATE_DATASET_TYPES)
        if intermediates:
            raise ValueError(
                f"{', '.join(intermediates)} is an intermediate product, not something "
                "to promote into a shared collection; pass allowIntermediates to copy it "
                "anyway")

    refs = findDatasets(butler, datasetTypes, inputCollections, where=where)

    # Resolve the source timespans before copying: the copies have new dataset
    # IDs, and only the originals are certified anywhere.
    sourceSpans: dict[DatasetRef, Timespan] | None = None
    if calibCollection is not None:
        if beginDate is None:
            sourceSpans = sourceTimespans(butler, refs, inputCollections)
        else:
            sourceSpans = dict.fromkeys(refs, Timespan(isoTai(beginDate), isoTai(endDate)))

    if dryRun:
        for ref in refs:
            where_ = ref.run if outputRun is None else outputRun
            logger.info("would %s %s %s -> %s", "link" if outputRun is None else "copy",
                        ref.datasetType.name, ref.dataId, where_)
            if calibCollection is not None:
                logger.info("    then certify into %s for %s", calibCollection, sourceSpans[ref])
        return []

    if outputRun is None:
        copied = linkDatasets(butler, refs)
    else:
        copied = copyDatasets(butler, refs, outputRun, skipExisting=skipExisting,
                              transferMode=transferMode)

    if calibCollection is not None:
        certifyDatasets(butler, calibCollection,
                        {copied[ref]: span for ref, span in sourceSpans.items()})

    return list(copied.values())
