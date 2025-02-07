from collections import defaultdict
from multiprocessing import Pool
from typing import Any, Dict, Set, Mapping, Tuple, Union

from astro_metadata_translator import ObservationInfo
from lsst.daf.butler import CollectionType, DataCoordinate, DatasetType, DimensionRecord, DimensionUniverse
from lsst.obs.base.ingest import RawIngestTask, RawFileDatasetInfo, RawExposureData
from lsst.obs.base.ingest import _log_msg_counter
from lsst.resources import ResourcePath


def quantizeDither(dither: float) -> int:
    """Quantize a dither value.

    Parameters
    ----------
    dither : `float`
        Dither value.

    Returns
    -------
    quantized : `int`
        Quantized dither value.
    """
    return int(dither * 1e4 + 0.5)


class PfsRawIngestTask(RawIngestTask):
    """Driver Task for ingesting raw data into Gen3 Butler repositories.

    This subclass is specialised for ingesting raw PFS data. This is required
    in order to read PFS-specific dimension data into the "visit" table.

    Parameters
    ----------
    config : `RawIngestConfig`
        Configuration for the task.
    butler : `~lsst.daf.butler.Butler`
        Writeable butler instance, with ``butler.run`` set to the appropriate
        `~lsst.daf.butler.CollectionType.RUN` collection for these raw
        datasets.
    on_success : `Callable`, optional
        A callback invoked when all of the raws associated with an exposure
        are ingested.  Will be passed a list of `FileDataset` objects, each
        containing one or more resolved `DatasetRef` objects.  If this callback
        raises it will interrupt the entire ingest process, even if
        `RawIngestConfig.failFast` is `False`.
    on_metadata_failure : `Callable`, optional
        A callback invoked when a failure occurs trying to translate the
        metadata for a file.  Will be passed the URI and the exception, in
        that order, as positional arguments.  Guaranteed to be called in an
        ``except`` block, allowing the callback to re-raise or replace (with
        ``raise ... from``) to override the task's usual error handling (before
        `RawIngestConfig.failFast` logic occurs).
    on_ingest_failure : `Callable`, optional
        A callback invoked when dimension record or dataset insertion into the
        database fails for an exposure.  Will be passed a `RawExposureData`
        instance and the exception, in that order, as positional arguments.
        Guaranteed to be called in an ``except`` block, allowing the callback
        to re-raise or replace (with ``raise ... from``) to override the task's
        usual error handling (before `RawIngestConfig.failFast` logic occurs).
    **kwargs
        Additional keyword arguments are forwarded to the `lsst.pipe.base.Task`
        constructor.

    Notes
    -----
    Each instance of `RawIngestTask` writes to the same Butler.  Each
    invocation of `RawIngestTask.run` ingests a list of files.
    """

    def getDatasetType(self) -> DatasetType:
        """Return the default DatasetType of the datasets ingested by this
        Task.

        This override sets the dimensions for the "raw" dataset type, replacing
        the default "detector" dimension with "spectrograph" and "arm", and uses
        our custom `PfsRaw` dataset type.

        Returns
        -------
        datasetType : `DatasetType`
            The default dataset type to use for the data being ingested. This
            is only used if the relevant `~lsst.pipe.base.Instrument` does not
            define an override.
        """
        return DatasetType(
            "raw",
            ("instrument", "visit", "spectrograph", "arm"),
            "PfsRaw",
            universe=self.butler.dimensions,
        )

    def groupByExposure(self, files):
        """Group an iterable of `RawFileData` by exposure.

        Parameters
        ----------
        files : iterable of `RawFileData`
            File-level information to group.

        Returns
        -------
        exposures : `list` of `RawExposureData`
            A list of structures that group the file-level information by
            exposure. All fields will be populated.  The
            `RawExposureData.dataId` attributes will be minimal (unexpanded)
            `~lsst.daf.butler.DataCoordinate` instances.
        """
        exposureDimensions = self.universe["visit"].dimensions
        byExposure = defaultdict(list)
        for f in files:
            # Assume that the first dataset is representative for the file.
            byExposure[f.datasets[0].dataId.subset(exposureDimensions)].append(f)

        return [
            RawExposureData(
                dataId=dataId,
                files=exposureFiles,
                universe=self.universe,
                record=self.makeExposureRecord(exposureFiles[0].datasets[0].obsInfo, self.universe),
                dependencyRecords=self.makeDependencyRecords(
                    exposureFiles[0].datasets[0].obsInfo, self.universe
                ),
            )
            for dataId, exposureFiles in byExposure.items()
        ]

    def expandDataIds(self, data: RawExposureData) -> RawExposureData:
        """Expand the data IDs associated with a raw exposure.

        This adds the metadata records.

        Parameters
        ----------
        exposure : `RawExposureData`
            A structure containing information about the exposure to be
            ingested.  Must have `RawExposureData.record` populated. Should
            be considered consumed upon return.

        Returns
        -------
        exposure : `RawExposureData`
            An updated version of the input structure, with
            `RawExposureData.dataId` and nested `RawFileData.dataId` attributes
            updated to data IDs for which
            `~lsst.daf.butler.DataCoordinate.hasRecords` returns `True`.
        """
        # We start by expanded the exposure-level data ID; we won't use that
        # directly in file ingest, but this lets us do some database lookups
        # once per exposure instead of once per file later.
        data.dataId = self.butler.registry.expandDataId(
            data.dataId,
            # We pass in the records we'll be inserting shortly so they aren't
            # looked up from the database.  We do expect instrument and filter
            # records to be retrieved from the database here (though the
            # Registry may cache them so there isn't a lookup every time).
            records={"visit": data.record},
        )
        # Now we expand the per-file (exposure+detector) data IDs.  This time
        # we pass in the records we just retrieved from the exposure data ID
        # expansion.
        for file in data.files:
            for dataset in file.datasets:
                dataset.dataId = self.butler.registry.expandDataId(
                    dataset.dataId, records=data.dataId.records
                )
        return data

    def ingestFiles(  # type:ignore
        self,
        files,
        *,
        pool,
        processes: int = 1,
        run: str | None = None,
        skip_existing_exposures: bool = False,
        update_exposure_records: bool = False,
        track_file_attrs: bool = True,
    ):
        """Ingest files into a Butler data repository.

        This creates any new exposure or visit Dimension entries needed to
        identify the ingested files, creates new Dataset entries in the
        Registry and finally ingests the files themselves into the Datastore.
        Any needed instrument, detector, and physical_filter Dimension entries
        must exist in the Registry before `run` is called.

        Parameters
        ----------
        files : iterable over `lsst.resources.ResourcePath`
            URIs to the files to be ingested.
        pool : `multiprocessing.Pool`, optional
            If not `None`, a process pool with which to parallelize some
            operations.
        processes : `int`, optional
            The number of processes to use.  Ignored if ``pool`` is not `None`.
        run : `str`, optional
            Name of a RUN-type collection to write to, overriding
            the default derived from the instrument name.
        skip_existing_exposures : `bool`, optional
            If `True` (`False` is default), skip raws that have already been
            ingested (i.e. raws for which we already have a dataset with the
            same data ID in the target collection, even if from another file).
            Note that this is much slower than just not passing
            already-ingested files as inputs, because we still need to read and
            process metadata to identify which exposures to search for.  It
            also will not work reliably if multiple processes are attempting to
            ingest raws from the same exposure concurrently, in that different
            processes may still attempt to ingest the same raw and conflict,
            causing a failure that prevents other raws from the same exposure
            from being ingested.
        update_exposure_records : `bool`, optional
            If `True` (`False` is default), update existing exposure records
            that conflict with the new ones instead of rejecting them.  THIS IS
            AN ADVANCED OPTION THAT SHOULD ONLY BE USED TO FIX METADATA THAT IS
            KNOWN TO BE BAD.  This should usually be combined with
            ``skip_existing_exposures=True``.
        track_file_attrs : `bool`, optional
            Control whether file attributes such as the size or checksum should
            be tracked by the datastore. Whether this parameter is honored
            depends on the specific datastore implentation.

        Returns
        -------
        refs : `list` of `lsst.daf.butler.DatasetRef`
            Dataset references for ingested raws.
        bad_files : `list` of `ResourcePath`
            Given paths that could not be ingested.
        n_exposures : `int`
            Number of exposures successfully ingested.
        n_exposures_failed : `int`
            Number of exposures that failed when inserting dimension data.
        n_ingests_failed : `int`
            Number of exposures that failed when ingesting raw datasets.
        """
        created_pool = False
        if pool is None and processes > 1:
            pool = Pool(processes)
            created_pool = True

        try:
            exposureData, bad_files = self.prep(files, pool=pool)
        finally:
            if created_pool and pool:
                # The pool is not needed any more so close it if we created
                # it to ensure we clean up resources.
                pool.close()
                pool.join()

        # Up to this point, we haven't modified the data repository at all.
        # Now we finally do that, with one transaction per exposure.  This is
        # not parallelized at present because the performance of this step is
        # limited by the database server.  That may or may not change in the
        # future once we increase our usage of bulk inserts and reduce our
        # usage of savepoints; we've tried to get everything but the database
        # operations done in advance to reduce the time spent inside
        # transactions.
        refs = []
        runs = set()
        datasetTypes: dict[str, DatasetType] = {}
        n_exposures = 0
        n_exposures_failed = 0
        n_ingests_failed = 0
        for exposure in self.progress.wrap(exposureData, desc="Ingesting raw exposures"):
            assert exposure.record is not None, "Should be guaranteed by prep()"
            self.log.debug(
                "Attempting to ingest %d file%s from exposure %s:%s",
                *_log_msg_counter(exposure.files),
                exposure.record.instrument,
                exposure.record.obs_id,
            )

            try:
                for name, record in exposure.dependencyRecords.items():
                    self.butler.registry.syncDimensionData(name, record, update=update_exposure_records)
                inserted_or_updated = self.butler.registry.syncDimensionData(
                    "visit",
                    exposure.record,
                    update=update_exposure_records,
                )
            except Exception as e:
                self._on_ingest_failure(exposure, e)
                n_exposures_failed += 1
                self.log.warning(
                    "Exposure %s:%s could not be registered: %s",
                    exposure.record.instrument,
                    exposure.record.obs_id,
                    e,
                )
                if self.config.failFast:
                    raise e
                continue

            if isinstance(inserted_or_updated, dict):
                # Exposure is in the registry and we updated it, so
                # syncDimensionData returned a dict.
                self.log.info(
                    "Exposure %s:%s was already present, but columns %s were updated.",
                    exposure.record.instrument,
                    exposure.record.obs_id,
                    str(list(inserted_or_updated.keys())),
                )

            # Determine the instrument so we can work out the dataset type.
            instrument = exposure.files[0].instrument
            assert (
                instrument is not None
            ), "file should have been removed from this list by prep if instrument could not be found"

            if raw_definition := getattr(instrument, "raw_definition", None):
                datasetTypeName, dimensions, storageClass = raw_definition
                if not (datasetType := datasetTypes.get(datasetTypeName)):
                    datasetType = DatasetType(
                        datasetTypeName, dimensions, storageClass, universe=self.butler.dimensions
                    )
            else:
                datasetType = self.datasetType
            if datasetType.name not in datasetTypes:
                self.butler.registry.registerDatasetType(datasetType)
                datasetTypes[datasetType.name] = datasetType

            # Override default run if nothing specified explicitly.
            if run is None:
                this_run = instrument.makeDefaultRawIngestRunName()
            else:
                this_run = run
            if this_run not in runs:
                self.butler.registry.registerCollection(this_run, type=CollectionType.RUN)
                runs.add(this_run)
            try:
                datasets_for_exposure = self.ingestExposureDatasets(
                    exposure,
                    datasetType=datasetType,
                    run=this_run,
                    skip_existing_exposures=skip_existing_exposures,
                    track_file_attrs=track_file_attrs,
                )
            except Exception as e:
                self._on_ingest_failure(exposure, e)
                n_ingests_failed += 1
                self.log.warning("Failed to ingest the following for reason: %s", e)
                for f in exposure.files:
                    self.log.warning("- %s", f.filename)
                if self.config.failFast:
                    raise e
                continue
            else:
                self._on_success(datasets_for_exposure)
                for dataset in datasets_for_exposure:
                    refs.extend(dataset.refs)

            # Success for this exposure.
            n_exposures += 1
            self.log.info(
                "Exposure %s:%s ingested successfully", exposure.record.instrument, exposure.record.obs_id
            )

        return refs, bad_files, n_exposures, n_exposures_failed, n_ingests_failed

    @classmethod
    def getObservationInfoSubsets(cls) -> Tuple[Set, Set]:
        """Return subsets of fields in the `ObservationInfo` that we care about

        These fields will be used in constructing an exposure record.

        Returns
        -------
        required : `set`
            Set of `ObservationInfo` field names that are required.
        optional : `set`
            Set of `ObservationInfo` field names we will use if they are
            available.
        """
        required, optional = super().getObservationInfoSubsets()
        required |= {"ext_arm", "ext_spectrograph", "ext_pfs_design_id", "ext_dither"}
        optional |= {"ext_shift", "ext_focus", "ext_lamps"}
        return required, optional

    def makeExposureRecord(
        self, obsInfo: ObservationInfo, universe: DimensionUniverse, **kwargs
    ) -> DimensionRecord:
        """Construct a registry record for an exposure

        This adds PFS-specific fields to the exposure record.
        Copied from lsst.obs.base._instrument.makeExposureRecordFromObsInfo.

        Parameters
        ----------
        obsInfo : `ObservationInfo`
            Observation details for (one of the components of) the exposure.
        universe : `DimensionUniverse`
            Set of all known dimensions.
        **kwargs
            Additional field values for this record.

        Returns
        -------
        record : `DimensionRecord`
            The exposure record that must be inserted into the
            `~lsst.daf.butler.Registry` prior to file-level ingest.
        """
        dimension = universe["visit"]

        # Some registries support additional items.
        supported = {meta.name for meta in dimension.metadata}

        ra, dec, sky_angle, azimuth, zenith_angle = (None, None, None, None, None)
        if obsInfo.tracking_radec is not None:
            icrs = obsInfo.tracking_radec.icrs
            ra = icrs.ra.degree
            dec = icrs.dec.degree
            if obsInfo.boresight_rotation_coord == "sky":
                sky_angle = obsInfo.boresight_rotation_angle.degree
        if obsInfo.altaz_begin is not None:
            zenith_angle = obsInfo.altaz_begin.zen.degree
            azimuth = obsInfo.altaz_begin.az.degree

        extras: dict[str, Any] = {}
        for meta_key, info_key in (
            ("has_simulated", "has_simulated_content"),
            ("seq_start", "group_counter_start"),
            ("seq_end", "group_counter_end"),
        ):
            if meta_key in supported:
                extras[meta_key] = getattr(obsInfo, info_key)

        if (k := "azimuth") in supported:
            extras[k] = azimuth

        return dimension.RecordClass(
            instrument=obsInfo.instrument,
            id=obsInfo.exposure_id,
            obs_id=obsInfo.observation_id,
            group_name=obsInfo.exposure_group,
            group_id=obsInfo.visit_id,
            datetime_begin=obsInfo.datetime_begin,
            datetime_end=obsInfo.datetime_end,
            exposure_time=obsInfo.exposure_time.to_value("s"),
            # we are not mandating that dark_time be calculable
            dark_time=obsInfo.dark_time.to_value("s") if obsInfo.dark_time is not None else None,
            observation_type=obsInfo.observation_type,
            observation_reason=obsInfo.observation_reason,
            day_obs=obsInfo.observing_day,
            seq_num=obsInfo.observation_counter,
            science_program=obsInfo.science_program,
            target_name=obsInfo.object,
            tracking_ra=ra,
            tracking_dec=dec,
            sky_angle=sky_angle,
            zenith_angle=zenith_angle,
            pfs_design_id=obsInfo.ext_pfs_design_id,  # type: ignore[attr-defined]
            dither=quantizeDither(obsInfo.ext_dither),  # type: ignore[attr-defined]
            shift=obsInfo.ext_shift,  # type: ignore[attr-defined]
            focus=obsInfo.ext_focus,  # type: ignore[attr-defined]
            lamps=obsInfo.ext_lamps,  # type: ignore[attr-defined]
            **extras,
            **kwargs,
        )

    def makeDependencyRecords(
        self, obsInfo: ObservationInfo, universe: DimensionUniverse
    ) -> Dict[str, DimensionRecord]:
        """Construct dependency records

        These dependency records will be inserted into the
        `~lsst.daf.butler.Registry` before the exposure records, because they
        are dependencies of the exposure. This allows an opportunity to satisfy
        foreign key constraints that exist because of dimensions related to the
        exposure.

        This implementation adds the required records for the ``dither`` and
        ``pfs_design_id`` dimension.

        Parameters
        ----------
        obsInfo : `ObservationInfo`
            Observation details for (one of the components of) the exposure.
        universe : `DimensionUniverse`
            Set of all known dimensions.

        Returns
        -------
        records : `dict` [`str`, `DimensionRecord`]
            The records to insert, indexed by dimension name.
        """
        arm = universe["arm"].RecordClass(
            name=obsInfo.ext_arm,  # type: ignore[attr-defined]
            instrument=obsInfo.instrument,
        )
        spectrograph = universe["spectrograph"].RecordClass(
            num=obsInfo.ext_spectrograph,  # type: ignore[attr-defined]
            instrument=obsInfo.instrument,
        )
        dither = universe["dither"].RecordClass(
            value=quantizeDither(obsInfo.ext_dither),  # type: ignore[attr-defined]
            instrument=obsInfo.instrument,
        )
        pfs_design_id = universe["pfs_design_id"].RecordClass(
            value=obsInfo.ext_pfs_design_id,  # type: ignore[attr-defined]
            instrument=obsInfo.instrument,
        )
        return dict(arm=arm, spectrograph=spectrograph, dither=dither, pfs_design_id=pfs_design_id)

    def _calculate_dataset_info(
        self, header: Union[Mapping[str, Any], ObservationInfo], filename: ResourcePath
    ) -> RawFileDatasetInfo:
        """Calculate a RawFileDatasetInfo from the supplied information.

        This naughty override of a private member function is required to get
        the PFS-specific observation information into the dataId.

        Parameters
        ----------
        header : Mapping or `astro_metadata_translator.ObservationInfo`
            Header from the dataset or previously-translated content.
        filename : `lsst.resources.ResourcePath`
            Filename to use for error messages.

        Returns
        -------
        dataset : `RawFileDatasetInfo`
            The dataId, and observation information associated with this
            dataset.
        """
        required, optional = self.getObservationInfoSubsets()
        if isinstance(header, ObservationInfo):
            obsInfo = header
            missing = []
            # Need to check the required properties are present.
            for property in required:
                # getattr does not need to be protected because it is using
                # the defined list above containing properties that must exist.
                value = getattr(obsInfo, property)
                if value is None:
                    missing.append(property)
            if missing:
                raise ValueError(
                    f"Requested required properties are missing from file {filename}: {missing} (via JSON)"
                )

        else:
            obsInfo = ObservationInfo(
                header,  # type:ignore
                pedantic=False,
                filename=str(filename),
                required=required,
                subset=required | optional,
            )

        dataId = DataCoordinate.standardize(
            instrument=obsInfo.instrument,
            visit=obsInfo.exposure_id,
            detector=obsInfo.detector_num,
            arm=obsInfo.ext_arm,  # type: ignore[attr-defined]
            spectrograph=obsInfo.ext_spectrograph,  # type: ignore[attr-defined]
            universe=self.universe,
        )
        return RawFileDatasetInfo(obsInfo=obsInfo, dataId=dataId)
