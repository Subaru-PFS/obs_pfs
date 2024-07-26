from typing import Any, Dict, Set, Mapping, Tuple, Union

from astro_metadata_translator import ObservationInfo
from lsst.daf.butler import DataCoordinate, DatasetType, DimensionRecord, DimensionUniverse
from lsst.obs.base.ingest import RawIngestTask, RawFileDatasetInfo
from lsst.resources import ResourcePath


class PfsRawIngestTask(RawIngestTask):
    """Driver Task for ingesting raw data into Gen3 Butler repositories.

    This subclass is specialised for ingesting raw PFS data. This is required
    in order to read PFS-specific dimension data.

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
            ("instrument", "spectrograph", "arm", "exposure", "detector"),
            "PfsRaw",
            universe=self.butler.dimensions,
        )

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
        return required, optional

    def makeExposureRecord(
        self, obsInfo: ObservationInfo, universe: DimensionUniverse, **kwargs
    ) -> DimensionRecord:
        """Construct a registry record for an exposure

        This adds PFS-specific fields to the exposure record.

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
        return super().makeExposureRecord(
            obsInfo,
            universe,
            arm=obsInfo.ext_arm,  # type: ignore[attr-defined]
            spectrograph=obsInfo.ext_spectrograph,  # type: ignore[attr-defined]
            pfs_design_id=obsInfo.ext_pfs_design_id,  # type: ignore[attr-defined]
            dither=obsInfo.ext_dither,  # type: ignore[attr-defined]
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
            value=obsInfo.ext_dither,  # type: ignore[attr-defined]
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
            exposure=obsInfo.exposure_id,
            detector=obsInfo.detector_num,
            arm=obsInfo.ext_arm,  # type: ignore[attr-defined]
            spectrograph=obsInfo.ext_spectrograph,  # type: ignore[attr-defined]
            universe=self.universe,
        )
        return RawFileDatasetInfo(obsInfo=obsInfo, dataId=dataId)
