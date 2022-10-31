import lsst.daf.persistence
import lsst.daf.persistence.butler
import lsst.daf.persistence.storage
import lsst.daf.persistence.posixStorage
from lsst.daf.persistence.posixStorage import readFitsStorage, writeFitsStorage
from lsst.daf.persistence.posixStorage import readParquetStorage, writeParquetStorage
from lsst.daf.persistence.posixStorage import readConfigStorage, writeConfigStorage
from lsst.daf.persistence.posixStorage import readPickleStorage, writePickleStorage
from lsst.daf.persistence.posixStorage import readFitsCatalogStorage, writeFitsCatalogStorage
from lsst.daf.persistence.posixStorage import readMatplotlibStorage, writeMatplotlibStorage
from lsst.daf.persistence.posixStorage import readYamlStorage, writeYamlStorage
from lsst.daf.persistence import (
    Butler, ButlerComposite, DataId, LogicalLocation, NoResults, PosixStorage, doImport
)


class PfsButler(Butler):
    def datasetExists(self, datasetType, dataId={}, write=False, **rest):
        """Determines if a dataset file exists.

        Parameters
        ----------
        datasetType - string
            The type of dataset to inquire about.
        dataId - DataId, dict
            The data id of the dataset.
        write - bool
            If True, look only in locations where the dataset could be written,
            and return True only if it is present in all of them.
        **rest keyword arguments for the data id.

        Returns
        -------
        exists - bool
            True if the dataset exists or is non-file-based.
        """
        datasetType = self._resolveDatasetTypeAlias(datasetType)
        dataId = DataId(dataId)
        dataId.update(**rest)
        locations = self._locate(datasetType, dataId, write=write, checkExistence=True)
        if not write:  # when write=False, locations is not a sequence
            if locations is None:
                return False
            locations = [locations]

        if not locations:  # empty list
            return False

        for location in locations:
            # If the location is a ButlerComposite (as opposed to a ButlerLocation),
            # verify the component objects exist.
            if isinstance(location, ButlerComposite):
                for name, componentInfo in location.componentInfo.items():
                    if componentInfo.subset:
                        subset = self.subset(datasetType=componentInfo.datasetType, dataId=location.dataId)
                        exists = all([obj.datasetExists() for obj in subset])
                    else:
                        exists = self.datasetExists(componentInfo.datasetType, location.dataId)
                    if exists is False:
                        return False
            else:
                if not location.repository.exists(location):
                    return False
        return True

    def _locate(self, datasetType, dataId, write, checkExistence=False):
        """Get one or more ButlerLocations and/or ButlercComposites.

        Parameters
        ----------
        datasetType : string
            The datasetType that is being searched for. The datasetType may be followed by a dot and
            a component name (component names are specified in the policy). IE datasetType.componentName

        dataId : dict or DataId class instance
            The dataId

        write : bool
            True if this is a search to write an object. False if it is a search to read an object. This
            affects what type (an object or a container) is returned.

        checkExistence : bool
            Only check if the file exists (useful with bypass functions)

        Returns
        -------
        If write is False, will return either a single object or None. If write is True, will return a list
        (which may be empty)
        """
        repos = self._repos.outputs() if write else self._repos.inputs()
        locations = []
        for repoData in repos:
            # enforce dataId & repository tags when reading:
            if not write and dataId.tag and len(dataId.tag.intersection(repoData.tags)) == 0:
                continue
            components = datasetType.split('.')
            datasetType = components[0]
            components = components[1:]
            try:
                location = repoData.repo.map(datasetType, dataId, write=write)
            except NoResults:
                continue
            if location is None:
                continue
            location.datasetType = datasetType  # todo is there a better way than monkey patching here?
            if len(components) > 0:
                if not isinstance(location, ButlerComposite):
                    raise RuntimeError("The location for a dotted datasetType must be a composite.")
                # replace the first component name with the datasetType
                components[0] = location.componentInfo[components[0]].datasetType
                # join components back into a dot-delimited string
                datasetType = '.'.join(components)
                location = self._locate(datasetType, dataId, write)
                # if a component location is not found, we can not continue with this repo, move to next repo.
                if location is None:
                    break
            # if reading, only one location is desired.
            if location:
                if not write:
                    # If there is a bypass function for this dataset type, we can't test to see if the object
                    # exists in storage, because the bypass function may not actually use the location
                    # according to the template.
                    # If the user provides a bypass_datasetExists_dataType method use that to set
                    # location.bypass to True/False.
                    # Otherwise execute the bypass function and include its results
                    # in the bypass attribute of the location. The bypass function may fail for any reason,
                    # the most common case being that a file does not exist. If it raises an exception
                    # indicating such, we ignore the bypass function and proceed as though it does not exist.

                    bypass = None
                    if checkExistence:
                        bypass = self._getBypassFunc(location, dataId, "bypass_datasetExists_")
                    if not bypass:
                        bypass = self._getBypassFunc(location, dataId)

                    if bypass:
                        try:
                            location.bypass = bypass()
                        except (NoResults, IOError):
                            self.log.debug("Continuing dataset search while evaluating "
                                           "bypass function for Dataset type:{} Data ID:{} at "
                                           "location {}".format(datasetType, dataId, location))
                    # If a location was found but the location does not exist, keep looking in input
                    # repositories (the registry may have had enough data for a lookup even thought the object
                    # exists in a different repository.)
                    if (isinstance(location, ButlerComposite) or hasattr(location, 'bypass')
                            or location.repository.exists(location)):
                        return location
                else:
                    try:
                        locations.extend(location)
                    except TypeError:
                        locations.append(location)
        if not write:
            return None
        return locations

    @staticmethod
    def _getBypassFunc(location, dataId, bypassBasename="bypass_"):
        """Return a bypass function for the location and dataId, or None

        location: `lsst.daf.persistence.butlerLocation.ButlerLocation`
        dataId: `lsst.daf.persistence.dataId.DataId`
        bypassBasename : `str` the prefix for the bypass method's name (default: "bypass_")
        """
        if not hasattr(location.mapper, bypassBasename + location.datasetType):
            return None

        pythonType = location.getPythonType()
        if pythonType is not None:
            if isinstance(pythonType, str):
                pythonType = doImport(pythonType)
        bypassFunc = getattr(location.mapper, bypassBasename + location.datasetType)
        return lambda: bypassFunc(location.datasetType, pythonType, location, dataId)


class PfsPosixStorage(PosixStorage):
    def butlerLocationExists(self, location):
        """Implementation of PosixStorage.exists for ButlerLocation objects.
        """
        storageName = location.getStorageName()
        if storageName not in ('FitsStorage',
                               'PickleStorage', 'ConfigStorage', 'FitsCatalogStorage',
                               'YamlStorage', 'ParquetStorage', 'MatplotlibStorage'):
            self.log.warn("butlerLocationExists for non-supported storage %s" % location)
            return False

        if hasattr(location, "bypass"):
            return location.bypass

        for locationString in location.getLocations():
            logLoc = LogicalLocation(locationString, location.getAdditionalData()).locString()
            obj = self.instanceSearch(path=logLoc)
            if obj:
                return True
        return False


# Put our modified versions in place of the originals
lsst.daf.persistence.Butler = PfsButler
lsst.daf.persistence.butler.Butler = PfsButler
lsst.daf.persistence.butler.PosixStorage = PfsPosixStorage
lsst.daf.persistence.PosixStorage = PfsPosixStorage
lsst.daf.persistence.posixStorage.PosixStorage = PfsPosixStorage

# Perform the same registrations that were done for the originals
readFormatters = PfsPosixStorage._readFormatters()
readFormatters["FitsStorage"] = readFitsStorage
readFormatters["ParquetStorage"] = readParquetStorage
readFormatters["ConfigStorage"] = readConfigStorage
readFormatters["PickleStorage"] = readPickleStorage
readFormatters["FitsCatalogStorage"] = readFitsCatalogStorage
readFormatters["MatplotlibStorage"] = readMatplotlibStorage
readFormatters["YamlStorage"] = readYamlStorage
writeFormatters = PfsPosixStorage._writeFormatters()
writeFormatters["FitsStorage"] = writeFitsStorage
writeFormatters["ParquetStorage"] = writeParquetStorage
writeFormatters["ConfigStorage"] = writeConfigStorage
writeFormatters["PickleStorage"] = writePickleStorage
writeFormatters["FitsCatalogStorage"] = writeFitsCatalogStorage
writeFormatters["MatplotlibStorage"] = writeMatplotlibStorage
writeFormatters["YamlStorage"] = writeYamlStorage

lsst.daf.persistence.storage.Storage.storages[""] = PfsPosixStorage
lsst.daf.persistence.storage.Storage.storages["file"] = PfsPosixStorage
