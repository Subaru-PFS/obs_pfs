import sqlite3
import psycopg2 as pgsql

import lsst.daf.persistence
import lsst.daf.persistence.registries
from lsst.daf.persistence import sequencify
from lsst.daf.persistence.registries import SqlRegistry, SqliteRegistry, PgsqlRegistry


class PfsSqlRegistry(SqlRegistry):
    def _lookup(self, lookupProperties, dataId, reference, checkColumns=False):
        """Perform a lookup in the registry.

        This is the worker code for cls.lookup with the added option of checking
        that all the columns being looked up are in the database.  The classic
        example here is adding a template with an hdu, where the hdu in the dataId
        prevents us looking up e.g. dateObs.  checkColumns results in a performance
        penalty, so is only invoked when a problem in the dataId keys has been seen

        Return values are refined by the values in dataId.
        Returns a list of values that match keys in lookupProperties.
        e.g. if the template is 'raw/raw_v%(visit)d_f%(filter)s.fits.gz', and
        dataId={'visit':1}, and lookupProperties is ['filter'], and the
        filesystem under self.root has exactly one file 'raw/raw_v1_fg.fits.gz'
        then the return value will be [('g',)]

        :param lookupProperties:
        :param dataId: must be a key/value iterable. Keys must be string.
        See `SqlRegistry.lookup` for further details
        :param reference: other data types that may be used to search for values.
        :param checkColumns: if True, check that keys are actually in the registry and ignore them if not
        :return: a list of values that match keys in lookupProperties.
        """
        cmd = "SELECT DISTINCT "
        cmd += ", ".join(lookupProperties)
        cmd += " FROM " + " NATURAL JOIN ".join(reference)
        valueList = []
        if dataId is not None and len(dataId) > 0:
            whereList = []
            for k, v in dataId.items():
                if checkColumns:        # check if k is in registry
                    try:
                        self.conn.cursor().execute(
                            f'SELECT {k} FROM {" NATURAL JOIN ".join(reference)} LIMIT 1')
                    except (pgsql.errors.UndefinedColumn, sqlite3.OperationalError):
                        self.conn.rollback()
                        continue

                if hasattr(k, '__iter__') and not isinstance(k, str):
                    if len(k) != 2:
                        raise RuntimeError("Wrong number of keys for range:%s" % (k,))
                    whereList.append("(%s BETWEEN %s AND %s)" % (self.placeHolder, k[0], k[1]))
                    valueList.append(v)
                else:
                    whereList.append("%s = %s" % (k, self.placeHolder))
                    valueList.append(v)
            cmd += " WHERE " + " AND ".join(whereList)
        cursor = self.conn.cursor()
        cursor.execute(cmd, valueList)
        return [row for row in cursor.fetchall()]

    def lookup(self, lookupProperties, reference, dataId, **kwargs):
        """Perform a lookup in the registry.

        Return values are refined by the values in dataId.
        Returns a list of values that match keys in lookupProperties.
        e.g. if the template is 'raw/raw_v%(visit)d_f%(filter)s.fits.gz', and
        dataId={'visit':1}, and lookupProperties is ['filter'], and the
        filesystem under self.root has exactly one file 'raw/raw_v1_fg.fits.gz'
        then the return value will be [('g',)]

        :param lookupProperties:
        :param dataId: must be a key/value iterable. Keys must be string.
        If value is a string then will look for elements in the repository that match value for value.
        If value is a 2-item iterable then will look for elements in the repository where the value is between
        the values of value[0] and value[1].
        :param reference: other data types that may be used to search for values.
        :param **kwargs: nothing needed for sqlite lookup
        :return: a list of values that match keys in lookupProperties.
        """
        if not self.conn:
            return None

        # input variable sanitization:
        reference = sequencify(reference)
        lookupProperties = sequencify(lookupProperties)

        try:
            return self._lookup(lookupProperties, dataId, reference)
        except (pgsql.errors.UndefinedColumn,
                sqlite3.OperationalError):  # try again, with extra checking of the dataId keys
            self.conn.rollback()
            return self._lookup(lookupProperties, dataId, reference, checkColumns=True)


class PfsSqliteRegistry(PfsSqlRegistry, SqliteRegistry):
    pass


class PfsPgsqlRegistry(PfsSqlRegistry, PgsqlRegistry):
    pass


# Put our modified versions in place of the originals
lsst.daf.persistence.SqlRegistry = PfsSqlRegistry
lsst.daf.persistence.SqliteRegistry = PfsSqliteRegistry
lsst.daf.persistence.PgsqlRegistry = PfsPgsqlRegistry
lsst.daf.persistence.registries.SqlRegistry = PfsSqlRegistry
lsst.daf.persistence.registries.SqliteRegistry = PfsSqliteRegistry
lsst.daf.persistence.registries.PgsqlRegistry = PfsPgsqlRegistry
