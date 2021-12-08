import os
import yaml
import sqlite3

# PostgreSQL support
try:
    import psycopg2 as pgsql
    havePgsql = True
except ImportError:
    havePgsql = False

from lsst.daf.persistence.utils import sequencify
from lsst.daf.persistence.registries import SqlRegistry

__all__ = []


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
                    except sqlite3.OperationalError:
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
        except sqlite3.OperationalError:  # try again, with extra checking of the dataId keys
            return self._lookup(lookupProperties, dataId, reference, checkColumns=True)


class PfsSqliteRegistry(PfsSqlRegistry):
    """A SQLite-based registry"""
    placeHolder = "?"  # Placeholder for parameter substitution

    def __init__(self, location):
        """Constructor

        Parameters
        ----------
        location : `str`
            Path to SQLite3 file
        """
        if os.path.exists(location):
            conn = sqlite3.connect(location)
            conn.text_factory = str
            self.root = location
        else:
            conn = None
        SqlRegistry.__init__(self, conn)


class PfsPgsqlRegistry(PfsSqlRegistry):
    """A PostgreSQL-based registry"""
    placeHolder = "%s"

    def __init__(self, location):
        """Constructor

        Parameters
        ----------
        location : `str`
            Path to PostgreSQL configuration file.
        """
        if not havePgsql:
            raise RuntimeError("Cannot use PgsqlRegistry: could not import psycopg2")
        config = self.readYaml(location)
        self._config = config
        conn = pgsql.connect(host=config["host"], port=config["port"], database=config["database"],
                             user=config["user"], password=config["password"])
        self.root = location
        SqlRegistry.__init__(self, conn)

    @staticmethod
    def readYaml(location):
        """Read YAML configuration file

        The YAML configuration file should contain:
        * host : host name for database connection
        * port : port for database connection
        * user : user name for database connection
        * database : database name

        It may also contain:
        * password : password for database connection

        The optional entries are set to `None` in the output configuration.

        Parameters
        ----------
        location : `str`
            Path to PostgreSQL YAML config file.

        Returns
        -------
        config : `dict`
            Configuration
        """
        try:
            # PyYAML >=5.1 prefers a different loader
            loader = yaml.UnsafeLoader
        except AttributeError:
            loader = yaml.Loader
        with open(location) as ff:
            data = yaml.load(ff, Loader=loader)
        requireKeys = set(["host", "port", "database", "user"])
        optionalKeys = set(["password"])
        haveKeys = set(data.keys())
        if haveKeys - optionalKeys != requireKeys:
            raise RuntimeError(
                "PostgreSQL YAML configuration (%s) should contain only %s, and may contain 'password', "
                "but this contains: %s" %
                (location, ",".join("'%s'" % key for key in requireKeys),
                 ",".join("'%s'" % key for key in data.keys()))
            )
        for key in optionalKeys:
            if key not in data:
                data[key] = None

        return data

    def lookup(self, *args, **kwargs):
        try:
            return SqlRegistry.lookup(self, *args, **kwargs)
        except Exception:
            self.conn.rollback()
            raise


import lsst.daf.persistence.registries  # noqa E402: monkey-patching
import lsst.daf.persistence  # noqa E402: monkey-patching
for module in (lsst.daf.persistence, lsst.daf.persistence.registries):
    module.SqlRegistry = PfsSqlRegistry
    module.SqliteRegistry = PfsSqliteRegistry
    module.PgsqlRegistry = PfsPgsqlRegistry
