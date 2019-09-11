import os
from collections.abc import Iterable
from collections import defaultdict

import yaml


class HeaderFixDatabase:
    """Database of header fixes to apply

    Header fixes are a set of `dict`s indexed by ``visit``. These get written
    in YAML format for the astro_metadata_translator package to read.

    The standard set of header fixes is defined in the ``addStandardFixes``
    method. Further fixes can be added using the ``add`` method. New fixes can
    clobber (or supplement) previous fixes.
    """
    def __init__(self):
        self.fixes = defaultdict(dict)  # visit --> header fixes to apply
        self.addStandardFixes()

    def add(self, visits, **fixes):
        """Add a set of fixes for a list of visits

        Parameters
        ----------
        visits : iterable or `int`
            List of visits (or a single visit) for which the fixes apply.
        **fixes : `dict` mapping `str` to `str`/`int`/`float`
            Header fixes to apply.
        """
        for vv in visits if isinstance(visits, Iterable) else [visits]:
            self.fixes[vv].update(fixes)

    def write(self, path):
        """Write the database

        The database is written as one YAML file for each visit; this is the
        format mandated by astro_metadata_translator.

        Parameters
        ----------
        path : `str`
            Path to which to write correction files.
        """
        for visit, fixes in self.fixes.items():
            filename = os.path.join(path, "PFS-%s.yaml" % (visit,))
            with open(filename, "w") as ff:
                yaml.dump(fixes, ff)

    def addStandardFixes(self):
        """Add the standard set of fixes

        This is the official list of header fixes to apply.
        """
        pass


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Build header correction files")
    parser.add_argument("path", help="Path to which to write correction files")
    args = parser.parse_args()
    HeaderFixDatabase().write(args.path)
