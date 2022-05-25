#!/usr/bin/env python

from argparse import ArgumentParser
from glob import glob
from typing import Iterable
from lsst.obs.pfs.checkRawHeader import checkRawHeader
from astropy import log


def main(filenames: Iterable[str], allowFix: bool = False):
    """Main driver to check/fix PFS headers

    Parameters
    ----------
    filenames : Iterable[str]
        List of files (or globs) for which to check/fix headers.
    allowFix : bool, optional
        Allow headers to be fixed?
    """
    numGood = 0
    numBad = 0
    for gg in filenames:
        for fn in glob(gg):
            try:
                checkRawHeader(fn, allowFix)
                numGood += 1
            except Exception as exc:
                log.warning(exc)
                numBad += 1
    log.info(f"Processed {numGood} good and {numBad} bad files.")


if __name__ == "__main__":
    parser = ArgumentParser(description="Check raw PFS files for header compliance")
    parser.add_argument("--fix", default=False, action="store_true", help="Fix bad/missing headers?")
    parser.add_argument("filenames", nargs="+", help="Names of files (or glob) to check")
    args = parser.parse_args()
    main(args.filenames, args.fix)
