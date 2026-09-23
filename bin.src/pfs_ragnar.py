#!/usr/bin/env python

from argparse import ArgumentParser

from lsst.obs.pfs.ragnar import readSweepFile


def main(sweepFile: str, outputFilename: str):
    """Convert a RAGNAR sweep file to a PFS reference line list.

    Parameters
    ----------
    sweepFile : `str`
        Path to a RAGNAR sweep file (``.pkl``).
    outputFilename : `str`
        Path to which to write the resulting line list (text format).
    """
    lines = readSweepFile(sweepFile)
    lines.writeLineList(outputFilename)


if __name__ == "__main__":
    parser = ArgumentParser(description="Convert a RAGNAR sweep file to a PFS reference line list")
    parser.add_argument("sweepFile", help="RAGNAR sweep file (.pkl) to convert")
    parser.add_argument("-o", "--output", required=True, help="Output filename for the reference line list")
    args = parser.parse_args()
    main(args.sweepFile, args.output)
