#!/usr/bin/env python

from argparse import ArgumentParser

from lsst.obs.pfs.ragnar import readSweepFile


def main(sweepFile: str, outputFilename: str, arm: str | None = None):
    """Convert a RAGNAR sweep file to a PFS reference line list.

    Parameters
    ----------
    sweepFile : `str`
        Path to a RAGNAR sweep file (``.pkl``).
    outputFilename : `str`
        Path to which to write the resulting line list (text format).
    arm : `str`, optional
        Spectrograph arm (``b``, ``r``, ``n``) for which to select lines;
        required for unified (multi-arm) sweep files.
    """
    lines = readSweepFile(sweepFile, arm)
    lines.writeLineList(outputFilename)


if __name__ == "__main__":
    parser = ArgumentParser(description="Convert a RAGNAR sweep file to a PFS reference line list")
    parser.add_argument("sweepFile", help="RAGNAR sweep file (.pkl) to convert")
    parser.add_argument("-o", "--output", required=True, help="Output filename for the reference line list")
    parser.add_argument(
        "--arm", choices=("b", "r", "n"), help="Arm to select (required for unified sweep files)"
    )
    args = parser.parse_args()
    main(args.sweepFile, args.output, args.arm)
