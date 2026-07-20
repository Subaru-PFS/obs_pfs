#!/usr/bin/env python
"""Copy or link datasets between collections in a butler repository.

A dataset's RUN collection is immutable, so there is no butler move or rename,
and ``butler transfer-datasets`` only transfers between distinct repositories.
This script copies each dataset's file and ingests it under a fresh dataset ID in
``--output-run``, which is the only way to re-home a dataset within one
repository.

Its intended use is promoting a calibration from a personal collection into the
shared one, e.g. from ``u/<user>/calib/<ticket>/...`` to ``PFS/calib/<ticket>/...``:

    copyDatasets.py /work/datastore \\
        --input u/cpl/calib/PIPE2D-1664/irp4/nirDark_irp4Gen.20260709a/put \\
        --dataset-type nirDark_irp4 \\
        --output-run PFS/calib/PIPE2D-1664/irp4/nirDark_irp4Gen.20260709a/put \\
        --certify PFS/calib/PIPE2D-1664/irp4/nirDark_irp4Gen.20260709a

With ``--certify`` and no ``--begin-date``, the validity period is inherited from
the source certification, so it need not be restated.

Pass ``--no-copy`` to certify the datasets where they already are, without
duplicating the pixels. That leaves the destination CALIBRATION collection
referencing the source RUN, which must then not be pruned.
"""

from __future__ import annotations

import logging
from argparse import ArgumentParser

from lsst.daf.butler import Butler

from lsst.obs.pfs import datasetCopy


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("repo", help="Path to the butler repository")
    parser.add_argument("--input", nargs="+", required=True, dest="inputCollections",
                        help="Collection(s) to take the datasets from")
    parser.add_argument("--dataset-type", nargs="+", required=True, dest="datasetTypes",
                        help="Dataset type name(s) to copy")
    parser.add_argument("--output-run", default=None, dest="outputRun",
                        help="RUN collection to copy the datasets into; created if absent")
    parser.add_argument("--no-copy", action="store_true", dest="noCopy",
                        help="Leave the datasets in their existing RUN instead of copying "
                             "them. Requires --certify, and makes the certified collection "
                             "depend on the source RUN surviving")
    parser.add_argument("--certify", default=None, dest="calibCollection",
                        help="CALIBRATION collection to certify the results into")
    parser.add_argument("--begin-date", default=None, dest="beginDate",
                        help="ISO-8601 TAI start of the validity period "
                             "(default: inherited from the source certification)")
    parser.add_argument("--end-date", default=None, dest="endDate",
                        help="ISO-8601 TAI end of the validity period (default: open-ended)")
    parser.add_argument("--where", default="",
                        help="Butler query expression restricting which datasets are copied")
    parser.add_argument("--transfer", default="copy", dest="transferMode",
                        choices=["copy", "hardlink", "symlink", "relsymlink", "link",
                                 "move", "direct", "auto"],
                        help="How to transfer the file artifact. Only 'copy' (the "
                             "default) leaves the destination independent of the source")
    parser.add_argument("--skip-existing", action="store_true", dest="skipExisting",
                        help="Skip datasets already present in the output run, so an "
                             "interrupted copy can be resumed")
    parser.add_argument("--allow-intermediates", action="store_true", dest="allowIntermediates",
                        help="Permit copying an intermediate product such as rawISRCube, "
                             "which is refused by default")
    parser.add_argument("--dry-run", action="store_true", dest="dryRun",
                        help="Report what would be done without writing anything")
    args = parser.parse_args()

    if args.noCopy:
        if args.outputRun is not None:
            parser.error("--no-copy leaves the datasets in place; do not pass --output-run")
        if args.calibCollection is None:
            parser.error("--no-copy does nothing without --certify")
    elif args.outputRun is None:
        parser.error("--output-run is required unless --no-copy is given")

    butler = Butler.from_config(args.repo, writeable=not args.dryRun)
    datasetCopy.transfer(butler, args.datasetTypes, args.inputCollections,
                         outputRun=args.outputRun, where=args.where,
                         calibCollection=args.calibCollection,
                         beginDate=args.beginDate, endDate=args.endDate,
                         skipExisting=args.skipExisting, transferMode=args.transferMode,
                         allowIntermediates=args.allowIntermediates, dryRun=args.dryRun)


if __name__ == "__main__":
    main()
