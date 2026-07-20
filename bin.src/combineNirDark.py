#!/usr/bin/env python
"""Combine N raw NIR dark ramps into a single master nirDark ramp.

Reads the per-visit ``rawISRCube`` ramps for the given dark visits from the
input collection, median-combines them read-by-read (correcting per-read offsets
between ramps), puts the result into a RUN collection and certifies it into a
CALIBRATION collection, so that the pipeline selects it by observation date. The
output dataset type is chosen from the IRP ratio of the input ramps: ``nirDark``
for the default (IRP1) or ``nirDark_irp<N>`` (e.g. ``nirDark_irp4``) otherwise.

Collections follow the DMTN-222 convention (https://dmtn-222.lsst.io), as the
pipeline-generated ``dark`` and ``bias`` do. Under the ``--output`` base, with
``{iteration}`` a ``YYYYMMDDv`` date string:

- ``{output}/{ticket}/{tag}/{product}.{iteration}`` -- CALIBRATION, certified
  into, and what the pipeline selects the dark from;
- ``{output}/{ticket}/{tag}/{product}Gen.{iteration}`` -- CHAINED, gathering this
  generation's run and its input;
- ``{output}/{ticket}/{tag}/{product}Gen.{iteration}/{YYYYMMDDTHHMMSSZ}`` -- the
  RUN holding the datasets. Each attempt gets its own timestamped RUN, so a retry
  after a failed write does not collide with the previous one.

Only those three collections are written. Shared chains such as ``PFS/calib`` are
maintained by hand and are never touched; promoting a dark into one is a separate
step (see ``copyDatasets.py``).

By default all four NIR spectrographs are reduced in parallel, one worker
process each. The output collections and dataset type are named and registered
before any combine starts, and each combined dark is written to ``--save-dir``
as a plain FITS file before it is put to the butler.
"""

from __future__ import annotations

import concurrent.futures
import datetime
import logging
import os
from argparse import ArgumentParser

from lsst.obs.pfs import nirSuperdark


def parseVisitRange(token: str) -> list[int]:
    """Expand a single visit token into a list of visit ids.

    Accepts a plain integer (``"1234"``) or an LSST-style inclusive range
    ``"BEGIN..END"`` with an optional step ``"BEGIN..END:STEP"`` (e.g.
    ``"144784..144808"`` or ``"144784..144808:2"``).
    """
    if ".." not in token:
        return [int(token)]
    rangePart, _, stepPart = token.partition(":")
    try:
        beginStr, endStr = rangePart.split("..")
        begin, end = int(beginStr), int(endStr)
        step = int(stepPart) if stepPart else 1
    except ValueError:
        raise ValueError(f"malformed visit range {token!r}; expected BEGIN..END[:STEP]")
    if step <= 0:
        raise ValueError(f"step must be positive in visit range {token!r}")
    if end < begin:
        raise ValueError(f"range end precedes begin in visit range {token!r}")
    return list(range(begin, end + 1, step))


def parseVisits(tokens: list[str]) -> list[int]:
    """Expand a list of visit tokens (ints and/or ranges) into visit ids."""
    visits: list[int] = []
    for token in tokens:
        visits.extend(parseVisitRange(token))
    return visits


def reduceSpectrographs(repo: str, inputRun: str, outputBase: str,
                        spectrographs: list[int], visits: list[int],
                        ticket: str, tag: str, iteration: str, timestamp: str,
                        startDate=None,
                        processes: int | None = None,
                        saveDir: str | None = None) -> None:
    """Reduce each spectrograph's NIR dark, in parallel across worker processes.

    Each spectrograph is an independent reduction (its own dataId and output
    dataset), so they run concurrently. ``processes`` defaults to one worker per
    spectrograph; pass a smaller number to cap memory use, since each reduction
    holds the full ramp stack in memory. Failures are collected so that a problem
    with one spectrograph does not abandon the others; a `RuntimeError`
    summarizing them is raised once all have been attempted.

    The collections, the dataset type, the RUN timestamp, and the validity start
    date are resolved once, here in the parent, before any worker starts: a combine
    takes hours, so it must not be a butler misconfiguration that discovers itself
    at the end. Doing it once also keeps the four workers from racing to register
    the same things, gives the whole dark set a single validity start rather than
    one per detector's timestamp, and puts all four spectrographs in one RUN.
    """
    dataIds = [dict(instrument="PFS", arm="n", spectrograph=s) for s in spectrographs]
    processes = len(dataIds) if processes is None else processes
    processes = max(1, min(processes, len(dataIds)))

    plan = nirSuperdark.preflight(inputRun, outputBase, dataIds[0], visits,
                                  ticket, tag, iteration, timestamp, repo_path=repo)
    if startDate is None:
        startDate = plan.startDate
    logger = logging.getLogger(__name__)
    logger.info("will write %s for spectrograph(s) %s", plan.datasetType,
                ", ".join(str(s) for s in spectrographs))
    logger.info("  run       %s", plan.outputRun)
    logger.info("  chain     %s", plan.genCollection)
    logger.info("  certify   %s from %s", plan.calibCollection, startDate)
    if saveDir is not None:
        os.makedirs(saveDir, exist_ok=True)

    kwargs = dict(startDate=startDate, repo_path=repo,
                  saveDir=saveDir, calibCollection=plan.calibCollection)
    outputRun = plan.outputRun

    errors = []
    if processes == 1:
        for dataId in dataIds:
            try:
                nirSuperdark.processMasterDark(inputRun, outputRun, dataId, visits,
                                               **kwargs)
            except Exception as exc:
                errors.append((dataId["spectrograph"], exc))
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as pool:
            futures = {}
            for dataId in dataIds:
                future = pool.submit(nirSuperdark.processMasterDark, inputRun,
                                     outputRun, dataId, visits, **kwargs)
                futures[future] = dataId["spectrograph"]
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    errors.append((futures[future], exc))

    if errors:
        summary = "; ".join(f"n{s}: {exc}" for s, exc in sorted(errors))
        raise RuntimeError(
            f"dark reduction failed for {len(errors)} of {len(dataIds)} "
            f"spectrograph(s): {summary}")


def main():
    # Configure logging before the worker pool is created so forked workers
    # inherit it and their per-read progress lines reach the terminal.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("repo", help="Path to the butler repository")
    parser.add_argument("--input", required=True,
                        help="Input collection holding the per-visit rawISRCubes")
    parser.add_argument("--output", required=True,
                        help="Collection base the DMTN-222 output collections are "
                             "composed under, e.g. u/<user>/calib or PFS/calib")
    parser.add_argument("--ticket", required=True, help="Ticket name, e.g. PIPE2D-1664")
    parser.add_argument("--tag", required=True, help="Tag/label for this calibration set")
    parser.add_argument("--iteration", default=None,
                        help="DMTN-222 rerun iteration (YYYYMMDDv); default: today + 'a'")
    parser.add_argument("--run-timestamp", default=None, dest="timestamp",
                        help="Timestamp naming the output RUN within the generation "
                             "chain (YYYYMMDDTHHMMSSZ); default: now, in UTC")
    parser.add_argument("--save-dir", default=".", dest="saveDir",
                        help="Directory to write the combined dark to as a plain FITS "
                             "file before putting it to the butler, so a butler failure "
                             "does not discard hours of work (default: current directory)")
    parser.add_argument("--spectrograph", type=int, nargs="+", default=[1, 2, 3, 4],
                        help="Spectrograph number(s) to process (default: all four)")
    parser.add_argument("--visits", nargs="+", required=True,
                        help="Dark visits to combine: integers and/or LSST-style "
                             "inclusive ranges BEGIN..END[:STEP] (e.g. 144784..144808)")
    parser.add_argument("--start-date", default=None, dest="startDate",
                        type=datetime.datetime.fromisoformat,
                        help="ISO-8601 start of the calibration's open-ended validity "
                             "period (default: the date of the first visit)")
    parser.add_argument("--processes", type=int, default=None,
                        help="Number of parallel worker processes "
                             "(default: one per spectrograph)")
    args = parser.parse_args()

    iteration = args.iteration
    if iteration is None:
        iteration = nirSuperdark.defaultIteration()
    timestamp = args.timestamp
    if timestamp is None:
        timestamp = nirSuperdark.defaultTimestamp()

    reduceSpectrographs(args.repo, args.input, args.output, args.spectrograph,
                        parseVisits(args.visits), args.ticket, args.tag, iteration,
                        timestamp, startDate=args.startDate, processes=args.processes,
                        saveDir=args.saveDir)


if __name__ == "__main__":
    main()
