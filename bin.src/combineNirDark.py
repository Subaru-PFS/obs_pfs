#!/usr/bin/env python
"""Combine N raw NIR dark ramps into a single master nirDark ramp.

Reads the per-visit ``rawISRCube`` ramps that ``makeNirRawCubes.py`` wrote to the
input collection, median-combines them read-by-read (correcting per-read offsets
between ramps), puts the result into a RUN collection and certifies it into a
CALIBRATION collection, so that the pipeline selects it by observation date. The
output dataset type is chosen from the IRP ratio of the input ramps: ``nirDark``
for the default (IRP1) or ``nirDark_irp<N>`` (e.g. ``nirDark_irp4``) otherwise.

Usually all that is needed is the input, which ``makeNirRawCubes.py`` prints::

    combineNirDark.py /work/datastore --input u/cpl/calib/PIPE2D-1664/irp1/scratchCubes

``--output``, ``--ticket`` and ``--tag`` default to the parts of an input named
``{output}/{ticket}/{tag}/scratchCubes`` (or a RUN within it), and ``--visits``
to every visit the input holds a ``rawISRCube`` for. Given ``--visits``, each
must have a ``rawISRCube`` for every spectrograph, or the combine stops before
starting; ``--skip-missing`` combines the rest instead, as it must when a visit
could not be processed (e.g. one with no SpS exposure).

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

The dark is valid from the date of its first visit (``--start-date``), and
open-ended unless ``--end-date`` is given. Close it when regenerating a dark for
older data, so it does not apply to data taken after it.

By default all four NIR spectrographs are reduced in parallel, one worker
process each. The visits, output collections and dataset type are resolved and
registered before any combine starts, and each combined dark is written to
``--save-dir`` as a plain FITS file before it is put to the butler.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import logging
import multiprocessing
import os
import time

from lsst.obs.pfs import nirSuperdark


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Keep the description's layout, and report each option's default."""

    def _get_help_string(self, action):
        text = action.help
        default = action.default
        if (default is not None and default is not argparse.SUPPRESS
                and not isinstance(default, bool) and "%(default)" not in text):
            text += " (default: %(default)s)"
        return text


def helpDescription(doc: str) -> str:
    """The module docstring, with the reStructuredText markup dropped."""
    return doc.replace("``", "").replace("::\n", ":\n")


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
    """Expand a list of visit tokens (ints and/or ranges) into visit ids.

    Each token may itself be a comma-separated list, e.g.
    ``"121023..121030,121032..121048"``.
    """
    visits: list[int] = []
    for token in tokens:
        for field in token.split(","):
            if field.strip():
                visits.extend(parseVisitRange(field.strip()))
    return visits


def reduceSpectrographs(repo: str, inputRun: str, outputBase: str,
                        spectrographs: list[int], visits: list[int] | None,
                        ticket: str, tag: str, iteration: str, timestamp: str,
                        startDate=None,
                        processes: int | None = None,
                        saveDir: str | None = None,
                        endDate=None,
                        skipMissing: bool = False) -> None:
    """Reduce each spectrograph's NIR dark, in parallel across worker processes.

    Each spectrograph is an independent reduction (its own dataId and output
    dataset), so they run concurrently. ``processes`` is the total number of
    worker processes, each reducing one spectrograph at a time; it defaults to
    one per spectrograph. Pass a smaller number to cap memory use, since each
    reduction holds the full ramp stack in memory. Failures are collected so that
    a problem with one spectrograph does not abandon the others; a `RuntimeError`
    summarizing them is raised once all have been attempted.

    The visits, the collections, the dataset type, the RUN timestamp, and the
    validity range are resolved once, here in the parent, before any worker
    starts: a combine takes hours, so it must not be a missing input or a butler
    misconfiguration that discovers itself at the end. Doing it once also keeps
    the workers from racing to register the same things, gives the whole dark set
    a single validity start rather than one per detector's timestamp, and puts
    all the spectrographs in one RUN.
    """
    logger = logging.getLogger(__name__)
    checkValidity(startDate, endDate)
    plan = nirSuperdark.preflight(inputRun, outputBase, spectrographs, visits,
                                  ticket, tag, iteration, timestamp, repo_path=repo,
                                  skipMissing=skipMissing)
    if startDate is None:
        startDate = plan.startDate
    checkValidity(startDate, endDate)
    cams = ", ".join(f"n{s}" for s in spectrographs)
    logger.info("will write %s for %s", plan.datasetType, cams)
    logger.info("  run       %s", plan.outputRun)
    logger.info("  chain     %s", plan.genCollection)
    logger.info("  certify   %s from %s to %s", plan.calibCollection, startDate,
                "open end" if endDate is None else endDate)
    if saveDir is not None:
        os.makedirs(saveDir, exist_ok=True)

    kwargs = dict(startDate=startDate, endDate=endDate, repo_path=repo,
                  saveDir=saveDir, calibCollection=plan.calibCollection)
    outputRun = plan.outputRun
    processes = len(spectrographs) if processes is None else processes
    processes = max(1, min(processes, len(spectrographs)))

    t0 = time.time()
    errors = []
    done = []

    def report(spectrograph, exc=None):
        """Log one spectrograph's outcome as it arrives."""
        done.append(spectrograph)
        progress = (f"{(time.time() - t0)/60:.1f} min "
                    f"({len(done)}/{len(spectrographs)} complete)")
        if exc is None:
            logger.info("n%d: done in %s", spectrograph, progress)
        else:
            errors.append((spectrograph, exc))
            logger.error("n%d: FAILED after %s: %s", spectrograph, progress, exc)

    logger.info("combining %s with %d worker process(es)", cams, processes)
    if processes == 1:
        for spectrograph in spectrographs:
            logger.info("n%d: starting", spectrograph)
            try:
                nirSuperdark.processMasterDark(inputRun, outputRun,
                                               nirSuperdark.darkDataId(spectrograph),
                                               plan.visits[spectrograph], **kwargs)
            except Exception as exc:
                report(spectrograph, exc)
            else:
                report(spectrograph)
    else:
        # A forkserver start method (rather than the default fork) keeps the
        # workers from being forked out of this multi-threaded parent, where an
        # inherited lock could deadlock an hours-long combine. Its workers start
        # from a clean interpreter, so the initializer restores their logging.
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=processes,
                mp_context=multiprocessing.get_context("forkserver"),
                initializer=nirSuperdark.configureLogging) as pool:
            futures = {}
            for spectrograph in spectrographs:
                future = pool.submit(nirSuperdark.processMasterDark, inputRun, outputRun,
                                     nirSuperdark.darkDataId(spectrograph),
                                     plan.visits[spectrograph], **kwargs)
                futures[future] = spectrograph
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    report(futures[future], exc)
                else:
                    report(futures[future])

    succeeded = [s for s in spectrographs if s not in dict(errors)]
    logger.info("%s for %s certified into %s in %.1f min%s", plan.datasetType,
                ", ".join(f"n{s}" for s in succeeded) or "no spectrograph",
                plan.calibCollection, (time.time() - t0)/60,
                f"; FAILED: {', '.join(f'n{s}' for s, _ in sorted(errors))}" if errors else "")
    if errors:
        summary = "; ".join(f"n{s}: {exc}" for s, exc in sorted(errors))
        raise RuntimeError(
            f"dark reduction failed for {len(errors)} of {len(spectrographs)} "
            f"spectrograph(s): {summary}")


def checkValidity(startDate, endDate) -> None:
    """Refuse a validity range that ends before it starts."""
    if startDate is None or endDate is None:
        return
    if nirSuperdark.asTime(endDate) <= nirSuperdark.asTime(startDate):
        raise ValueError(f"validity end {endDate} is not after its start {startDate}")


def main():
    nirSuperdark.configureLogging()
    parser = argparse.ArgumentParser(description=helpDescription(__doc__),
                                     formatter_class=HelpFormatter)
    parser.add_argument("repo", help="Path to the butler repository")
    parser.add_argument("--input", required=True,
                        help="Input collection holding the per-visit rawISRCubes: "
                             "normally the scratchCubes collection makeNirRawCubes wrote")
    parser.add_argument("--output", default=None,
                        help="Collection base the DMTN-222 output collections are "
                             "composed under, e.g. u/<user>/calib or PFS/calib "
                             "(default: inferred from --input)")
    parser.add_argument("--ticket", default=None,
                        help="Ticket name, e.g. PIPE2D-1664 (default: inferred from --input)")
    parser.add_argument("--tag", default=None,
                        help="Tag/label for this calibration set "
                             "(default: inferred from --input)")
    parser.add_argument("--iteration", default=None,
                        help="DMTN-222 rerun iteration (YYYYMMDDv); default: today + 'a'")
    parser.add_argument("--run-timestamp", default=None, dest="timestamp",
                        help="Timestamp naming the output RUN within the generation "
                             "chain (YYYYMMDDTHHMMSSZ); default: now, in UTC")
    parser.add_argument("--save-dir", default=".", dest="saveDir",
                        help="Directory to write the combined dark to as a plain FITS "
                             "file before putting it to the butler, so a butler failure "
                             "does not discard hours of work")
    parser.add_argument("--spectrograph", type=int, nargs="+", default=[1, 2, 3, 4],
                        help="Spectrograph number(s) to process")
    parser.add_argument("--visits", nargs="+", default=None,
                        help="Dark visits to combine: integers and/or LSST-style "
                             "inclusive ranges BEGIN..END[:STEP], space- or "
                             "comma-separated (e.g. 121023..121030,121032..121048); "
                             "default: every visit with a rawISRCube in --input")
    parser.add_argument("--skip-missing", action="store_true", dest="skipMissing",
                        help="Combine the rest when some --visits have no rawISRCube, "
                             "rather than stopping")
    parser.add_argument("--start-date", default=None, dest="startDate",
                        type=datetime.datetime.fromisoformat,
                        help="ISO-8601 start of the calibration's validity period; "
                             "default: the date of the first visit")
    parser.add_argument("--end-date", default=None, dest="endDate",
                        type=datetime.datetime.fromisoformat,
                        help="ISO-8601 end of the calibration's validity period; "
                             "default: open-ended")
    parser.add_argument("-j", "--processes", type=int, default=None,
                        help="Total number of worker processes, each reducing one "
                             "spectrograph at a time, so more than the number of "
                             "spectrographs gains nothing; default: one per spectrograph")
    args = parser.parse_args()

    inferred = nirSuperdark.parseScratchCollection(args.input)
    names = dict(output=args.output, ticket=args.ticket, tag=args.tag)
    if inferred is not None:
        for key, value in zip(names, inferred):
            if names[key] is None:
                names[key] = value
    unset = [f"--{key}" for key, value in names.items() if value is None]
    if unset:
        parser.error(f"{', '.join(unset)} cannot be inferred from --input {args.input} "
                     f"(not named {{output}}/{{ticket}}/{{tag}}/"
                     f"{nirSuperdark.SCRATCH_COLLECTION}); pass them explicitly")

    iteration = args.iteration
    if iteration is None:
        iteration = nirSuperdark.defaultIteration()
    timestamp = args.timestamp
    if timestamp is None:
        timestamp = nirSuperdark.defaultTimestamp()
    visits = None if args.visits is None else parseVisits(args.visits)

    reduceSpectrographs(args.repo, args.input, names["output"], args.spectrograph,
                        visits, names["ticket"], names["tag"], iteration,
                        timestamp, startDate=args.startDate, endDate=args.endDate,
                        processes=args.processes, saveDir=args.saveDir,
                        skipMissing=args.skipMissing)


if __name__ == "__main__":
    main()
