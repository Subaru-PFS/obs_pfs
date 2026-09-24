#!/usr/bin/env python
"""Generate the raw NIR ramp cubes (``rawISRCube``) that combineNirDark consumes.

Runs the ``isr`` step of ``reduceExposure.yaml`` over the given NIR visits with
the corrections that a dark must *not* have applied: no dark subtraction, no
flat, no linearity, and no CR/glitch correction. The point of a
nirDark is to capture the full raw instrument dark signature, which is subtracted
from an exposure ramp before any corrections, so the ramps it is built from must
themselves be raw. Per-ramp CR repair is redundant besides: combineNirDark's
per-read median across ramps rejects temporally independent events (CRs).

The overridden configuration, and the default it replaces:

=========================  ========  =====
Config                     Default   Here
=========================  ========  =====
isr.doDark                 True      False
isr.doFlat                 True      False
isr.h4.doWriteRawCube      False     True
isr.h4.doLinearize         True      False
isr.h4.doCR                True      False
isr.h4.quickCDS            False     False
=========================  ========  =====

Reference-pixel handling is selectable, and must match the exposures the dark
will be subtracted from: ``--irp-filter`` sets ``isr.h4.IRPfilter`` (default -1
= per-channel, per-column median, as ISR itself defaults to; 0 = no smoothing;
odd 15..31 = Hann-smoothed), and ``--no-irp`` sets ``isr.h4.useIRP=False`` to
bypass the IRP planes entirely and border-correct instead.

The input defaults to ``PFS/defaults``, which chains the raws (``PFS/raw/sps``),
the pfsConfigs and the calibrations (``PFS/calib``) the ISR step needs;
``PFS/defaults`` is appended to any ``--input`` given.

Like ``combineNirDark.py``, ``--output`` is a collection base beneath which the
DMTN-222 layout is composed: the cubes go to a ``scratchCubes`` output collection
at ``{output}/{ticket}/{tag}/scratchCubes``. They are an intermediate product, not
a calibration, so they get no ``Gen.{iteration}`` name and are not certified.
Pipetask makes that a CHAINED collection and writes a timestamped RUN inside it,
so re-running does not collide.

The output collection is reported before pipetask starts, and once it succeeds
the ``combineNirDark.py`` command that combines the cubes is printed, ready to
paste. It passes ``--skip-missing``, because any requested visit without a cube
by then is one pipetask had nothing to process for (e.g. no SpS exposure).

Example::

    makeNirRawCubes.py /work/datastore \\
        --output u/cpl/calib --ticket PIPE2D-1664 --tag irp1 \\
        --visits 144587..144636

writes to ``u/cpl/calib/PIPE2D-1664/irp1/scratchCubes``, with one pipetask
process per NIR camera.

``--show config`` (or ``--show uri``) is passed through to pipetask, and reports
the resolved configuration without running anything.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys

from lsst.obs.pfs import nirSuperdark
from lsst.utils import getPackageDir

# The corrections a raw dark ramp must not have had applied, independent of how
# the reference pixels are handled. quickCDS is already the default, but is set
# explicitly: these ramps are defined by their processing.
ISR_CONFIG = (
    "isr:doDark=False",
    "isr:doFlat=False",
    "isr:h4.doWriteRawCube=True",
    "isr:h4.quickCDS=False",
    "isr:h4.doLinearize=False",
    "isr:h4.doCR=False",
)

# The default IRP filter, matching the ISR default the exposures are processed
# with: -1 = per-channel, per-column median; 0 = no smoothing; odd 15..31 =
# Hann-smoothed. `useIRP=False` bypasses IRP for border correction.
DEFAULT_IRP_FILTER = -1

DEFAULT_INPUTS = ("PFS/defaults",)

# The NIR spectrographs, one pipetask process each by default.
NIR_SPECTROGRAPHS = (1, 2, 3, 4)


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


def irpConfig(irpFilter: int, useIRP: bool) -> tuple[str, ...]:
    """The reference-pixel config overrides.

    ``useIRP`` is always stated explicitly; ``IRPfilter`` is emitted only when IRP
    is in use, since it does nothing when the reference planes are bypassed.
    """
    overrides = [f"isr:h4.useIRP={useIRP}"]
    if useIRP:
        overrides.append(f"isr:h4.IRPfilter={irpFilter}")
    return tuple(overrides)


def visitQuery(visits: list[str], spectrographs: list[int] | None = None) -> str:
    """Build the pipetask data-id expression selecting the NIR dark visits.

    ``visits`` are passed through to the butler query verbatim, so both plain
    integers and LSST-style ranges (``144587..144636``, optionally ``:STEP``)
    work, as the butler understands both.
    """
    query = f"visit in ({', '.join(visits)}) and arm='n'"
    if spectrographs:
        query += f" and spectrograph in ({', '.join(str(s) for s in spectrographs)})"
    return query


def pipelinePath() -> str:
    """The ``isr`` subset of drp_stella's reduceExposure pipeline."""
    return f"{getPackageDir('drp_stella')}/pipelines/reduceExposure.yaml#isr"


def buildCommand(repo: str, inputs: list[str], output: str, visits: list[str],
                 spectrographs: list[int] | None = None,
                 processes: int = 1, logLevel: str | None = None,
                 irpFilter: int = DEFAULT_IRP_FILTER, useIRP: bool = True,
                 config: list[str] | None = None,
                 show: list[str] | None = None) -> list[str]:
    """Assemble the ``pipetask run`` command line."""
    command = ["pipetask", "--long-log"]
    if logLevel:
        command += ["--log-level", logLevel]
    command += [
        "run",
        "-b", repo,
        "-i", ",".join(inputs),
        "-o", output,
        "-p", pipelinePath(),
        "--fail-fast",
        "-j", str(processes),
        "-d", visitQuery(visits, spectrographs),
    ]
    for override in ISR_CONFIG + irpConfig(irpFilter, useIRP) + tuple(config or ()):
        command += ["-c", override]
    for item in show or ():
        command += ["--show", item]
    return command


def combineCommand(repo: str, output: str, visits: list[str],
                   spectrographs: list[int] | None,
                   outputBase: str, ticket: str, tag: str) -> list[str]:
    """The ``combineNirDark.py`` command that combines the cubes in ``output``.

    combineNirDark infers ``--output``, ``--ticket`` and ``--tag`` from a
    ``scratchCubes`` collection, so they are only spelled out for an overridden
    output collection. The visits are passed on because the chain may also hold
    cubes from earlier runs for other visits.
    """
    command = ["combineNirDark.py", repo, "--input", output]
    if nirSuperdark.parseScratchCollection(output) is None:
        command += ["--output", outputBase, "--ticket", ticket, "--tag", tag]
    command += ["--visits", ",".join(visits), "--skip-missing"]
    if spectrographs:
        command += ["--spectrograph", *(str(s) for s in spectrographs)]
    return command


def main():
    parser = argparse.ArgumentParser(description=helpDescription(__doc__),
                                     formatter_class=HelpFormatter)
    parser.add_argument("repo", help="Path to the butler repository")
    parser.add_argument("--input", nargs="+", default=[], dest="inputs",
                        help="Input collection(s), to which PFS/defaults is appended "
                             "if absent; default: PFS/defaults alone")
    parser.add_argument("--output", required=True,
                        help="Collection base the output is composed under, "
                             "e.g. u/<user>/calib or PFS/calib")
    parser.add_argument("--ticket", required=True, help="Ticket name, e.g. PIPE2D-1664")
    parser.add_argument("--tag", required=True, help="Tag/label for this set of cubes")
    parser.add_argument("--output-collection", default=None, dest="outputCollection",
                        help="Write here instead of {output}/{ticket}/{tag}/scratchCubes")
    parser.add_argument("--visits", nargs="+", required=True,
                        help="Dark visits: integers and/or LSST-style inclusive ranges "
                             "BEGIN..END[:STEP] (e.g. 144587..144636)")
    parser.add_argument("--spectrograph", type=int, nargs="+", default=None,
                        help="Spectrograph number(s); default: all")
    parser.add_argument("-j", "--processes", type=int, default=None,
                        help="Number of pipetask worker processes; default: one per "
                             "spectrograph processed")
    parser.add_argument("--irp-filter", type=int, default=DEFAULT_IRP_FILTER, dest="irpFilter",
                        help="IRP filter: -1=per-channel, per-column median, 0=no smoothing, "
                             "odd 15..31=Hann-smoothed. Ignored with --no-irp")
    parser.add_argument("--no-irp", action="store_false", dest="useIRP",
                        help="Bypass the interleaved reference pixels entirely and "
                             "border-correct instead (isr:h4.useIRP=False)")
    parser.add_argument("--log-level", default=None, dest="logLevel",
                        help="pipetask --log-level, e.g. .=DEBUG")
    parser.add_argument("-c", "--config", action="append", default=None,
                        help="Extra pipetask config override, e.g. isr:h4.doCR=True. "
                             "Repeatable; applied after the raw-cube overrides")
    parser.add_argument("--show", action="append", default=None,
                        help="Passed to pipetask, e.g. 'config' or 'uri'. Reports the "
                             "resolved configuration without processing anything. "
                             "Repeatable")
    parser.add_argument("--dry-run", action="store_true", dest="dryRun",
                        help="Print the pipetask command without running it")
    args = parser.parse_args()

    inputs = list(args.inputs)
    for default in DEFAULT_INPUTS:
        if default not in inputs:
            inputs.append(default)

    output = args.outputCollection
    if output is None:
        output = nirSuperdark.scratchCollectionName(args.output, args.ticket, args.tag)
    processes = args.processes
    if processes is None:
        processes = len(args.spectrograph or NIR_SPECTROGRAPHS)

    command = buildCommand(args.repo, inputs, output, args.visits,
                           spectrographs=args.spectrograph, processes=processes,
                           irpFilter=args.irpFilter, useIRP=args.useIRP,
                           logLevel=args.logLevel, config=args.config, show=args.show)
    if args.dryRun:
        # Only the command goes to stdout, so that it can be pasted or piped;
        # shlex.join, since the data-id expression contains spaces and quotes.
        print(f"output collection: {output}", file=sys.stderr)
        print(shlex.join(command))
        return
    print(f"output collection: {output}", flush=True)
    status = subprocess.call(command)
    if status == 0 and not args.show:
        print(f"\nrawISRCubes written to {output}; combine them with:\n\n    "
              + shlex.join(combineCommand(args.repo, output, args.visits, args.spectrograph,
                                          args.output, args.ticket, args.tag)))
    sys.exit(status)


if __name__ == "__main__":
    main()
