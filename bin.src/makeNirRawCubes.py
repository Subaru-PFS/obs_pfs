#!/usr/bin/env python
"""Generate the raw NIR ramp cubes (``rawISRCube``) that combineNirDark consumes.

Runs the ``isr`` step of ``reduceExposure.yaml`` over the given NIR visits with
the corrections that a dark must *not* have applied: no dark subtraction, no
flat, no linearity, and no IRP smoothing. The point of a nirDark is to capture
the full raw instrument dark signature, which is subtracted from an exposure ramp
before any corrections, so the ramps it is built from must themselves be raw.

The overridden configuration, and the default it replaces:

===============================  ========  =====
Config                           Default   Here
===============================  ========  =====
``isr.doDark``                   True      False
``isr.doFlat``                   True      False
``isr.h4.doWriteRawCube``        False     True
``isr.h4.doLinearize``           True      False
``isr.h4.quickCDS``              False     False
===============================  ========  =====

Reference-pixel handling is selectable: ``--irp-filter`` sets ``isr.h4.IRPfilter``
(default 0 = use IRP with no smoothing; -1 = per-column median; odd 15..31 =
Hann-smoothed), and ``--no-irp`` sets ``isr.h4.useIRP=False`` to bypass the IRP
planes entirely and border-correct instead. The dark must be built the same way
as the exposures it will be subtracted from.

Like ``combineNirDark.py``, ``--output`` is a collection base beneath which the
DMTN-222 layout is composed: the cubes go to a ``scratchCubes`` output collection
at ``{output}/{ticket}/{tag}/scratchCubes``. They are an intermediate product, not
a calibration, so they get no ``Gen.{iteration}`` name and are not certified.
Pipetask makes that a CHAINED collection and writes a timestamped RUN inside it,
so re-running does not collide; pass the chain straight to ``combineNirDark
--input``.

Example, reproducing the IRP1 nirDark inputs::

    makeNirRawCubes.py /work/datastore \\
        --input u/cpl/calib/PIPE2D-1858/PIPE2D-1857/badRefPixelsGen.20000101a \\
        --output u/cpl/calib --ticket PIPE2D-1664 --tag irp1 \\
        --visits 144587..144636 -j 12

writes to ``u/cpl/calib/PIPE2D-1664/irp1/scratchCubes``.

``--show config`` (or ``--show uri``) is passed through to pipetask, and reports
the resolved configuration without running anything.
"""

from __future__ import annotations

import shlex
import subprocess
import sys
from argparse import ArgumentParser

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
)

# The default IRP smoothing: use IRP planes with no smoothing, matching the
# nirDark recipe. `IRPfilter=0` = no smoothing; -1 = per-column median; odd
# 15..31 = Hann-smoothed. `useIRP=False` bypasses IRP for border correction.
DEFAULT_IRP_FILTER = 0

DEFAULT_INPUTS = ("PFS/defaults",)


def irpConfig(irpFilter: int, useIRP: bool) -> tuple[str, ...]:
    """The reference-pixel config overrides.

    ``useIRP`` is always stated explicitly; ``IRPfilter`` is emitted only when IRP
    is in use, since it does nothing when the reference planes are bypassed.
    """
    overrides = [f"isr:h4.useIRP={useIRP}"]
    if useIRP:
        overrides.append(f"isr:h4.IRPfilter={irpFilter}")
    return tuple(overrides)


def scratchCollection(outputBase: str, ticket: str, tag: str) -> str:
    """The output collection for the raw cubes, within the DMTN-222 layout.

    The cubes are an intermediate product rather than a calibration, so they get
    no ``{product}Gen.{iteration}`` name. Pipetask makes this a CHAINED collection
    and writes a timestamped RUN into it, so repeated runs do not collide.
    """
    return f"{outputBase}/{ticket}/{tag}/scratchCubes"


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


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("repo", help="Path to the butler repository")
    parser.add_argument("--input", nargs="+", required=True, dest="inputs",
                        help="Input collection(s); PFS/defaults is appended if absent")
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
                        help="Spectrograph number(s) (default: all)")
    parser.add_argument("-j", "--processes", type=int, default=1,
                        help="Number of pipetask worker processes (default: 1)")
    parser.add_argument("--irp-filter", type=int, default=DEFAULT_IRP_FILTER, dest="irpFilter",
                        help="IRP smoothing window: 0=none (default), -1=per-column median, "
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
        output = scratchCollection(args.output, args.ticket, args.tag)

    command = buildCommand(args.repo, inputs, output, args.visits,
                           spectrographs=args.spectrograph, processes=args.processes,
                           irpFilter=args.irpFilter, useIRP=args.useIRP,
                           logLevel=args.logLevel, config=args.config, show=args.show)
    if args.dryRun:
        # shlex.join, so the printed command can be pasted into a shell: the
        # data-id expression contains spaces and quotes.
        print(shlex.join(command))
        return
    sys.exit(subprocess.call(command))


if __name__ == "__main__":
    main()
