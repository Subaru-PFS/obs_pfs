import re

import numpy as np
import pandas as pd

from pfs.drp.stella.referenceLine import (
    ReferenceLineSet,
    ReferenceLineSource,
    ReferenceLineStatus,
)

__all__ = ("readSweepFile",)

# Threshold on ``sigma_wav`` (nm) below which a clean, robust line is
# considered "super-high confidence" for wavelength calibration.
SIGMA_WAV_THRESHOLD = 0.001

_SPECIES_RE = re.compile(r"^[A-Za-z0-9]+")


def _speciesOf(component: str) -> str:
    """Identify the species of a single ``label`` component.

    Parameters
    ----------
    component : `str`
        A single (non-``|``-separated) contributor from the RAGNAR
        ``label`` column, e.g. ``"(9-0)pP2e(11.5)"`` or ``"NaI(D1)"``.

    Returns
    -------
    species : `str`
        Species identification, e.g. ``"OH"``, ``"NaI"``, ``"O2"``, ``"OI"``.
    """
    if component.startswith("("):
        return "OH"
    match = _SPECIES_RE.match(component)
    return match.group(0) if match else component


def _descriptionOf(label: str) -> str:
    """Derive the ``description`` field from a RAGNAR ``label`` string.

    Parameters
    ----------
    label : `str`
        ``|``-separated list of contributing transitions, from the RAGNAR
        ``label`` column.

    Returns
    -------
    description : `str`
        Sorted, deduplicated, ``|``-joined set of species names
        contributing to this line (e.g. ``"OH"``, ``"O2|OH"``).
    """
    species = sorted(set(_speciesOf(component) for component in label.split("|")))
    return "|".join(species)


def _statusOf(clean: bool, robust: bool, sigmaWav: float) -> ReferenceLineStatus:
    """Determine the ``status`` of a line from RAGNAR's quality columns.

    A line is considered good for "super-high confidence" wavelength
    calibration only if it is ``clean``, ``robust``, and has
    ``sigma_wav`` below `SIGMA_WAV_THRESHOLD`; otherwise it is flagged
    as a blend.

    Parameters
    ----------
    clean : `bool`
        Do all contributing transitions share the same upper level?
    robust : `bool`
        Was the peak found in all simulated spectra?
    sigmaWav : `float`
        Standard deviation of the peak position (nm) across simulated
        spectra.

    Returns
    -------
    status : `ReferenceLineStatus`
        ``GOOD`` or ``BLEND``.
    """
    if clean and robust and sigmaWav < SIGMA_WAV_THRESHOLD:
        return ReferenceLineStatus.GOOD
    return ReferenceLineStatus.BLEND


def readSweepFile(filename: str) -> ReferenceLineSet:
    """Read a RAGNAR sweep file and convert it to a `ReferenceLineSet`.

    Parameters
    ----------
    filename : `str`
        Path to a RAGNAR sweep file (``.pkl``), containing a
        ``pandas.DataFrame`` with (at least) the columns ``wav``,
        ``intensity``, ``label``, ``clean``, ``robust`` and ``sigma_wav``.

    Returns
    -------
    lines : `ReferenceLineSet`
        Reference lines converted from the sweep file, in the same
        (wavelength-sorted) order as the input.
    """
    data: pd.DataFrame = pd.read_pickle(filename)

    description = np.array([_descriptionOf(label) for label in data["label"]])
    status = np.array(
        [
            _statusOf(clean, robust, sigmaWav)
            for clean, robust, sigmaWav in zip(data["clean"], data["robust"], data["sigma_wav"])
        ]
    )

    return ReferenceLineSet.fromColumns(
        wavelength=data["wav"].to_numpy(),
        intensity=data["intensity"].to_numpy(),
        description=description,
        status=status,
        transition=data["label"].to_numpy(),
        source=np.full(len(data), ReferenceLineSource.RAGNAR),
    )
