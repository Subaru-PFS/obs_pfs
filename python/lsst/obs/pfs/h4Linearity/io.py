"""FITS persistence for LinearityCorrection objects."""

from __future__ import annotations

import datetime as _dt
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
from astropy.io import fits

from .models import MODEL_REGISTRY
from .types import Diagnostics, LinearityCorrection


def _packageVersion() -> str:
    try:
        return version("nirLinearity")
    except PackageNotFoundError:
        return "unknown"


def saveFits(path: str | Path, correction: LinearityCorrection) -> None:
    """Write a ``LinearityCorrection`` to a FITS file."""
    path = Path(path)

    # Build PRIMARY header.
    primaryHeader = fits.Header()
    primaryHeader["MODEL"] = (
        correction.model.modelName, "model form identifier"
    )
    primaryHeader["FITDATE"] = (
        _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "ISO-8601 fit timestamp",
    )
    primaryHeader["RELINVER"] = (_packageVersion(), "nirLinearity package version")
    # Scalar summary fields. Long keys become HIERARCH cards (case-preserved);
    # short keys are uppercased by FITS and the original Python key is stored in
    # the comment so ``loadFits`` can reconstruct the dict without collisions.
    for key, value in correction.diagnostics.summary.items():
        if isinstance(value, (int, float, bool, str)):
            if len(key) > 8:
                # HIERARCH card: keyword IS the Python key; no comment needed.
                primaryHeader["HIERARCH " + key] = (value, "")
            else:
                # Short key: FITS uppercases it; store original in comment.
                primaryHeader[key] = (value, key)
    primary = fits.PrimaryHDU(header=primaryHeader)

    # Model-specific HDUs.
    modelHdus = list(correction.model.toFitsHdus(correction))

    # Standard HDUs for the non-model-specific arrays.
    fitMinHdu = fits.ImageHDU(data=correction.fitMin, name="FITMIN")
    fitMaxHdu = fits.ImageHDU(data=correction.fitMax, name="FITMAX")
    bpHdu = fits.ImageHDU(data=correction.badPixelMask, name="BPMASK")
    bpHdu.header["COMMENT"] = "Bit flags: MASKED_BY_INPUT=0x01 INSUFFICIENT_POINTS=0x02"
    bpHdu.header["COMMENT"] = "          FIT_FAILED=0x04 NON_MONOTONIC=0x08 BORDER_PIX=0x10"
    bpHdu.header["COMMENT"] = "          BELOW_VALID_RANGE=0x20 ABOVE_VALID_RANGE=0x40"
    resRmsHdu = fits.ImageHDU(
        data=correction.diagnostics.residualRms, name="RESRMS"
    )
    resMaxHdu = fits.ImageHDU(
        data=correction.diagnostics.maxAbsResidual, name="RESMAX"
    )
    nPtsHdu = fits.ImageHDU(
        data=correction.diagnostics.nPointsUsed, name="NPOINTS"
    )
    monoHdu = fits.ImageHDU(
        data=correction.diagnostics.monotonic.astype(np.uint8), name="MONOTON"
    )
    condHdu = fits.ImageHDU(
        data=correction.diagnostics.conditionNumber, name="CONDNUM"
    )

    hdul = fits.HDUList(
        [primary, *modelHdus, fitMinHdu, fitMaxHdu, bpHdu,
         resRmsHdu, resMaxHdu, nPtsHdu, monoHdu, condHdu]
    )

    # Add CHECKSUM/DATASUM to every image HDU.
    for hdu in hdul[1:]:
        hdu.add_checksum()

    hdul.writeto(path, overwrite=True)


def loadFits(path: str | Path) -> LinearityCorrection:
    """Read a FITS file written by :func:`saveFits`."""
    path = Path(path)
    with fits.open(path) as hdul:
        primary = hdul[0]
        modelName = primary.header["MODEL"]
        if modelName not in MODEL_REGISTRY:
            raise ValueError(
                f"Unknown model {modelName!r}; known: {sorted(MODEL_REGISTRY)}"
            )
        modelClass = MODEL_REGISTRY[modelName]

        # Collect model HDUs (anything the model classmethod consumes) and the
        # fixed non-model HDUs by name.
        allHdus = [hdu for hdu in hdul]
        model, coefficients = modelClass.fromFitsHdus(allHdus)

        fitMin = _arrayByName(hdul, "FITMIN")
        fitMax = _arrayByName(hdul, "FITMAX")
        badPixelMask = _arrayByName(hdul, "BPMASK").astype(np.uint8)
        residualRms = _arrayByName(hdul, "RESRMS")
        maxAbsResidual = _arrayByName(hdul, "RESMAX")
        nPointsUsed = _arrayByName(hdul, "NPOINTS").astype(np.int32)
        monotonic = _arrayByName(hdul, "MONOTON").astype(bool)
        conditionNumber = _arrayByName(hdul, "CONDNUM")

        # Rebuild summary from primary header (best-effort; drops non-scalar keys).
        # HIERARCH cards preserve the Python key in card.keyword (mixed case /
        # underscores); for short keys that FITS uppercased, the original Python
        # key is stored in the comment.
        summary: dict = {}
        _skipKeys = {"SIMPLE", "BITPIX", "NAXIS", "EXTEND", "MODEL",
                     "FITDATE", "RELINVER"}
        for card in primary.header.cards:
            key = card.keyword
            if key in _skipKeys or key.startswith("NAXIS"):
                continue
            if key == key.upper():
                # Standard FITS key (uppercased); recover original via comment.
                if card.comment:
                    originalKey = card.comment
                    summary[originalKey] = card.value
            else:
                # HIERARCH card — keyword IS the original Python key.
                summary[key] = card.value

    diagnostics = Diagnostics(
        residualRms=residualRms,
        maxAbsResidual=maxAbsResidual,
        nPointsUsed=nPointsUsed,
        monotonic=monotonic,
        conditionNumber=conditionNumber,
        summary=summary,
    )
    return LinearityCorrection(
        model=model,
        coefficients=coefficients,
        fitMin=fitMin,
        fitMax=fitMax,
        badPixelMask=badPixelMask,
        diagnostics=diagnostics,
    )


def _arrayByName(hdul: fits.HDUList, name: str) -> np.ndarray:
    for hdu in hdul:
        if getattr(hdu, "name", "") == name:
            return np.asarray(hdu.data)
    raise ValueError(f"HDU {name!r} not found in FITS file")
