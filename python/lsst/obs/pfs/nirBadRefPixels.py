from __future__ import annotations

import numpy as np
import astropy.io.fits

from lsst.daf.base import PropertyList

from pfs.datamodel.utils import astropyHeaderFromDict

__all__ = ("NirBadRefPixels",)


def _applyProvenance(metadata: PropertyList, provenance: dict) -> None:
    """Record scan provenance in ``metadata`` as FITS-safe header cards.

    Maps the ``badRefPixels.yaml`` ``metadata`` block (who/when/visits/dates/
    thresholds) to <=8-char header keywords, stringifying list/dict values so
    they survive the FITS round-trip.
    """
    def add(key, value):
        if value is not None:
            metadata.set(key, value)

    add("CALGENBY", provenance.get("generatedBy"))
    add("CALGENAT", provenance.get("generatedAt"))
    add("CALARM", provenance.get("arm"))
    visits = provenance.get("visits")
    if visits:
        add("CALVIS", ",".join(str(v) for v in visits))
        add("CALNVIS", len(visits))
    dates = sorted(d for d in (provenance.get("visitDates") or {}).values() if d)
    if dates:
        add("CALDBEG", dates[0])
        add("CALDEND", dates[-1])
    add("CALTHR", provenance.get("threshold"))
    add("CALMAXTH", provenance.get("maxThreshold"))
    add("CALMINVI", provenance.get("minVisits"))
    nreads = provenance.get("nreads")
    add("CALNRD", "all" if nreads is None else nreads)


class NirBadRefPixels:
    """Bad IRP reference-row pixel list for a single NIR (H4RG) detector.

    Consists of a header-only primary HDU and a ``BADPIX`` binary table HDU
    with a single integer column ``PIXEL`` holding the row-pixel indices that
    are bad in the Interleaved Reference Pixel (IRP) data.

    Parameters
    ----------
    pixels : `np.ndarray`, 1-D int
        Row-pixel indices that are bad in the IRP.
    metadata : `lsst.daf.base.PropertyList`, optional
        Metadata (FITS header); carries e.g. the detector name (``DETNAME``).
    """

    def __init__(self, pixels: np.ndarray, metadata: "PropertyList | None" = None) -> None:
        self.pixels = np.asarray(pixels, dtype=np.int32).ravel()
        self.metadata = metadata if metadata is not None else PropertyList()

    @classmethod
    def fromList(cls, pixels, detectorName: "str | None" = None,
                 provenance: "dict | None" = None) -> "NirBadRefPixels":
        """Construct from a list of pixel indices

        Parameters
        ----------
        pixels : iterable of `int`
            Row-pixel indices that are bad in the IRP.
        detectorName : `str`, optional
            Detector name, recorded in the metadata as ``DETNAME``.
        provenance : `dict`, optional
            Scan provenance (the ``badRefPixels.yaml`` ``metadata`` block); its
            who/when/visits/dates/thresholds are recorded in the FITS header.

        Returns
        -------
        self : `NirBadRefPixels`
            The bad reference-pixel list.
        """
        metadata = PropertyList()
        if detectorName is not None:
            metadata.set("DETNAME", detectorName)
        if provenance:
            _applyProvenance(metadata, provenance)
        return cls(np.array(pixels, dtype=np.int32), metadata)

    @classmethod
    def fromFits(cls, fits: astropy.io.fits.HDUList) -> "NirBadRefPixels":
        """Construct from a FITS file

        Parameters
        ----------
        fits : `astropy.io.fits.HDUList`
            The FITS file.

        Returns
        -------
        self : `NirBadRefPixels`
            The bad reference-pixel list.
        """
        metadata = PropertyList.from_mapping(fits[0].header)
        hdu = fits["BADPIX"]
        if hdu.data is None:
            pixels = np.zeros(0, dtype=np.int32)
        else:
            pixels = np.asarray(hdu.data["PIXEL"], dtype=np.int32)
        return cls(pixels, metadata)

    def toFits(self) -> astropy.io.fits.HDUList:
        """Write to a FITS file

        Returns
        -------
        fits : `astropy.io.fits.HDUList`
            The FITS file.
        """
        fits = astropy.io.fits.HDUList()
        header = astropyHeaderFromDict(self.metadata)
        fits.append(astropy.io.fits.PrimaryHDU(data=None, header=header))
        column = astropy.io.fits.Column(name="PIXEL", format="J", array=self.pixels)
        fits.append(astropy.io.fits.BinTableHDU.from_columns([column], name="BADPIX"))
        return fits

    def writeFits(self, path: str) -> None:
        """Write the bad reference-pixel list to a FITS file

        Parameters
        ----------
        path : `str`
            Path to the output FITS file.
        """
        with open(path, "wb") as fd:
            self.toFits().writeto(fd)

    @classmethod
    def readFits(cls, path: str) -> "NirBadRefPixels":
        """Read the bad reference-pixel list from a FITS file

        Parameters
        ----------
        path : `str`
            Path to the FITS file.

        Returns
        -------
        self : `NirBadRefPixels`
            The bad reference-pixel list.
        """
        with astropy.io.fits.open(path) as fits:
            return cls.fromFits(fits)
