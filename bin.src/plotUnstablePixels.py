#!/usr/bin/env python

import os
from argparse import ArgumentParser

from astropy import log
from lsst.daf.butler import Butler
from lsst.obs.pfs.h4Linearity.rateStabilityDiagnostic import (
    runRateStabilityDiagnostic,
    makeUnstablePixelReport,
)


def main(repo, inputs, dataQuery, outputDir, *,
         threshold, rateFloorADU, minDeltasPerSegment,
         faintMaxRate, brightMinRate,
         nFaintPages, nBrightPages,
         panelsPerPage, ncols):
    """Write one rate-stability PDF report per matched NIR exposure.

    The report's pixel pages are split into a "fainter" band
    (``rateFloorADU <= |rate| < faintMaxRate``) and a "brighter" band
    (``|rate| >= brightMinRate``), with up to ``nFaintPages`` /
    ``nBrightPages`` of ``panelsPerPage`` pixels each.

    Parameters
    ----------
    repo : str
        Butler repository root.
    inputs : list of str
        Input collections.
    dataQuery : str
        pipetask-style data-query string.
    outputDir : str
        Directory for the output PDFs (created if absent).
    """
    os.makedirs(outputDir, exist_ok=True)
    butler = Butler(repo, collections=list(inputs))
    dataIds = sorted(set(butler.registry.queryDataIds(
        ["visit", "arm", "spectrograph", "detector"],
        where=dataQuery, datasets="raw", collections=list(inputs))))
    if not dataIds:
        log.warning(f"No raw exposures matched query: {dataQuery!r}")
        return
    for dataId in dataIds:
        if dataId["arm"] != "n":
            log.warning(f"Skipping non-NIR exposure {dataId}")
            continue
        data = runRateStabilityDiagnostic(
            butler, dataId,
            threshold=threshold, rateFloorADU=rateFloorADU,
            minDeltasPerSegment=minDeltasPerSegment)
        outPath = os.path.join(
            outputDir, f"unstable-{data.cam}-{data.visit}.pdf")
        makeUnstablePixelReport(
            data, outPath,
            faintMaxRate=faintMaxRate,
            brightMinRate=brightMinRate,
            nFaintPages=nFaintPages,
            nBrightPages=nBrightPages,
            panelsPerPage=panelsPerPage,
            ncols=ncols,
        )
        log.info(f"{data.cam} visit={data.visit}: "
                 f"{data.result.nRejected} RATE_UNSTABLE pixels -> {outPath}")


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Plot rate-stability rejects to per-exposure PDFs")
    parser.add_argument("repo", help="Butler repository root")
    parser.add_argument("-i", "--input", required=True,
                        help="Input collection(s), comma-separated")
    parser.add_argument("-d", "--data-query", required=True,
                        help="Data query, e.g. \"instrument='PFS' AND "
                             "visit IN (1,2) AND arm='n' AND spectrograph=3\"")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Directory for the output PDF files")
    parser.add_argument("--threshold", type=float, default=0.20,
                        help="Rate-stability fractional rejection threshold "
                             "(default 0.20)")
    parser.add_argument("--rate-floor-adu", type=float, default=5.0,
                        help="Denominator floor for the fractional metric "
                             "(default 5.0)")
    parser.add_argument("--min-deltas-per-segment", type=int, default=3,
                        help="Min un-flagged deltas for a testable segment")
    parser.add_argument("--faint-max-rate", type=float, default=20.0,
                        help="Upper edge of the fainter band (ADU/read, "
                             "default 20.0)")
    parser.add_argument("--bright-min-rate", type=float, default=20.0,
                        help="Lower edge of the brighter band (ADU/read, "
                             "default 20.0)")
    parser.add_argument("--n-faint-pages", type=int, default=4,
                        help="Number of fainter-band pages (default 4)")
    parser.add_argument("--n-bright-pages", type=int, default=4,
                        help="Number of brighter-band pages (default 4)")
    parser.add_argument("--panels-per-page", type=int, default=32,
                        help="Pixel cells per page (default 32)")
    parser.add_argument("--ncols", type=int, default=4,
                        help="Pixel cells per row (default 4)")
    args = parser.parse_args()
    main(
        args.repo, args.input.split(","), args.data_query, args.output_dir,
        threshold=args.threshold,
        rateFloorADU=args.rate_floor_adu,
        minDeltasPerSegment=args.min_deltas_per_segment,
        faintMaxRate=args.faint_max_rate,
        brightMinRate=args.bright_min_rate,
        nFaintPages=args.n_faint_pages,
        nBrightPages=args.n_bright_pages,
        panelsPerPage=args.panels_per_page,
        ncols=args.ncols,
    )
