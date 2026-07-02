#!/usr/bin/env python
"""Create a Gen3 butler repository for PFS, ready for use.

The repository is created with the PFS dimension universe and storage classes,
and the instrument (together with its calibration dataset types, e.g.
``nirLinearity``, ``badRefPixels`` and ``defects``) is registered. The result
can be used directly with ``makePfsCalibs.py`` to install calibration products.
"""

from __future__ import annotations

import os
from argparse import ArgumentParser

from lsst.daf.butler import Butler, DimensionConfig
from lsst.utils import getPackageDir

from lsst.obs.pfs.PrimeFocusSpectrograph import (
    PfsDevelopment,
    PfsSimulator,
    PrimeFocusSpectrograph,
)

INSTRUMENTS = {
    "PFS": PrimeFocusSpectrograph,
    "PFS-F": PfsSimulator,
    "PFS-L": PfsDevelopment,
}


def makePfsTestRepo(repo: str, instrument: str = "PFS") -> Butler:
    """Create and populate a PFS butler repository.

    Parameters
    ----------
    repo : `str`
        Path at which to create the repository.
    instrument : `str`, optional
        Instrument to register (``PFS``, ``PFS-F`` for the simulator, or
        ``PFS-L`` for LAM development data).

    Returns
    -------
    butler : `lsst.daf.butler.Butler`
        A writeable butler for the newly created repository.
    """
    obsPfsDir = getPackageDir("obs_pfs")
    dimensionConfig = DimensionConfig(os.path.join(obsPfsDir, "gen3", "dimensions.yaml"))
    Butler.makeRepo(
        repo,
        config=os.path.join(obsPfsDir, "gen3", "butler.yaml"),
        dimensionConfig=dimensionConfig,
    )
    butler = Butler(repo, writeable=True)
    INSTRUMENTS[instrument]().register(butler.registry)
    return butler


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("repo", help="Path to the butler repository to create")
    parser.add_argument(
        "--instrument",
        default="PFS",
        choices=sorted(INSTRUMENTS),
        help="Instrument to register (default: PFS)",
    )
    args = parser.parse_args()
    makePfsTestRepo(args.repo, args.instrument)
    print(f"Created {args.instrument} butler repo at {args.repo}")


if __name__ == "__main__":
    main()
