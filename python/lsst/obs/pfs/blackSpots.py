import os
from dataclasses import dataclass
from typing import Optional, overload

import astropy.io.ascii
import numpy as np
from lsst.pex.config import Config, Field
from lsst.pipe.base import Struct
from lsst.utils import getPackageDir
from numpy.typing import ArrayLike
from scipy.spatial import KDTree


@dataclass
class BlackSpots:
    """Black spot positions

    Parameters
    ----------
    spotId : `numpy.ndarray`
        Spot ID.
    x : `numpy.ndarray`
        x-coordinate of spot (mm).
    y : `numpy.ndarray`
        y-coordinate of spot (mm).
    r : `numpy.ndarray`
        Radius of spot (mm).
    tree : `scipy.spatial.KDTree`
        KDTree for finding nearest black spot.
    """

    spotId: np.ndarray
    x: np.ndarray
    y: np.ndarray
    r: np.ndarray
    tree: KDTree

    @classmethod
    def read(cls, filename: str, leafsize: int = 10) -> "BlackSpots":
        """Read black spot positions from file

        Parameters
        ----------
        filename : `str`
            Fully-qualified path for file containing black spot positions.

        Returns
        -------
        blackSpots : `BlackSpot`
            Black spot positions.
        """
        blackSpots = astropy.io.ascii.read(filename, format="csv")
        tree = KDTree(np.array([blackSpots["x"], blackSpots["y"]]).T, leafsize=leafsize)
        return cls(
            blackSpots["spotId"],
            blackSpots["x"],
            blackSpots["y"],
            blackSpots["r"],
            tree,
        )

    @overload
    def findNearest(self, x: ArrayLike, y: ArrayLike) -> Struct:
        """Find black spots nearest to a point

        Parameters
        ----------
        x : array_like
            x-coordinates of point (mm).
        y : array_like
            y-coordinates of point (mm).

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``distance``: distance to nearest black spot (mm).
            - ``spotId``: spotId of nearest black spot; -1 if no spot was found.
            - ``x``: x-coordinate of nearest black spot (mm).
            - ``y``: y-coordinate of nearest black spot (mm).
            - ``r``: radius of nearest black spot (mm).
        """
        ...

    @overload
    def findNearest(self, x: ArrayLike) -> Struct:
        """Find black spots nearest to a point

        Parameters
        ----------
        xy : array_like
            x- and y- coordinates of point (mm).

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``distance``: distance to nearest black spot (mm).
            - ``spotId``: spotId of nearest black spot; -1 if no spot was found.
            - ``x``: x-coordinate of nearest black spot (mm).
            - ``y``: y-coordinate of nearest black spot (mm).
            - ``r``: radius of nearest black spot (mm).
        """
        ...

    def findNearest(self, x: ArrayLike, y: Optional[ArrayLike] = None) -> Struct:
        distance, index = self.tree.query(x if y is None else np.array([x, y]).T, k=1)

        # Deal with bad values
        if np.isscalar(x):
            if index == self.tree.n:
                return Struct(
                    distance=np.nan, spotId=-1, x=np.nan, y=np.nan, r=np.nan
                )
            else:
                return Struct(
                    distance=distance,
                    spotId=self.spotId[index],
                    x=self.x[index],
                    y=self.y[index],
                    r=self.r[index],
                )
        good = index != self.tree.n
        if np.all(good):
            return Struct(
                distance=distance,
                spotId=self.spotId[index],
                x=self.x[index],
                y=self.y[index],
                r=self.r[index],
            )
        spotId = np.full(index.shape, -1, dtype=int)
        xx = np.full(index.shape, np.nan, dtype=float)
        yy = np.full(index.shape, np.nan, dtype=float)
        rr = np.full(index.shape, np.nan, dtype=float)
        goodIndex = index[good]
        spotId[good] = self.spotId[goodIndex]
        xx[good] = self.x[goodIndex]
        yy[good] = self.y[goodIndex]
        rr[good] = self.r[goodIndex]
        return Struct(distance=distance, spotId=spotId, x=xx, y=yy, r=rr)


class BlackSpotsConfig(Config):
    """Configuration for black spot positions"""

    filename = Field(
        dtype=str,
        default="pfs/black_dots_mm-2022-09-23.csv",
        doc="File containing black spot positions",
    )
    leafsize = Field(dtype=int, default=10, doc="Leaf size for KDTree")

    def read(self):
        """Read black spot positions from file

        Returns
        -------
        blackSpots : `BlackSpots`
            Black spot positions.
        """
        filename = self.filename
        if not os.path.isabs(filename):
            filename = os.path.join(getPackageDir("obs_pfs"), filename)
        return BlackSpots.read(filename, leafsize=self.leafsize)
