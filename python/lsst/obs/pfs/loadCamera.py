import os
from functools import lru_cache

from lsst.afw.cameraGeom import Camera
from lsst.obs.base.yamlCamera import makeCamera
from lsst.utils import getPackageDir

__all__ = ("loadCamera",)


@lru_cache
def loadCamera(category: str = None) -> Camera:
    """Retrieve the cameraGeom representation of an instrument

    Parameters
    ----------
    category : `str`, optional
        Category of PFS. This should be one of ``S`` (Subaru), ``L`` (LAM) or
        ``F`` (fake/simulator), or ``None``.

    Returns
    -------
    camera : `Camera`
        Description of instrument.
    """
    if category is not None:
        yamlFile = os.path.join(
            getPackageDir("obs_pfs"), "pfs-" + category, "camera", "camera.yaml"
        )
        if os.path.exists(yamlFile):
            return makeCamera(yamlFile)
    yamlFile = os.path.join(getPackageDir("obs_pfs"), "pfs", "camera", "camera.yaml")
    return makeCamera(yamlFile)
