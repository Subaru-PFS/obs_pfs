# -*- python -*-
import os
from lsst.sconsUtils import scripts
scripts.BasicSConstruct(
    "obs_pfs",
    defaultTargets=scripts.DEFAULT_TARGETS + ("corrections",),
    disableCc=True,
    subDirList=([path for path in os.listdir(".") if os.path.isdir(path) and not path.startswith(".")] +
                ["bin"]),
)
