try:
    from .version import *
except ImportError:
    print("WARNING: unable to import version.py in obs_pfs; did you build with scons?")
    __version__ = "unknown"

from .PrimeFocusSpectrograph import *
