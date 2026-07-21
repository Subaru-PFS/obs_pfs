"""Shared helpers for the obs_pfs tests.

Some tests need packages that ``obs_pfs`` cannot declare as dependencies:

- ``drp_stella`` declares ``setupRequired(obs_pfs)``, so depending on it back
  would be circular. ``lsst.obs.pfs.isrTask`` and ``lsst.obs.pfs.nirSuperdark``
  nevertheless import from it, so tests touching them can only be skipped.
- ``drp_pfs_data`` is a bulky data package that a self-contained build does not
  set up.

Tests requiring either are skipped rather than failed, so ``scons`` passes on a
plain ``setup -r .`` checkout.
"""

from __future__ import annotations

import importlib.util
import os
import unittest

from lsst.utils import getPackageDir

OBS_PFS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def loadScript(name):
    """Import one of the ``bin.src`` scripts as a module."""
    path = os.path.join(OBS_PFS_DIR, "bin.src", f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hasPackage(name: str) -> bool:
    """Is the named EUPS product set up?"""
    try:
        getPackageDir(name)
    except LookupError:
        return False
    return True


def hasModule(name: str) -> bool:
    """Is the named Python module importable?"""
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


HAS_DRP_STELLA = hasModule("pfs.drp.stella")
HAS_DRP_PFS_DATA = hasPackage("drp_pfs_data")

requireDrpStella = unittest.skipUnless(
    HAS_DRP_STELLA,
    "drp_stella is not set up (it depends on obs_pfs, so obs_pfs cannot require it)")
requireDrpPfsData = unittest.skipUnless(
    HAS_DRP_PFS_DATA, "drp_pfs_data is not set up")
