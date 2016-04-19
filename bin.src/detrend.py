#!/usr/bin/env python
#/Users/azuri/stella-git/obs_pfs/bin.src/detrend.py '/Users/azuri/spectra/pfs/PFS' --calib='/Users/azuri/spectra/pfs/PFS/CALIB' --id site='S' category='A' filter='PFS-M' spectrograph=2 dateObs='2015-12-22' -C myConfig.py

from lsst.obs.pfs.detrendTask import DetrendTask
DetrendTask.parseAndRun()

