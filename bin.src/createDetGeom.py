#!/usr/bin/env python
# example:
# createDetGeom.py 
#
# Note that CreateDetGeomTask is a Task and not a CmdLineTask,
# so it does not have the standard CmdLineTask arguments.
#
from lsst.obs.pfs.createDetGeom import CreateDetGeomTask
task = CreateDetGeomTask()
task.run()
