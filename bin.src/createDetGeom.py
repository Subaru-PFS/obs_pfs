#!/usr/bin/env python
# example:
# createDetGeom.py
#
# Note that CreateDetGeomTask is a Task and not a CmdLineTask,
# so it does not have the standard CmdLineTask arguments.The reason
# for this is that CmdLineTasks expects arguments like dataId which
# are of no use here.
#
from lsst.obs.pfs.createDetGeom import CreateDetGeomTask
task = CreateDetGeomTask()
task.run()
