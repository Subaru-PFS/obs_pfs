"""Tests for the ``makeNirRawCubes`` command-line wrapper.

The wrapper assembles a ``pipetask run`` command that produces the ``rawISRCube``
ramps ``combineNirDark`` consumes. The tests check the assembled command line, the
data-id expression, and -- most importantly -- that every config override names a
field that really exists on `PfsIsrTask.ConfigClass`, since ``pipetask`` would
otherwise reject them only at run time, hours into a submission.
"""

import sys
import unittest

import lsst.utils.tests

from testUtils import HAS_DRP_STELLA, loadScript, requireDrpStella

if HAS_DRP_STELLA:
    # isrTask imports pfs.drp.stella.crosstalk, and the pipeline this script runs
    # lives in drp_stella, so pipelinePath() needs the product set up too.
    from lsst.obs.pfs.isrTask import PfsIsrTask


@requireDrpStella
class BuildCommandTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.script = loadScript("makeNirRawCubes")

    def testVisitQueryPassesRangesThrough(self):
        # The butler understands BEGIN..END itself, so no expansion is needed.
        self.assertEqual(self.script.visitQuery(["144587..144636"]),
                         "visit in (144587..144636) and arm='n'")
        self.assertEqual(self.script.visitQuery(["1", "3..5"]),
                         "visit in (1, 3..5) and arm='n'")

    def testVisitQuerySpectrographs(self):
        self.assertEqual(
            self.script.visitQuery(["7"], [1, 3]),
            "visit in (7) and arm='n' and spectrograph in (1, 3)")

    def testCommandLine(self):
        command = self.script.buildCommand(
            "/work/datastore", ["u/me/badRefPixels", "PFS/defaults"], "u/me/out",
            ["144587..144636"], processes=12)
        self.assertEqual(command[0], "pipetask")
        self.assertIn("run", command)
        self.assertIn("--fail-fast", command)
        for flag, value in (("-b", "/work/datastore"),
                            ("-o", "u/me/out"),
                            ("-i", "u/me/badRefPixels,PFS/defaults"),
                            ("-j", "12"),
                            ("-d", "visit in (144587..144636) and arm='n'")):
            self.assertEqual(command[command.index(flag) + 1], value)
        self.assertTrue(command[command.index("-p") + 1].endswith(
            "/pipelines/reduceExposure.yaml#isr"))

    def overrides(self, command):
        return [command[i + 1] for i, a in enumerate(command) if a == "-c"]

    def testRawCubeConfigApplied(self):
        command = self.script.buildCommand("/repo", ["in"], "out", ["7"])
        # Fixed corrections, then the default IRP config (use IRP, no smoothing).
        self.assertEqual(self.overrides(command),
                         list(self.script.ISR_CONFIG)
                         + ["isr:h4.useIRP=True", "isr:h4.IRPfilter=0"])

    def testIrpFilterSelectable(self):
        command = self.script.buildCommand("/repo", ["in"], "out", ["7"], irpFilter=15)
        self.assertIn("isr:h4.IRPfilter=15", self.overrides(command))
        self.assertNotIn("isr:h4.IRPfilter=0", self.overrides(command))

    def testPerColumnMedianMode(self):
        command = self.script.buildCommand("/repo", ["in"], "out", ["7"], irpFilter=-1)
        self.assertIn("isr:h4.IRPfilter=-1", self.overrides(command))

    def testNoIrpBypassesReferencePlanes(self):
        command = self.script.buildCommand("/repo", ["in"], "out", ["7"], useIRP=False)
        overrides = self.overrides(command)
        self.assertIn("isr:h4.useIRP=False", overrides)
        # With IRP bypassed, no IRPfilter is emitted -- it does nothing.
        self.assertFalse([o for o in overrides if o.startswith("isr:h4.IRPfilter")])

    def testExtraConfigAppliedLast(self):
        command = self.script.buildCommand("/repo", ["in"], "out", ["7"],
                                           config=["isr:h4.doCR=True"])
        self.assertEqual(self.overrides(command)[-1], "isr:h4.doCR=True")

    def testLogLevelOmittedByDefault(self):
        self.assertNotIn("--log-level", self.script.buildCommand("/repo", ["in"], "out", ["7"]))
        command = self.script.buildCommand("/repo", ["in"], "out", ["7"], logLevel=".=DEBUG")
        self.assertEqual(command[command.index("--log-level") + 1], ".=DEBUG")

    def testShowPassedThrough(self):
        self.assertNotIn("--show", self.script.buildCommand("/repo", ["in"], "out", ["7"]))
        command = self.script.buildCommand("/repo", ["in"], "out", ["7"], show=["config", "uri"])
        shown = [command[i + 1] for i, a in enumerate(command) if a == "--show"]
        self.assertEqual(shown, ["config", "uri"])

    def testScratchCollectionLayout(self):
        self.assertEqual(self.script.scratchCollection("u/cpl/calib", "PIPE2D-1664", "irp4"),
                         "u/cpl/calib/PIPE2D-1664/irp4/scratchCubes")


@requireDrpStella
class MainTestCase(lsst.utils.tests.TestCase):
    """``main`` marshals its arguments the way combineNirDark does."""

    def setUp(self):
        self.script = loadScript("makeNirRawCubes")
        self.origArgv = sys.argv

    def tearDown(self):
        sys.argv = self.origArgv

    def run_main(self, argv):
        sys.argv = ["makeNirRawCubes"] + argv
        self.script.main()

    def capture(self, argv):
        """Run main() in --dry-run and return the printed command as a list."""
        import io
        import contextlib
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.run_main(argv + ["--dry-run"])
        import shlex
        return shlex.split(out.getvalue().strip())

    def testDefaultOutputIsScratchCubes(self):
        command = self.capture(["/repo", "--input", "in", "--output", "u/me/calib",
                                "--ticket", "PIPE2D-1664", "--tag", "irp4",
                                "--visits", "7"])
        self.assertEqual(command[command.index("-o") + 1],
                         "u/me/calib/PIPE2D-1664/irp4/scratchCubes")

    def testOutputCollectionOverride(self):
        command = self.capture(["/repo", "--input", "in", "--output", "u/me/calib",
                                "--ticket", "T", "--tag", "g", "--visits", "7",
                                "--output-collection", "u/me/explicit"])
        self.assertEqual(command[command.index("-o") + 1], "u/me/explicit")

    def testDefaultsAppendedToInputs(self):
        command = self.capture(["/repo", "--input", "in", "--output", "u/me/calib",
                                "--ticket", "T", "--tag", "g", "--visits", "7"])
        self.assertEqual(command[command.index("-i") + 1], "in,PFS/defaults")

    def testDefaultsNotDuplicated(self):
        command = self.capture(["/repo", "--input", "in", "PFS/defaults",
                                "--output", "u/me/calib", "--ticket", "T", "--tag", "g",
                                "--visits", "7"])
        self.assertEqual(command[command.index("-i") + 1], "in,PFS/defaults")

    def testSpectrographRestriction(self):
        command = self.capture(["/repo", "--input", "in", "--output", "u/me/calib",
                                "--ticket", "T", "--tag", "g", "--visits", "7",
                                "--spectrograph", "1", "2"])
        self.assertEqual(command[command.index("-d") + 1],
                         "visit in (7) and arm='n' and spectrograph in (1, 2)")

    def testTicketAndTagRequired(self):
        for missing in (["--ticket", "T"], ["--tag", "g"]):
            with self.assertRaises(SystemExit):
                self.run_main(["/repo", "--input", "in", "--output", "u/me/calib",
                               "--visits", "7", *missing, "--dry-run"])

    def testIrpFlagsReachCommand(self):
        command = self.capture(["/repo", "--input", "in", "--output", "u/me/calib",
                                "--ticket", "T", "--tag", "g", "--visits", "7",
                                "--irp-filter", "-1"])
        overrides = [command[i + 1] for i, a in enumerate(command) if a == "-c"]
        self.assertIn("isr:h4.IRPfilter=-1", overrides)

    def testNoIrpFlagReachesCommand(self):
        command = self.capture(["/repo", "--input", "in", "--output", "u/me/calib",
                                "--ticket", "T", "--tag", "g", "--visits", "7", "--no-irp"])
        overrides = [command[i + 1] for i, a in enumerate(command) if a == "-c"]
        self.assertIn("isr:h4.useIRP=False", overrides)
        self.assertFalse([o for o in overrides if o.startswith("isr:h4.IRPfilter")])


@requireDrpStella
class IsrConfigFieldsTestCase(lsst.utils.tests.TestCase):
    """Every override must name a real config field, and actually change it."""

    def setUp(self):
        self.script = loadScript("makeNirRawCubes")
        self.config = PfsIsrTask.ConfigClass()

    def assertFieldExists(self, override):
        label, _, assignment = override.partition(":")
        self.assertEqual(label, "isr")
        name, _, _value = assignment.partition("=")
        target = self.config
        *parents, leaf = name.split(".")
        for parent in parents:
            target = getattr(target, parent)
        self.assertTrue(hasattr(target, leaf), f"{name} is not a PfsIsrTask config field")

    def testOverriddenFieldsExist(self):
        # Every override -- fixed, and both IRP branches -- must name a real field,
        # or pipetask would reject it only at run time, hours into a submission.
        allOverrides = (self.script.ISR_CONFIG
                        + self.script.irpConfig(0, True)
                        + self.script.irpConfig(-1, True)
                        + self.script.irpConfig(0, False))
        for override in allOverrides:
            self.assertFieldExists(override)

    def testRawCubesRequireUncorrectedRamps(self):
        """The darks must be combined without linearity or a dark subtracted."""
        overrides = dict(
            override.split(":", 1)[1].split("=", 1) for override in self.script.ISR_CONFIG)
        self.assertEqual(overrides["doDark"], "False")
        self.assertEqual(overrides["h4.doLinearize"], "False")
        self.assertEqual(overrides["h4.doWriteRawCube"], "True")

    def testDefaultIrpFilterNoSmoothing(self):
        overrides = dict(o.split(":", 1)[1].split("=", 1)
                         for o in self.script.irpConfig(self.script.DEFAULT_IRP_FILTER, True))
        self.assertEqual(overrides["h4.useIRP"], "True")
        self.assertEqual(overrides["h4.IRPfilter"], "0")

    def testOverridesDifferFromDefaults(self):
        """A no-op override would mean the default drifted under us."""
        self.assertTrue(self.config.doDark)
        self.assertTrue(self.config.h4.doLinearize)
        self.assertFalse(self.config.h4.doWriteRawCube)
        self.assertNotEqual(self.config.h4.IRPfilter, self.script.DEFAULT_IRP_FILTER)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
