"""Tests for the ``combineNirDark`` command-line wrapper.

The wrapper is a thin driver over ``lsst.obs.pfs.nirSuperdark.processMasterDark``
that expands visit ranges and reduces the requested spectrographs in parallel.
These tests confirm argument marshalling, visit-range parsing, the parallel
fan-out, and per-spectrograph error aggregation.
"""

import datetime
import os
import shutil
import sys
import tempfile
import unittest

import lsst.utils.tests

from testUtils import HAS_DRP_STELLA, loadScript, requireDrpStella

if HAS_DRP_STELLA:
    # nirSuperdark imports pfs.drp.stella.calibs.setCalibHeader.
    from lsst.obs.pfs import nirSuperdark


def markerWorker(inputRun, outputRun, dataId, visits, startDate=None, repo_path=None,
                 saveDir=None, calibCollection=None):
    """Stand-in for ``processMasterDark`` that records a run as a file.

    Module-level (hence picklable) so it can be dispatched to worker processes
    in the parallel-execution test.
    """
    with open(os.path.join(repo_path, f"n{dataId['spectrograph']}.done"), "w") as fd:
        fd.write(repr((inputRun, outputRun, sorted(visits), startDate)))


VISIT_START = datetime.datetime(2026, 6, 25, 5, 23, 17)
TIMESTAMP = "20260709T123456Z"


def fakePreflight(inputRun, outputBase, dataId, visits, ticket, tag, iteration,
                  timestamp, repo_path=None):
    """Stand-in for ``preflight``: no butler to name collections against here."""
    gen = nirSuperdark.genCollectionName(outputBase, ticket, tag, "nirDark", iteration)
    return nirSuperdark.Plan(
        datasetType="nirDark",
        calibCollection=nirSuperdark.calibCollectionName(outputBase, ticket, tag,
                                                         "nirDark", iteration),
        genCollection=gen,
        outputRun=nirSuperdark.runName(gen, timestamp),
        startDate=VISIT_START,
    )


@requireDrpStella
class CombineNirDarkTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.script = loadScript("combineNirDark")
        self.calls = []
        self.origArgv = sys.argv
        self.origProcess = nirSuperdark.processMasterDark
        self.origPreflight = nirSuperdark.preflight

        def recorder(*args, **kwargs):
            self.calls.append((args, kwargs))

        nirSuperdark.processMasterDark = recorder
        nirSuperdark.preflight = fakePreflight

    def tearDown(self):
        nirSuperdark.processMasterDark = self.origProcess
        nirSuperdark.preflight = self.origPreflight
        sys.argv = self.origArgv

    def run_main(self, argv):
        # Force inline execution so the monkeypatched recorder (a closure, not
        # picklable) runs in this process and records the calls.
        extra = ["--processes", "1"]
        if "--ticket" not in argv:
            extra += ["--ticket", "PIPE2D-1664", "--tag", "tag",
                      "--iteration", "20260709a"]
        if "--run-timestamp" not in argv:
            extra += ["--run-timestamp", TIMESTAMP]
        sys.argv = ["combineNirDark"] + argv + extra
        self.script.main()

    def testArgsMarshalled(self):
        self.run_main([
            "/some/repo",
            "--input", "u/me/in",
            "--output", "u/me/out",
            "--spectrograph", "1",
            "--visits", "10", "20", "30",
        ])
        self.assertEqual(len(self.calls), 1)
        args, kwargs = self.calls[0]
        inputRun, outputRun, dataId, visits = args
        self.assertEqual(inputRun, "u/me/in")
        # Names composed under the --output base: the certified collection drops
        # the Gen suffix, and the datasets land in a timestamped run beneath it.
        base = "u/me/out/PIPE2D-1664/tag"
        self.assertEqual(outputRun, f"{base}/nirDarkGen.20260709a/{TIMESTAMP}")
        self.assertEqual(kwargs["calibCollection"], f"{base}/nirDark.20260709a")
        self.assertEqual(dataId, dict(instrument="PFS", arm="n", spectrograph=1))
        self.assertEqual(visits, [10, 20, 30])
        self.assertEqual(kwargs["repo_path"], "/some/repo")
        # Resolved from the first visit by the parent, not left for the worker.
        self.assertEqual(kwargs["startDate"], VISIT_START)

    def testAllSpectrographsShareOneStartDate(self):
        """One dark set gets one validity start, not one per detector's timestamp."""
        self.run_main(["/repo", "--input", "in", "--output", "out", "--visits", "7"])
        self.assertEqual(len(self.calls), 4)
        startDates = {kwargs["startDate"] for _, kwargs in self.calls}
        self.assertEqual(startDates, {VISIT_START})

    def testStartDateParsedToDatetime(self):
        """The start date reaches saveNirDark as a datetime, not a string."""
        self.run_main([
            "/repo", "--input", "in", "--output", "out", "--spectrograph", "1",
            "--visits", "7", "--start-date", "2026-01-02T03:04:05",
        ])
        self.assertEqual(self.calls[0][1]["startDate"],
                         datetime.datetime(2026, 1, 2, 3, 4, 5))

    def testIterationDefaultsToToday(self):
        self.run_main([
            "/repo", "--input", "in", "--output", "out", "--spectrograph", "1",
            "--visits", "7", "--ticket", "PIPE2D-1664", "--tag", "tag",
        ])
        expected = nirSuperdark.defaultIteration()
        self.assertTrue(self.calls[0][1]["calibCollection"].endswith(f"nirDark.{expected}"))

    def testRunTimestampDefaultsToNow(self):
        self.run_main([
            "/repo", "--input", "in", "--output", "out", "--spectrograph", "1",
            "--visits", "7", "--run-timestamp", "20991231T235959Z",
        ])
        self.assertTrue(self.calls[0][0][1].endswith("/20991231T235959Z"))

    def testAllSpectrographsShareOneRun(self):
        """A rerun gets a fresh RUN, but one run holds the whole set."""
        self.run_main(["/repo", "--input", "in", "--output", "out", "--visits", "7"])
        outputRuns = {args[1] for args, _ in self.calls}
        self.assertEqual(len(outputRuns), 1)
        self.assertTrue(outputRuns.pop().endswith(f"/{TIMESTAMP}"))

    def testDefaultsToAllSpectrographs(self):
        self.run_main([
            "/repo", "--input", "in", "--output", "out", "--visits", "7",
        ])
        self.assertEqual(len(self.calls), 4)
        spectrographs = [args[2]["spectrograph"] for args, _ in self.calls]
        self.assertEqual(spectrographs, [1, 2, 3, 4])
        for args, _ in self.calls:
            _, _, dataId, visits = args
            self.assertEqual(dataId["arm"], "n")
            self.assertEqual(visits, [7])

    def testStartDateAndSpectrographSubset(self):
        self.run_main([
            "/repo",
            "--input", "in", "--output", "out",
            "--spectrograph", "3", "4",
            "--visits", "7",
            "--start-date", "2026-01-02T03:04:05",
        ])
        self.assertEqual([args[2]["spectrograph"] for args, _ in self.calls], [3, 4])
        for _, kwargs in self.calls:
            self.assertEqual(kwargs["startDate"], datetime.datetime(2026, 1, 2, 3, 4, 5))

    def testSaveDirDefaultsToCwd(self):
        self.run_main(["/repo", "--input", "in", "--output", "out",
                       "--spectrograph", "1", "--visits", "7"])
        self.assertEqual(self.calls[0][1]["saveDir"], ".")

    def testSaveDirCreated(self):
        root = tempfile.mkdtemp(prefix="combineNirDark-save-")
        saveDir = os.path.join(root, "nested", "out")
        try:
            self.run_main(["/repo", "--input", "in", "--output", "out",
                           "--spectrograph", "1", "--visits", "7",
                           "--save-dir", saveDir])
            self.assertTrue(os.path.isdir(saveDir))
            self.assertEqual(self.calls[0][1]["saveDir"], saveDir)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def testPreflightFailureAbortsBeforeAnyWork(self):
        """A bad output collection or dataset type must not cost a combine."""
        def boom(*args, **kwargs):
            raise ValueError("output collection 'out' is of type CHAINED")

        nirSuperdark.preflight = boom
        with self.assertRaises(ValueError):
            self.script.reduceSpectrographs("/repo", "in", "out", [1, 2, 3, 4], [7],
                                            "PIPE2D-1664", "tag", "20260709a", TIMESTAMP,
                                            processes=1)
        self.assertEqual(self.calls, [])

    def testVisitsRequired(self):
        with self.assertRaises(SystemExit):
            self.run_main(["/repo", "--input", "in", "--output", "out",
                           "--spectrograph", "1"])

    def testVisitRangesExpanded(self):
        self.run_main([
            "/repo", "--input", "in", "--output", "out", "--spectrograph", "1",
            "--visits", "144784..144788", "200", "300..306:2",
        ])
        _, _, _, visits = self.calls[0][0]
        self.assertEqual(
            visits,
            [144784, 144785, 144786, 144787, 144788, 200, 300, 302, 304, 306])

    def testErrorsAggregated(self):
        """One spectrograph failing does not abandon the others."""
        def failing(inputRun, outputRun, dataId, visits, **kwargs):
            self.calls.append(dataId["spectrograph"])
            if dataId["spectrograph"] == 2:
                raise ValueError("boom")

        nirSuperdark.processMasterDark = failing
        with self.assertRaises(RuntimeError) as cm:
            self.script.reduceSpectrographs(
                "/repo", "in", "out", [1, 2, 3], [7],
                "PIPE2D-1664", "tag", "20260709a", TIMESTAMP, processes=1)
        self.assertIn("n2", str(cm.exception))
        self.assertEqual(sorted(self.calls), [1, 2, 3])

    def testParallelExecutionRunsAll(self):
        """The parallel path dispatches every spectrograph to a worker process."""
        nirSuperdark.processMasterDark = markerWorker
        root = tempfile.mkdtemp(prefix="combineNirDark-par-")
        try:
            self.script.reduceSpectrographs(
                root, "in", "out", [1, 2, 3, 4], [7],
                "PIPE2D-1664", "tag", "20260709a", TIMESTAMP, processes=2)
            self.assertEqual(sorted(os.listdir(root)),
                             ["n1.done", "n2.done", "n3.done", "n4.done"])
        finally:
            shutil.rmtree(root, ignore_errors=True)


@requireDrpStella
class ParseVisitsTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.script = loadScript("combineNirDark")

    def testPlainInts(self):
        self.assertEqual(self.script.parseVisits(["10", "20"]), [10, 20])

    def testInclusiveRange(self):
        self.assertEqual(self.script.parseVisitRange("5..8"), [5, 6, 7, 8])

    def testRangeWithStep(self):
        self.assertEqual(self.script.parseVisitRange("1..10:3"), [1, 4, 7, 10])

    def testMixed(self):
        self.assertEqual(
            self.script.parseVisits(["1", "3..5", "9"]), [1, 3, 4, 5, 9])

    def testBadRangesRaise(self):
        for bad in ("8..5", "1..10:0", "1..2..3", "1..x"):
            with self.assertRaises(ValueError):
                self.script.parseVisitRange(bad)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules["__main__"])
    unittest.main(failfast=True)
