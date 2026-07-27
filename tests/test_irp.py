import types
import unittest

import numpy as np
import lsst.utils.tests

from lsst.obs.pfs.isrTask import PfsIsrTask
from lsst.obs.pfs.h4utils import irp4
from lsst.obs.pfs.nirBadRefPixels import NirBadRefPixels


def makePfsRawStub(irpN, irpOffset, nchan=4, badPixels=None, readOrder=(0, 0)):
    """Build a minimal stand-in for `PfsRaw` for the IRP geometry tests.

    Only the attributes touched by the IRP interpolation/bad-pixel methods are
    provided.
    """
    detector = types.SimpleNamespace(getName=lambda: "stub")
    return types.SimpleNamespace(
        irpN=irpN,
        irpOffset=irpOffset,
        nchan=nchan,
        h4Size=4096,
        detector=detector,
        getH4channelReadOrder=lambda: readOrder,
        _badPixels=np.array([] if badPixels is None else badPixels, dtype=int),
    )


class IrpTestCase(lsst.utils.tests.TestCase):
    """Geometry of IRPn -> IRP1 interpolation and bad-pixel reinterpretation."""

    def setUp(self):
        self.task = PfsIsrTask(config=PfsIsrTask.ConfigClass())
        self.task.config.h4.doIRPbadPixels = True
        self.task.config.h4.IRPfilter = 0

    def testIrpPhaseOffsets(self):
        """W_H4IRPO -> forward phase irpOffset-1, reversed phase irpN-irpOffset."""
        # n1: irpN=4, W_H4IRPO=4 -> forward 3, reversed 0
        self.assertEqual(self.task.irpPhaseOffsets(makePfsRawStub(irpN=4, irpOffset=4)), (3, 0))
        # a mid-group reference: irpOffset=2 -> forward 1, reversed 2
        self.assertEqual(self.task.irpPhaseOffsets(makePfsRawStub(irpN=4, irpOffset=2)), (1, 2))
        # IRP1: both phases 0
        self.assertEqual(self.task.irpPhaseOffsets(makePfsRawStub(irpN=1, irpOffset=1)), (0, 0))

    def testBadRefLocalSamples(self):
        """Bad pixels are restricted to a channel and mapped to local samples by phase."""
        # channel 30 spans IRP1 rows [3840, 3968); 3899 = 3840 + 56 + 3 (offset 3)
        bad = np.array([3899])
        self.assertFloatsEqual(
            self.task.badRefLocalSamples(bad, 30 * 128, 128, 4, 3), np.array([14]))
        # at the reversed phase (0) that pixel is not sampled
        self.assertEqual(len(self.task.badRefLocalSamples(bad, 30 * 128, 128, 4, 0)), 0)
        # a pixel in another channel is ignored
        self.assertEqual(len(self.task.badRefLocalSamples(np.array([100]), 30 * 128, 128, 4, 3)), 0)
        # IRP1: the in-channel row maps straight through
        self.assertFloatsEqual(
            self.task.badRefLocalSamples(np.array([3850]), 30 * 128, 128, 1, 0), np.array([10]))

    def testReplaceBadRefSamplesBeforeElseAfter(self):
        """Bad samples take the previous neighbour, or the next at a channel edge."""
        chan = np.arange(8, dtype="f4").reshape(8, 1) * 10.0  # 0,10,20,...,70
        # sample 0 (channel start -> use after=1), sample 4 (-> use before=3)
        self.task.replaceBadRefSamples(chan, np.array([0, 4]))
        expected = np.array([10, 10, 20, 30, 30, 50, 60, 70], dtype="f4").reshape(8, 1)
        self.assertFloatsEqual(chan, expected)

    def testReplaceRunPropagatesForward(self):
        """A run of bad samples is filled from the nearest good sample below it."""
        chan = np.arange(6, dtype="f4").reshape(6, 1) * 10.0  # 0,10,...,50
        # samples 2,3 bad: both take sample 1 (nearest good below), skipping the run
        self.task.replaceBadRefSamples(chan, np.array([2, 3]))
        expected = np.array([0, 10, 10, 10, 40, 50], dtype="f4").reshape(6, 1)
        self.assertFloatsEqual(chan, expected)

    def testReplaceBadRefSamplesEdges(self):
        """Runs at the channel start/end, and an all-bad channel, are handled."""
        # adjacent bad at the start [0,1] -> both take the first good sample (2)
        chan = np.arange(4, dtype="f4").reshape(4, 1) * 10.0  # 0,10,20,30
        self.task.replaceBadRefSamples(chan, np.array([0, 1]))
        self.assertFloatsEqual(chan, np.array([20, 20, 20, 30], dtype="f4").reshape(4, 1))
        # bad at the very end -> takes the good sample below
        chan = np.arange(4, dtype="f4").reshape(4, 1) * 10.0
        self.task.replaceBadRefSamples(chan, np.array([3]))
        self.assertFloatsEqual(chan, np.array([0, 10, 20, 20], dtype="f4").reshape(4, 1))
        # a run spanning a good pixel: each bad takes its nearest good below
        chan = np.arange(5, dtype="f4").reshape(5, 1) * 10.0  # 0,10,20,30,40
        self.task.replaceBadRefSamples(chan, np.array([1, 3]))
        self.assertFloatsEqual(chan, np.array([0, 0, 20, 20, 40], dtype="f4").reshape(5, 1))
        # every sample bad -> nothing to copy from, left untouched
        chan = np.arange(3, dtype="f4").reshape(3, 1) * 10.0
        self.task.replaceBadRefSamples(chan, np.array([0, 1, 2]))
        self.assertFloatsEqual(chan, np.array([0, 10, 20], dtype="f4").reshape(3, 1))

    def testConstructFullIrp4Expands(self):
        """IRP4 image expands to full 4096 by repeating each sample irpN times."""
        nchan = 4
        irpN = 4
        height = 4096 // irpN  # 1024
        width = 8
        # distinct value per subsampled row so we can check the expansion pattern
        refImg = (np.arange(height, dtype="f4")[:, None] * np.ones(width)).copy()
        pfsRaw = makePfsRawStub(irpN=irpN, irpOffset=1, nchan=nchan)
        self.task._badRefPixels = NirBadRefPixels.fromList([], "stub")

        full = self.task.constructFullIrp(pfsRaw, refImg)
        self.assertEqual(full.shape, (4096, width))
        # each source row r maps to output rows [r*irpN, r*irpN+irpN)
        for r in (0, 1, 100, height - 1):
            block = full[r * irpN:(r + 1) * irpN, 0]
            self.assertFloatsEqual(block, np.full(irpN, float(r)))

    def testConstructFullIrp4BadSampleReplaced(self):
        """A bad IRP4 sample is replaced before expansion, across its whole block."""
        nchan = 4
        irpN = 4
        height = 4096 // irpN
        width = 4
        refImg = (np.arange(height, dtype="f4")[:, None] * np.ones(width)).copy()
        # W_H4IRPO=1 -> forward phase 0; forward channel 0, pixel 20 -> sample 5.
        pfsRaw = makePfsRawStub(irpN=irpN, irpOffset=1, nchan=nchan, badPixels=[20])
        self.task._badRefPixels = NirBadRefPixels.fromList([20], "stub")

        full = self.task.constructFullIrp(pfsRaw, refImg)
        # sample 5 should have been replaced by sample 4, then repeated x4
        block = full[5 * irpN:6 * irpN, 0]
        self.assertFloatsEqual(block, np.full(irpN, 4.0))
        # neighbours untouched
        self.assertFloatsEqual(full[4 * irpN:5 * irpN, 0], np.full(irpN, 4.0))
        self.assertFloatsEqual(full[6 * irpN:7 * irpN, 0], np.full(irpN, 6.0))

    def testConstructFullIrp4Idempotent(self):
        """A full-frame IRP4 image (already expanded, e.g. re-read via fromArray in
        the UTR path) is returned as-is, not re-expanded x irpN, even with a
        badRefPixels calib present."""
        pfsRaw = makePfsRawStub(irpN=4, irpOffset=4, nchan=32)
        self.task._badRefPixels = NirBadRefPixels.fromList([1064], "stub")
        full = (np.arange(4096, dtype="f4")[:, None] * np.ones(4)).copy()  # already 4096
        out = self.task.constructFullIrp(pfsRaw, full)
        self.assertEqual(out.shape, full.shape)  # not (16384, 4)
        self.assertFloatsEqual(out, full)

    def testConstructFullIrp4ReversedChannelPhase(self):
        """A reversed (odd) channel uses the mirrored phase irpN-irpOffset."""
        nchan = 4
        irpN = 4
        height = 4096 // irpN
        width = 4
        refChanHeight = height // nchan  # 256
        refImg = (np.arange(height, dtype="f4")[:, None] * np.ones(width)).copy()
        # n1-like: W_H4IRPO=4 -> forward 3, reversed 0; readOrder makes odd channels reversed.
        # channel 1 (odd) spans IRP1 rows [1024, 2048); pick reversed-phase pixel
        # p = 1024 + 4*10 + 0 = 1064 -> local sample 10.
        pfsRaw = makePfsRawStub(irpN=irpN, irpOffset=4, nchan=nchan,
                                badPixels=[1064], readOrder=(0, 1))
        self.task._badRefPixels = NirBadRefPixels.fromList([1064], "stub")

        full = self.task.constructFullIrp(pfsRaw, refImg)
        # local sample 10 (value refChanHeight+10) replaced by sample 9 (value refChanHeight+9)
        expected = float(refChanHeight + 9)
        self.assertFloatsEqual(full[1064:1068, 0], np.full(irpN, expected))
        # had we (wrongly) used the forward phase 3, pixel 1064 would not be sampled at all
        self.assertNotEqual(expected, float(refChanHeight + 10))


class IrpAlignmentTestCase(lsst.utils.tests.TestCase):
    """Offset alignment from observed noisy IRP rows (h4utils.irp4)."""

    def testMakeQuickIsrTask(self):
        """The inspection-task factory applies the IRP-friendly config."""
        task = irp4.makeQuickIsrTask()
        self.assertTrue(task.config.h4.quickCDS)
        self.assertEqual(task.config.h4.IRPfilter, 0)
        self.assertFalse(task.config.doDark)
        self.assertFalse(task.config.h4.doLinearize)
        # overrides land on the right (possibly nested) config attribute
        task2 = irp4.makeQuickIsrTask(**{"h4.IRPfilter": 15})
        self.assertEqual(task2.config.h4.IRPfilter, 15)

    def testIqrStd(self):
        """iqrStd recovers the sigma of a normal sample (IQR/1.349)."""
        rng = np.random.RandomState(1)
        a = rng.normal(0.0, 7.0, size=(4096,))
        self.assertFloatsAlmostEqual(irp4.iqrStd(a), 7.0, atol=0.5)

    def testCollapseToIrp1Pixels(self):
        """IRP4 noisy-row blocks collapse to the sampled IRP1 pixel per channel."""
        ro = (0, 1)  # even forward, odd reversed; n1 W_H4IRPO=4 -> fwd 3, rev 0
        # even channel (8) block -> phase 3 (block end)
        self.assertFloatsEqual(
            irp4.collapseToIrp1Pixels([1084, 1085, 1086, 1087], 4, 4, ro), np.array([1087]))
        # odd channel (31) block -> phase 0 (block start)
        self.assertFloatsEqual(
            irp4.collapseToIrp1Pixels([3980, 3981, 3982, 3983], 4, 4, ro), np.array([3980]))
        # several blocks across channels -> sorted, de-duplicated
        self.assertFloatsEqual(
            irp4.collapseToIrp1Pixels([2536, 2537, 2538, 2539, 3556, 3557, 3558, 3559],
                                      4, 4, ro),
            np.array([2536, 3556]))
        # IRP1: pass-through (unique)
        self.assertFloatsEqual(irp4.collapseToIrp1Pixels([5, 5, 7], 1, 1, ro), np.array([5, 7]))

    def testExcessRowNoise(self):
        """excessRowNoise subtracts each channel's own median baseline."""
        std = np.zeros(256, dtype="f4")
        std[:128] = 5.0   # channel 0 pedestal
        std[10] = 13.0
        std[128:] = 3.0   # channel 1 pedestal
        std[200] = 9.0
        excess = irp4.excessRowNoise(std, nchan=2)
        self.assertAlmostEqual(excess[10], 8.0, places=5)
        self.assertAlmostEqual(excess[200], 6.0, places=5)
        self.assertAlmostEqual(excess[0], 0.0, places=5)

    def testExcessRowNoiseSubtractsSlope(self):
        """A sloped channel is detrended, so its ends are not spuriously flagged."""
        x = np.arange(128)
        std = np.zeros(256, dtype="f4")
        std[:128] = 1.0 + 0.1 * x   # channel 0 ramps along its length
        std[128:] = 2.0             # channel 1 flat
        std[50] += 6.0              # spike sitting on the slope
        std[200] += 6.0             # spike on the flat channel
        excess = irp4.excessRowNoise(std, nchan=2)
        # on-slope non-spike rows are ~0 (would be ~6 at the end without detrending)
        self.assertLess(abs(excess[10]), 0.5)
        self.assertLess(abs(excess[120]), 0.5)
        # spikes survive at their full height
        self.assertAlmostEqual(excess[50], 6.0, delta=0.5)
        self.assertAlmostEqual(excess[200], 6.0, delta=0.5)

    def testCorrectedRowsForOffset(self):
        """A bad IRP1 pixel repairs its whole repeat-expanded group at one offset."""
        bad = np.array([3899])
        # 3899 = 974*4 + 3 -> only offset 3 maps it to a sample
        self.assertFloatsEqual(
            irp4.correctedRowsForOffset(bad, 4, 3),
            np.array([3896, 3897, 3898, 3899]),
        )
        for offset in (0, 1, 2):
            self.assertEqual(len(irp4.correctedRowsForOffset(bad, 4, offset)), 0)
        # IRP1: every bad pixel is its own row
        self.assertFloatsEqual(irp4.correctedRowsForOffset(bad, 1, 0), bad)

    def testCorrectedRowsPerChannelParity(self):
        """correctedRows uses forward phase on even channels, reversed on odd."""
        # n1-like: irpN=4, W_H4IRPO=4 -> forward 3, reversed 0; readOrders=(0,1)
        ro = (0, 1)
        # 3899 in channel 30 (even, forward phase 3) -> repairs its block
        self.assertFloatsEqual(irp4.correctedRows(np.array([3899]), 4, 4, ro),
                               np.array([3896, 3897, 3898, 3899]))
        # 3840 in channel 30 (even) is at the wrong phase -> not repaired
        self.assertEqual(len(irp4.correctedRows(np.array([3840]), 4, 4, ro)), 0)
        # a forward-channel pixel: 7 in channel 0 (even, phase 3) -> block [4,7]
        self.assertFloatsEqual(irp4.correctedRows(np.array([7]), 4, 4, ro),
                               np.array([4, 5, 6, 7]))
        # an odd channel uses the reversed phase 0: 200 in channel 1 -> block [200,203]
        self.assertFloatsEqual(irp4.correctedRows(np.array([200]), 4, 4, ro),
                               np.array([200, 201, 202, 203]))

    def testAlignmentPicksOffsetThatFixesBlock(self):
        """The n1/144800 case: rows 3896-3899 noisy, only 3899 known -> offset 3."""
        noisyRows = np.array([3896, 3897, 3898, 3899])
        badPixels = np.array([3899])
        coverage = irp4.offsetCoverage(noisyRows, badPixels, irpN=4)
        best = max(coverage, key=lambda o: coverage[o]["nCovered"])
        self.assertEqual(best, 3)
        self.assertEqual(coverage[3]["nCovered"], 4)
        self.assertEqual(coverage[3]["uncovered"], [])
        for offset in (0, 1, 2):
            self.assertEqual(coverage[offset]["nCovered"], 0)

    def testOffsetCoverageByParity(self):
        """Even/odd channels can prefer different offsets (readout-flip signature)."""
        # channel 0 (even): bad pixel at offset 3 -> noisy block 20..23
        # channel 1 (odd):  bad pixel at offset 0 -> noisy block 148..151
        noisyRows = np.array([20, 21, 22, 23, 148, 149, 150, 151])
        badPixels = np.array([23, 148])
        byParity = irp4.offsetCoverageByParity(noisyRows, badPixels, irpN=4, nchan=32)
        evenBest = max(byParity["even"], key=lambda o: byParity["even"][o]["nCovered"])
        oddBest = max(byParity["odd"], key=lambda o: byParity["odd"][o]["nCovered"])
        self.assertEqual(evenBest, 3)
        self.assertEqual(oddBest, 0)

    def testRampExcessAcrossReads(self):
        """rampExcess aggregates per-read excess and flags a persistently-noisy row."""
        scatter = np.array([-20, 20, -20, 20, -20, 20, -20, 20], dtype="f4")
        # 4 reads (256x8): row 10 only noisy between reads 0 and 1
        imgs = np.zeros((4, 256, 8), dtype="f4")
        imgs[1, 10, :] = scatter   # diff(0,1) and diff(1,2) make row 10 noisy

        class FakeTask:
            def __init__(self):
                self.config = types.SimpleNamespace(
                    h4=types.SimpleNamespace(doIRPbadPixels=True))

            def makeRawIrpArray(self, raw, rn):
                return raw._imgs[rn]

        raw = types.SimpleNamespace(nchan=2, irpN=1, _imgs=imgs,
                                    getNumReads=lambda: len(imgs),
                                    detector=types.SimpleNamespace(getName=lambda: "n1"))
        excess = irp4.rampExcess(FakeTask(), raw)
        self.assertGreater(excess[10], 10.0)
        self.assertEqual(np.where(excess > 2.0)[0].tolist(), [10])

    def testScanBadRefPixelsMinVisits(self):
        """A pixel must be flagged in >= minVisits darks to be proposed."""
        scatter = np.array([-20, 20, -20, 20, -20, 20, -20, 20], dtype="f4")

        def makeRaw(noisyRows):
            imgs = np.zeros((3, 256, 8), dtype="f4")
            for r in noisyRows:
                imgs[1, r, :] = scatter
            return types.SimpleNamespace(nchan=2, irpN=1, _imgs=imgs,
                                         getNumReads=lambda: 3,
                                         detector=types.SimpleNamespace(getName=lambda: "n1"))

        class FakeTask:
            def __init__(self):
                self.config = types.SimpleNamespace(
                    h4=types.SimpleNamespace(doIRPbadPixels=True))

            def makeRawIrpArray(self, raw, rn):
                return raw._imgs[rn]

        raws = [makeRaw([10, 40]), makeRaw([10, 50])]  # 10 in both, 40/50 in one each
        task = FakeTask()
        self.assertEqual(irp4.scanBadRefPixels(task, raws, minVisits=2).pixels.tolist(), [10])
        self.assertEqual(irp4.scanBadRefPixels(task, raws, minVisits=1).pixels.tolist(),
                         [10, 40, 50])

    def testScanCatchesIntermittentByMax(self):
        """The max pass flags a row bad in a single read that the mean dilutes."""
        step = np.array([-4, 4, -4, 4, -4, 4, -4, 4], dtype="f4")  # iqrStd ~5.9
        imgs = np.zeros((7, 256, 8), dtype="f4")
        imgs[3:, 30, :] = step                             # one noisy diff -> intermittent
        for k in range(7):
            imgs[k, 60, :] = step if k % 2 else -step      # noisy every diff -> persistent
        raw = types.SimpleNamespace(nchan=2, irpN=1, _imgs=imgs, getNumReads=lambda: 7,
                                    detector=types.SimpleNamespace(getName=lambda: "n1"))

        class FakeTask:
            def __init__(self):
                self.config = types.SimpleNamespace(
                    h4=types.SimpleNamespace(doIRPbadPixels=True))

            def makeRawIrpArray(self, raw, rn):
                return raw._imgs[rn]

        task = FakeTask()
        scan = irp4.scanBadRefPixels(task, [raw], threshold=2.0, maxThreshold=4.0)
        self.assertIn(30, scan.pixels.tolist())   # caught by max pass
        self.assertIn(60, scan.pixels.tolist())   # caught by mean
        # disabling the max pass (huge maxThreshold) drops the intermittent row
        scan2 = irp4.scanBadRefPixels(task, [raw], threshold=2.0, maxThreshold=1e9)
        self.assertNotIn(30, scan2.pixels.tolist())
        self.assertIn(60, scan2.pixels.tolist())

    def testRunBadRefPixelSurvey(self):
        """The survey runner writes a plot per detector and the proposed YAML."""
        import os
        import tempfile
        step = np.array([-4, 4, -4, 4, -4, 4, -4, 4], dtype="f4")
        imgs = np.zeros((5, 256, 8), dtype="f4")
        imgs[1, 20, :] = step

        class FakeTask:
            def __init__(self):
                self.config = types.SimpleNamespace(
                    h4=types.SimpleNamespace(doIRPbadPixels=True))

            def makeRawIrpArray(self, raw, rn):
                return raw._imgs[rn]

        raw = types.SimpleNamespace(nchan=2, irpN=1, _imgs=imgs, getNumReads=lambda: 5,
                                    detector=types.SimpleNamespace(getName=lambda: "n2"))

        class FakeButler:
            def get(self, datasetType, dataId=None):
                return raw

        with tempfile.TemporaryDirectory() as outdir:
            scans = irp4.runBadRefPixelSurvey(
                FakeButler(), [1, 2, 3], outdir=outdir, spectrographs=(2,),
                minVisits=1, task=FakeTask())
            self.assertEqual(len(scans), 1)
            self.assertTrue(os.path.exists(os.path.join(outdir, "scan_n2.png")))
            self.assertTrue(os.path.exists(os.path.join(outdir, "badRefPixels_proposed.yaml")))

    def testProposedBadRefPixelsYaml(self):
        """The YAML formatter emits a detector-keyed list of ints, plus optional metadata."""
        import yaml as _yaml
        scan = irp4.BadRefPixelScan("n2", 2.0, 4.0, 1, 1, {3: 1, 7: 1},
                                    np.array([3, 7]), np.zeros(4096), np.zeros(4096))
        loaded = _yaml.safe_load(irp4.proposedBadRefPixelsYaml([scan]))
        self.assertEqual(loaded, {"n2": [3, 7]})

        # with metadata: a top-level key the pixel-list loaders ignore
        meta = {"generatedBy": "tester", "visits": [144587], "threshold": 2.0}
        loaded = _yaml.safe_load(irp4.proposedBadRefPixelsYaml([scan], metadata=meta))
        self.assertEqual(loaded["n2"], [3, 7])
        self.assertEqual(loaded["metadata"]["generatedBy"], "tester")
        self.assertEqual(loaded["metadata"]["visits"], [144587])

    def testValidateCorrection(self):
        """validateCorrection separates repaired rows from candidate new bad pixels."""
        scatter = np.array([-20, 20, -20, 20, -20, 20, -20, 20], dtype="f4")  # iqrStd ~30

        class FakeTask:
            def __init__(self):
                self.config = types.SimpleNamespace(
                    h4=types.SimpleNamespace(doIRPbadPixels=True))

            def makeRawIrpNcube(self, raw, r0=0, r1=1):
                cube = np.zeros((2, 64, 8), dtype="f4")
                cube[1, 20, :] = scatter  # row 20 stays noisy even when corrected
                if not self.config.h4.doIRPbadPixels:
                    cube[1, 10, :] = scatter  # row 10 is repaired by the correction
                return cube

        raw = types.SimpleNamespace(nchan=4,
                                    detector=types.SimpleNamespace(getName=lambda: "n1"))
        task = FakeTask()
        result = irp4.validateCorrection(task, raw, threshold=2.0)
        self.assertFloatsEqual(np.sort(result.noisyRows), np.array([10, 20]))
        self.assertFloatsEqual(result.repairedRows, np.array([10]))
        self.assertFloatsEqual(result.remainingRows, np.array([20]))
        self.assertTrue(task.config.h4.doIRPbadPixels)  # original value restored

    def testAlignmentGuardsIrp1(self):
        """alignByCorrection refuses IRP1 data, where there is no offset to align."""
        raw = types.SimpleNamespace(
            irpN=1, irpOffset=0,
            detector=types.SimpleNamespace(getName=lambda: "n1"),
        )
        with self.assertRaises(ValueError):
            irp4.alignByCorrection(task=None, raw=raw)

    def testAlignmentReportsUncoveredAsCandidates(self):
        """Noisy rows with no matching known pixel surface as uncovered."""
        noisyRows = np.array([3896, 3897, 3898, 3899, 1000])  # 1000 not in list
        badPixels = np.array([3899])
        coverage = irp4.offsetCoverage(noisyRows, badPixels, irpN=4)
        self.assertEqual(coverage[3]["nCovered"], 4)
        self.assertEqual(coverage[3]["uncovered"], [1000])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
