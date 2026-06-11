"""Tests for fit() threading: heuristic, parallel path, errors, determinism."""

from __future__ import annotations

import sys
import threading

import numpy as np
import pytest

from lsst.obs.pfs.h4Linearity.fit import _resolveWorkerCount, fit
from lsst.obs.pfs.h4Linearity.models import PolynomialModel


def test_resolveWorkerCountExplicitIntIsReturnedAsIs():
    # Explicit wins over heuristic — no clamping, no size check.
    assert _resolveWorkerCount(1, 10, 10) == 1
    assert _resolveWorkerCount(4, 10, 10) == 4
    assert _resolveWorkerCount(16, 10, 10) == 16
    # Even on a "large" frame, explicit 1 is honored.
    assert _resolveWorkerCount(1, 5000, 5000) == 1


def test_resolveWorkerCountSmallFrameDefaultsToOne():
    # With H*W < _SMALL_FRAME_PIXEL_LIMIT (1_000_000), None → 1 worker
    # regardless of os.cpu_count().
    assert _resolveWorkerCount(None, 100, 100) == 1
    assert _resolveWorkerCount(None, 1000, 999) == 1  # 999_000 < 1_000_000


def test_resolveWorkerCountLargeFrameCapsAtEight(monkeypatch):
    # When H*W >= 1_000_000 and os.cpu_count() > 8, cap at 8.
    fitModule = sys.modules["lsst.obs.pfs.h4Linearity.fit"]
    monkeypatch.setattr(fitModule.os, "cpu_count", lambda: 16)
    assert _resolveWorkerCount(None, 2000, 500) == 8  # 1_000_000 exactly
    assert _resolveWorkerCount(None, 4096, 4096) == 8


def test_resolveWorkerCountLargeFrameUncappedBelowEight(monkeypatch):
    # When H*W is large but os.cpu_count() < 8, use cpu_count.
    fitModule = sys.modules["lsst.obs.pfs.h4Linearity.fit"]
    monkeypatch.setattr(fitModule.os, "cpu_count", lambda: 4)
    assert _resolveWorkerCount(None, 4096, 4096) == 4


def test_resolveWorkerCountHandlesNoneCpuCount(monkeypatch):
    # os.cpu_count() can return None on some platforms; fall back to 1.
    fitModule = sys.modules["lsst.obs.pfs.h4Linearity.fit"]
    monkeypatch.setattr(fitModule.os, "cpu_count", lambda: None)
    assert _resolveWorkerCount(None, 4096, 4096) == 1


def test_resolveWorkerCountInvalidRaises():
    with pytest.raises(ValueError, match="workers"):
        _resolveWorkerCount(0, 10, 10)
    with pytest.raises(ValueError, match="workers"):
        _resolveWorkerCount(-3, 10, 10)


def _arraysEqual(a, b, label):
    """np.array_equal with a helpful assertion message."""
    assert a.shape == b.shape, f"{label}: shape {a.shape} != {b.shape}"
    assert a.dtype == b.dtype, f"{label}: dtype {a.dtype} != {b.dtype}"
    assert np.array_equal(a, b), f"{label}: arrays differ"


def test_fitWorkers4ProducesByteIdenticalOutputToWorkers1(smallSyntheticRamp):
    """The threaded path must produce output byte-identical to the sequential
    path, because each tile writes to disjoint output slices on the main
    thread and the tile-assembly and fit arithmetic are deterministic."""
    ramp, _ = smallSyntheticRamp
    # blockSize=(2, 3) on a 4x5 frame produces 4 tiles — enough to exercise
    # parallelism when workers=4.
    serial = fit([ramp], blockSize=(2, 3), workers=1)
    threaded = fit([ramp], blockSize=(2, 3), workers=4)

    _arraysEqual(threaded.coefficients, serial.coefficients, "coefficients")
    _arraysEqual(threaded.fitMin, serial.fitMin, "fitMin")
    _arraysEqual(threaded.fitMax, serial.fitMax, "fitMax")
    _arraysEqual(threaded.badPixelMask, serial.badPixelMask, "badPixelMask")
    _arraysEqual(
        threaded.diagnostics.residualRms,
        serial.diagnostics.residualRms,
        "residualRms",
    )
    _arraysEqual(
        threaded.diagnostics.maxAbsResidual,
        serial.diagnostics.maxAbsResidual,
        "maxAbsResidual",
    )
    _arraysEqual(
        threaded.diagnostics.nPointsUsed,
        serial.diagnostics.nPointsUsed,
        "nPointsUsed",
    )
    _arraysEqual(
        threaded.diagnostics.monotonic,
        serial.diagnostics.monotonic,
        "monotonic",
    )
    _arraysEqual(
        threaded.diagnostics.conditionNumber,
        serial.diagnostics.conditionNumber,
        "conditionNumber",
    )
    # Summary dicts must be equal (same float values, same keys).
    assert threaded.diagnostics.summary == serial.diagnostics.summary


def test_fitWorkers1DoesNotConstructExecutor(monkeypatch, smallSyntheticRamp):
    """The sequential fast path must not touch the executor factory at all."""
    fitModule = sys.modules["lsst.obs.pfs.h4Linearity.fit"]

    def _fail(*args, **kwargs):
        raise AssertionError(
            f"_executorFactory was called for workers=1 path: "
            f"args={args} kwargs={kwargs}"
        )

    monkeypatch.setattr(fitModule, "_executorFactory", _fail)
    ramp, _ = smallSyntheticRamp
    fit([ramp], blockSize=(2, 3), workers=1)  # must not raise


def test_fitWorkers4ConstructsExecutorWithMaxWorkers4(
    monkeypatch, smallSyntheticRamp
):
    """When workers=4, the factory must be called with max_workers=4."""
    fitModule = sys.modules["lsst.obs.pfs.h4Linearity.fit"]
    from concurrent.futures import ThreadPoolExecutor

    recorded = {}

    def _recordingFactory(*args, **kwargs):
        recorded["args"] = args
        recorded["kwargs"] = dict(kwargs)
        # Return the real executor so the fit still runs to completion.
        return ThreadPoolExecutor(*args, **kwargs)

    monkeypatch.setattr(fitModule, "_executorFactory", _recordingFactory)
    ramp, _ = smallSyntheticRamp
    fit([ramp], blockSize=(2, 3), workers=4)
    assert recorded["kwargs"].get("max_workers") == 4


def test_fitAutoWorkersUsesResolvedCount(monkeypatch, smallSyntheticRamp):
    """With workers=None and a small frame, resolved count is 1 (no executor
    call). This mirrors test_fitWorkers1DoesNotConstructExecutor but via the
    default code path rather than explicit workers=1."""
    fitModule = sys.modules["lsst.obs.pfs.h4Linearity.fit"]

    called = []
    monkeypatch.setattr(
        fitModule, "_executorFactory", lambda *a, **k: called.append(1)
    )
    ramp, _ = smallSyntheticRamp  # H=4 W=5, well below 1_000_000.
    fit([ramp], blockSize=(2, 3))  # workers=None
    assert called == []


def test_fitInvalidWorkersRaises(smallSyntheticRamp):
    ramp, _ = smallSyntheticRamp
    with pytest.raises(ValueError, match="workers"):
        fit([ramp], workers=0)
    with pytest.raises(ValueError, match="workers"):
        fit([ramp], workers=-1)


def test_fitWorkerExceptionIncludesTileCoords(smallSyntheticRamp):
    """If any fitBlock call raises on a worker thread, the exception must
    be re-raised as a RuntimeError whose message identifies the offending
    tile's row/col slice and whose __cause__ is the original exception."""
    ramp, _ = smallSyntheticRamp
    pm = PolynomialModel(order=2)
    originalFitBlock = pm.fitBlock
    failLock = threading.Lock()
    failedOnce = threading.Event()

    def failingFitBlock(m, t, valid, conditionNumberLimit):
        # Atomically check-and-set: without the lock, two concurrent workers
        # could both see is_set() == False and both raise. The lock guarantees
        # exactly one thread takes the failure branch.
        with failLock:
            if not failedOnce.is_set():
                failedOnce.set()
                raise RuntimeError("injected failure")
        return originalFitBlock(
            m=m, t=t, valid=valid, conditionNumberLimit=conditionNumberLimit
        )

    # PolynomialModel is a frozen dataclass; bypass the frozen __setattr__
    # to shadow the bound method with an instance attribute.
    object.__setattr__(pm, "fitBlock", failingFitBlock)

    with pytest.raises(
        RuntimeError,
        match=r"fitBlock failed at tile \[rows \d+:\d+, cols \d+:\d+\]",
    ) as excInfo:
        fit([ramp], model=pm, blockSize=(2, 3), workers=2)
    # __cause__ carries the original exception.
    cause = excInfo.value.__cause__
    assert isinstance(cause, RuntimeError)
    assert str(cause) == "injected failure"
