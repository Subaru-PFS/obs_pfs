"""Chebyshev polynomial nonlinearity model: t = Σ c_k T_k(x) (per pixel)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.io import fits

from .base import BlockFitResult
from ..types import FIT_FAILED, INSUFFICIENT_POINTS, NON_MONOTONIC


@dataclass(frozen=True)
class PolynomialModel:
    """Pluggable Chebyshev polynomial-fit model. Default 4th order.

    ``fitMinMargin`` (default 100 DN) lowers the per-pixel ``fitMin`` by
    an absolute DN amount so that apply() does not flag BELOW_VALID_RANGE
    for read noise, kTC, or modest bias-level offsets that push slightly
    below the measured range. Sized by detector physics (noise & bias),
    not signal level. The polynomial near read 0 is locally linear
    (anchored by the implicit zero read), so extrapolating through this
    margin remains accurate.
    """

    order: int = 4
    modelName: str = "CHEBYSHEV"
    # Preliminary guess; needs a better-justified value derived from
    # detector noise/bias measurements.
    fitMinMargin: float = 100.0

    def __post_init__(self) -> None:
        if not isinstance(self.order, int) or isinstance(self.order, bool):
            raise ValueError(f"order must be an int, got {type(self.order).__name__}")
        if self.order < 1:
            raise ValueError(f"order must be >= 1, got {self.order}")
        if self.fitMinMargin < 0:
            raise ValueError(
                f"fitMinMargin must be >= 0, got {self.fitMinMargin}"
            )

    @staticmethod
    def _chebToMonomialMatrix(order: int) -> np.ndarray:
        """Return the (order+1, order+1) matrix that converts Chebyshev
        coefficients to monomial coefficients: ``mon = M @ cheb``.

        For order <= 5 the conversion is exact in float64 and numerically
        benign (cheb-to-poly amplification factors are O(2^order)).
        """
        M = np.zeros((order + 1, order + 1), dtype=np.float64)
        for k in range(order + 1):
            ek = np.zeros(order + 1)
            ek[k] = 1.0
            # cheb2poly returns a TRIMMED 1-D array (length 1 for T_0,
            # length 2 for T_1, ...), so slice the destination to match.
            polyCoefs = np.polynomial.chebyshev.cheb2poly(ek)
            M[:len(polyCoefs), k] = polyCoefs
        return M

    def chebToMonomial(self, chebCoefs: np.ndarray) -> np.ndarray:
        """Convert Chebyshev-basis coefficients to monomial (standard
        polynomial) coefficients on axis 0.

        Input shape ``(order+1, ...)``; output same shape. Caller does
        this once at apply time so the per-chunk Horner evaluation can
        skip the conversion — see :func:`evaluateMonomial`.
        """
        chebCoefs = np.asarray(chebCoefs)
        order = chebCoefs.shape[0] - 1
        if order == 0:
            return chebCoefs.copy()
        M = self._chebToMonomialMatrix(order).astype(chebCoefs.dtype)
        flat = chebCoefs.reshape(order + 1, -1)
        return (M @ flat).reshape(chebCoefs.shape)

    def evaluateMonomial(
        self, monCoefs: np.ndarray, x: np.ndarray
    ) -> np.ndarray:
        """Evaluate a monomial polynomial via Horner's method.

        ``t = mon[0] + mon[1]*x + mon[2]*x^2 + ...``

        Uses one ``(..., H, W)`` accumulator buffer (plus the input
        ``x``) — two cubes total, vs three for Clenshaw on the
        Chebyshev basis. For polynomial order <= 5 there's no numerical
        downside to switching to the monomial form, and we save one
        full-cube-sized buffer per chunk.

        Parameters
        ----------
        monCoefs
            Shape ``(order+1, H, W)``. ``monCoefs[k]`` is the coefficient
            of ``x^k``. Usually produced by :func:`chebToMonomial`.
        x
            Shape ``(..., H, W)``. Mapped input in [-1, 1].

        Returns
        -------
        t
            Same shape as ``x``.
        """
        monCoefs = np.asarray(monCoefs)
        x = np.asarray(x)
        order = monCoefs.shape[0] - 1
        if order == 0:
            return np.broadcast_to(monCoefs[0], x.shape).astype(x.dtype).copy()

        # Horner: acc = mon[order]; for k = order-1..0: acc = acc*x + mon[k]
        acc = np.broadcast_to(monCoefs[order, :, :], x.shape).astype(x.dtype).copy()
        for k in range(order - 1, -1, -1):
            acc *= x
            acc += monCoefs[k, :, :]
        return acc

    def evaluate(self, coefficients: np.ndarray, x: np.ndarray) -> np.ndarray:
        """Evaluate the per-pixel Chebyshev series.

        Internally converts Chebyshev coefficients to monomial form once
        (cheap: an ``(order+1, order+1)`` matrix multiply on the
        ``(order+1, H, W)`` coefficient cube) and evaluates with Horner.
        Saves one full-cube buffer vs. Clenshaw at the cost of a small
        upfront conversion. For tight loops over chunks of ``x`` with
        the same coefficients, call :func:`chebToMonomial` once and
        :func:`evaluateMonomial` per chunk to skip the conversion each
        time.

        Parameters
        ----------
        coefficients
            Shape ``(order+1, H, W)``, float32. ``coefficients[k]`` is the
            coefficient of T_k(x).
        x
            Shape ``(..., H, W)``, float32. Mapped input in [-1, 1].

        Returns
        -------
        t
            Same shape as ``x``.
        """
        monCoefs = self.chebToMonomial(coefficients).astype(x.dtype, copy=False)
        return self.evaluateMonomial(monCoefs, x)

    def isMonotonic(
        self,
        coefficients: np.ndarray,
        mMin: np.ndarray,
        mMax: np.ndarray,
        nSamples: int = 32,
    ) -> np.ndarray:
        """Return an ``(H, W)`` boolean map: ``True`` if the fit is monotonically
        increasing on ``[mMin, mMax]`` per pixel.

        Computes the derivative of the Chebyshev series, evaluates it at
        ``nSamples`` evenly-spaced points on the mapped interval ``[-1, 1]``,
        and checks that all sampled derivatives (in original m-space) are
        non-negative.
        """
        coefficients = np.asarray(coefficients, dtype=np.float64)
        order = coefficients.shape[0] - 1
        if order < 1:
            return np.ones(coefficients.shape[1:], dtype=bool)

        # Derivative of Chebyshev series: d/dx [Σ c_k T_k(x)] = Σ d_k T_k(x)
        # where the derivative coefficients satisfy the recurrence:
        #   d_{p-1} = 2 * p * c_p
        #   d_k     = d_{k+2} + 2 * (k+1) * c_{k+1}   for k = p-2, ..., 1
        #   d_0     = d_2 / 2 + c_1
        derivCoefs = np.zeros((order, *coefficients.shape[1:]), dtype=np.float64)
        derivCoefs[order - 1] = 2.0 * order * coefficients[order]
        for k in range(order - 2, 0, -1):
            dKplus2 = derivCoefs[k + 2] if k + 2 < order else np.zeros_like(derivCoefs[0])
            derivCoefs[k] = dKplus2 + 2.0 * (k + 1) * coefficients[k + 1]
        # d_0: handle the k+2 index carefully (it's 0 if order < 3)
        d2 = derivCoefs[2] if order >= 3 else np.zeros_like(derivCoefs[0])
        derivCoefs[0] = d2 / 2.0 + coefficients[1]

        H, W = mMin.shape
        # Sample x in [-1, 1]
        xSamples = np.linspace(-1.0, 1.0, nSamples, dtype=np.float64)[:, None, None]
        xSamples = np.broadcast_to(xSamples, (nSamples, H, W)).copy()

        # Evaluate derivative Chebyshev series at sample points via Clenshaw.
        derivOrder = order - 1
        if derivOrder == 0:
            d = np.broadcast_to(derivCoefs[0], (nSamples, H, W)).copy()
        else:
            bNext = np.zeros((nSamples, H, W), dtype=np.float64)
            bCurr = np.broadcast_to(
                derivCoefs[derivOrder], (nSamples, H, W)
            ).copy().astype(np.float64)
            for k in range(derivOrder - 1, 0, -1):
                bPrev = 2.0 * xSamples * bCurr - bNext + derivCoefs[k]
                bNext = bCurr
                bCurr = bPrev
            d = xSamples * bCurr - bNext + derivCoefs[0]

        # Chain rule: dt/dm = (dt/dx) * (dx/dm) = (dt/dx) * 2/(fitMax - fitMin)
        # For monotonicity we only care about sign, and 2/(fitMax - fitMin) > 0,
        # so we can just check dt/dx >= 0.
        allNonNegative = (d >= 0).all(axis=0)
        degenerate = mMax <= mMin
        return allNonNegative | degenerate

    def fitBlock(
        self,
        m: np.ndarray,
        t: np.ndarray,
        valid: np.ndarray,
        conditionNumberLimit: float,
    ) -> BlockFitResult:
        """Fit a Chebyshev polynomial at every pixel in the block.

        Maps ``m`` to ``x ∈ [-1, 1]`` via ``x = 2*(m - fitMin)/(fitMax - fitMin) - 1``
        before forming normal equations in the Chebyshev basis ``T_k(x)``.
        """
        nPoints, H, W = m.shape
        p = self.order
        nCoefs = p + 1

        mD = m.astype(np.float64)
        v64 = valid.astype(np.float64)
        t64 = t.astype(np.float64)

        # Count valid points per pixel
        nPointsUsed = valid.sum(axis=0).astype(np.int32)  # (H, W)
        badMask = np.zeros((H, W), dtype=np.uint8)

        # fitMin / fitMax: min/max of m over valid reads per pixel
        mMasked = np.where(valid, mD, np.nan)
        with np.errstate(invalid="ignore"):
            fitMin = np.nanmin(mMasked, axis=0)
            fitMax = np.nanmax(mMasked, axis=0)
        fitMin = np.where(np.isnan(fitMin), 0.0, fitMin)
        fitMax = np.where(np.isnan(fitMax), 0.0, fitMax)
        # Extend fitMin downward by an absolute DN margin so apply()
        # does not flag noise / bad bias-level offsets below the
        # measured range. The size of the margin is driven by detector
        # physics (read noise, kTC, bias variation), not signal level.
        # The polynomial near read 0 (anchored by the implicit zero
        # read) is locally linear, so extrapolating through this margin
        # stays accurate.
        fitMin = fitMin - self.fitMinMargin

        # Affine map m → x ∈ [-1, 1]: x = 2*(m - fitMin)/(fitMax - fitMin) - 1
        denom = fitMax - fitMin
        denom = np.where(denom > 0, denom, 1.0)  # avoid /0 for degenerate pixels
        x = 2.0 * (mD - fitMin[None]) / denom[None] - 1.0  # (N, H, W)

        # Flag insufficient-points pixels now.
        insufficientPixels = nPointsUsed < (nCoefs + 1)
        badMask[insufficientPixels] |= INSUFFICIENT_POINTS

        # Compute Chebyshev basis values T_k(x) via three-term recurrence.
        # nCoefs is small (typically 5), so storing all of them is fine.
        tCheb = []
        for k in range(nCoefs):
            if k == 0:
                tk = np.ones_like(x)
            elif k == 1:
                tk = x.copy()
            else:
                tk = 2.0 * x * tCheb[k - 1] - tCheb[k - 2]
            tCheb.append(tk)

        # Accumulate normal equations AtA and Atb
        AtA = np.zeros((H, W, nCoefs, nCoefs), dtype=np.float64)
        Atb = np.zeros((H, W, nCoefs), dtype=np.float64)

        for i in range(nCoefs):
            vTi = v64 * tCheb[i]  # (N, H, W)
            Atb[..., i] = (vTi * t64[:, None, None]).sum(axis=0)
            for j in range(i, nCoefs):
                val = (vTi * tCheb[j]).sum(axis=0)  # (H, W)
                AtA[..., i, j] = val
                if i != j:
                    AtA[..., j, i] = val

        # Condition number check
        with np.errstate(divide="ignore", invalid="ignore"):
            conditionNumber = np.linalg.cond(AtA)
        conditionNumber = np.nan_to_num(conditionNumber, nan=np.inf, posinf=np.inf)

        fitFailed = (~insufficientPixels) & (conditionNumber > conditionNumberLimit)
        badMask[fitFailed] |= FIT_FAILED

        skip = insufficientPixels | fitFailed
        identityBlock = np.eye(nCoefs, dtype=np.float64)
        AtA[skip] = identityBlock
        Atb[skip] = 0.0

        # Batched solve
        try:
            sol = np.linalg.solve(AtA, Atb[..., None])[..., 0]  # (H, W, nCoefs)
        except np.linalg.LinAlgError:
            sol = np.zeros((H, W, nCoefs), dtype=np.float64)
            for hi in range(H):
                for wi in range(W):
                    if skip[hi, wi]:
                        continue
                    try:
                        sol[hi, wi] = np.linalg.solve(AtA[hi, wi], Atb[hi, wi])
                    except np.linalg.LinAlgError:
                        badMask[hi, wi] |= FIT_FAILED
                        sol[hi, wi] = 0.0
            skip = skip | (badMask & FIT_FAILED != 0)

        # No unscaling needed — coefficients are in Chebyshev basis directly.
        coefficients = np.zeros((p + 1, H, W), dtype=np.float32)
        for k in range(nCoefs):
            coefficients[k] = sol[..., k].astype(np.float32)
        coefficients[:, skip] = 0.0

        # Residuals: evaluate fit at each read and compare to t.
        tPred = self.evaluate(coefficients, x.astype(np.float32))  # (N, H, W)
        residuals = (t[:, None, None].astype(np.float32) - tPred) * valid
        nForDiv = np.where(nPointsUsed > 0, nPointsUsed, 1).astype(np.float32)
        residualRms = np.sqrt((residuals ** 2).sum(axis=0) / nForDiv).astype(np.float32)
        maxAbsResidual = np.abs(residuals).max(axis=0).astype(np.float32)

        # Monotonicity check — retry non-monotonic pixels at lower orders.
        monotonic = self.isMonotonic(
            coefficients, fitMin.astype(np.float32), fitMax.astype(np.float32)
        )
        monotonic[skip] = False
        nonMono = (~skip) & (~monotonic)

        for retryOrder in range(p - 1, 0, -1):
            retryIdx = np.flatnonzero(nonMono.ravel())
            if len(retryIdx) == 0:
                break
            retryNCoefs = retryOrder + 1
            # Need at least retryNCoefs + 1 valid points.
            canRetry = nonMono & (nPointsUsed >= retryNCoefs + 1)
            retryIdx = np.flatnonzero(canRetry.ravel())
            if len(retryIdx) == 0:
                continue
            rr = retryIdx // W
            rc = retryIdx % W

            # Build normal equations for these pixels at the reduced order.
            rAtA = np.zeros((len(retryIdx), retryNCoefs, retryNCoefs), dtype=np.float64)
            rAtb = np.zeros((len(retryIdx), retryNCoefs), dtype=np.float64)
            for i in range(retryNCoefs):
                vTi = v64[:, rr, rc] * tCheb[i][:, rr, rc]  # (N, nRetry)
                rAtb[:, i] = (vTi * t64[:, None]).sum(axis=0)
                for j in range(i, retryNCoefs):
                    val = (vTi * tCheb[j][:, rr, rc]).sum(axis=0)
                    rAtA[:, i, j] = val
                    if i != j:
                        rAtA[:, j, i] = val

            with np.errstate(divide="ignore", invalid="ignore"):
                rCond = np.linalg.cond(rAtA)
            goodCond = np.isfinite(rCond) & (rCond <= conditionNumberLimit)

            # Solve good-condition pixels.
            solveIdx = np.flatnonzero(goodCond)
            if len(solveIdx) == 0:
                continue
            try:
                rSol = np.linalg.solve(rAtA[solveIdx], rAtb[solveIdx, :, None])[..., 0]
            except np.linalg.LinAlgError:
                continue

            # Vectorized monotonicity check across all good-conditioned
            # retry pixels. Layout the candidate batch as (nCoefs, 1, nGood)
            # so isMonotonic / evaluate can run once instead of per-pixel.
            nGood = len(solveIdx)
            goodPxR = rr[solveIdx]
            goodPxC = rc[solveIdx]
            trialCoefsBatch = np.zeros((nCoefs, 1, nGood), dtype=np.float32)
            trialCoefsBatch[:retryNCoefs, 0, :] = rSol.T.astype(np.float32, copy=False)

            fmBatch = fitMin[goodPxR, goodPxC].astype(np.float32, copy=False).reshape(1, nGood)
            fMBatch = fitMax[goodPxR, goodPxC].astype(np.float32, copy=False).reshape(1, nGood)
            goodMono = self.isMonotonic(trialCoefsBatch, fmBatch, fMBatch)[0]

            if not goodMono.any():
                continue

            accPxR = goodPxR[goodMono]
            accPxC = goodPxC[goodMono]
            acceptedCoefs = trialCoefsBatch[:, 0, goodMono]  # (nCoefs, nAcc)

            coefficients[:, accPxR, accPxC] = acceptedCoefs
            monotonic[accPxR, accPxC] = True
            nonMono[accPxR, accPxC] = False
            conditionNumber[accPxR, accPxC] = rCond[solveIdx[goodMono]]

            # Recompute residuals for accepted pixels in one batched evaluate.
            xAcc = x[:, accPxR, accPxC].astype(np.float32, copy=False)[:, None, :]
            coefForEval = acceptedCoefs[:, None, :]                     # (nCoefs, 1, nAcc)
            tPredAcc = self.evaluate(coefForEval, xAcc)[:, 0, :]         # (N, nAcc)
            vAcc = valid[:, accPxR, accPxC]
            resAcc = (t.astype(np.float32, copy=False)[:, None] - tPredAcc) * vAcc
            nAccPts = np.maximum(nPointsUsed[accPxR, accPxC], 1).astype(np.float32, copy=False)
            residualRms[accPxR, accPxC] = np.sqrt(
                (resAcc ** 2).sum(axis=0) / nAccPts
            )
            maxAbsResidual[accPxR, accPxC] = np.abs(resAcc).max(axis=0)

        badMask[nonMono] |= NON_MONOTONIC

        return BlockFitResult(
            coefficients=coefficients,
            fitMin=fitMin.astype(np.float32),
            fitMax=fitMax.astype(np.float32),
            residualRms=residualRms,
            maxAbsResidual=maxAbsResidual,
            nPointsUsed=nPointsUsed,
            conditionNumber=conditionNumber.astype(np.float32),
            monotonic=monotonic,
            badPixelMask=badMask,
        )

    def toFitsHdus(self, correction) -> list[fits.ImageHDU]:
        """Serialize model coefficients to a single ImageHDU named COEFFS."""
        hdu = fits.ImageHDU(data=correction.coefficients, name="COEFFS")
        hdu.header["ORDER"] = (self.order, "polynomial order")
        hdu.header["COMMENT"] = "COEFFS axis 0 is the Chebyshev coefficient index; C0 (T_0) first."
        return [hdu]

    @classmethod
    def fromFitsHdus(cls, hdus) -> tuple["PolynomialModel", np.ndarray]:
        """Reconstruct a PolynomialModel + coefficients from HDUs written by toFitsHdus."""
        coeffsHdu = None
        for hdu in hdus:
            if getattr(hdu, "name", "") == "COEFFS":
                coeffsHdu = hdu
                break
        if coeffsHdu is None:
            raise ValueError("No COEFFS HDU found in provided hdus")
        order = int(coeffsHdu.header["ORDER"])
        coefficients = np.asarray(coeffsHdu.data, dtype=np.float32)
        return cls(order=order), coefficients
