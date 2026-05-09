"""Sanity tests for ``PolynomialModel.{chebToMonomial,evaluateMonomial,evaluate}``.

Production reads Chebyshev coefficients from disk and evaluates them
per pixel inside ``h4Linearity.apply``. ``apply`` calls
``chebToMonomial`` once and ``evaluateMonomial`` per chunk — a
Horner-on-monomial that uses one fewer (N, H, W) buffer than the
previous Clenshaw path. The trade-off is numerical (monomial form is
ill-conditioned at high order). Order <= 5 is fine; these tests cover
that range and assert agreement with ``numpy.polynomial.chebyshev.chebval``
to float32 precision.
"""
import unittest

import numpy as np

import lsst.utils.tests
from lsst.obs.pfs.h4Linearity.models.polynomial import PolynomialModel


def _refChebval(coefs, x):
    """Per-pixel reference: numpy's chebval with coefficients on axis 0."""
    H, W = coefs.shape[1], coefs.shape[2]
    out = np.empty(x.shape, dtype=x.dtype)
    for h in range(H):
        for w in range(W):
            out[..., h, w] = np.polynomial.chebyshev.chebval(
                x[..., h, w], coefs[:, h, w]
            )
    return out


class ChebToMonomialTestCase(lsst.utils.tests.TestCase):

    def _matchesReference(self, order):
        rng = np.random.default_rng(order)
        H, W = 3, 4
        cheb = rng.standard_normal((order + 1, H, W)).astype(np.float32)
        mon = PolynomialModel(order=order).chebToMonomial(cheb)

        # mon should match per-pixel cheb2poly applied 1-D.
        for h in range(H):
            for w in range(W):
                refMon = np.polynomial.chebyshev.cheb2poly(cheb[:, h, w])
                actMon = mon[:len(refMon), h, w]
                np.testing.assert_allclose(actMon, refMon, rtol=1e-5, atol=1e-5)
                # Any trailing entries (when cheb2poly trims) should be zero.
                if mon.shape[0] > len(refMon):
                    np.testing.assert_array_equal(
                        mon[len(refMon):, h, w],
                        np.zeros(mon.shape[0] - len(refMon), dtype=mon.dtype),
                    )

    def testOrders1To5MatchNumpy(self):
        for order in (1, 2, 3, 4, 5):
            with self.subTest(order=order):
                self._matchesReference(order)

    def testOrder1IsIdentityForLinearTerm(self):
        # T_0 = 1, T_1 = x → for order=1, mon[0]=cheb[0], mon[1]=cheb[1].
        cheb = np.array(
            [[[1.5, 2.5]], [[3.0, 4.0]]], dtype=np.float32
        )  # (2, 1, 2)
        mon = PolynomialModel(order=1).chebToMonomial(cheb)
        np.testing.assert_array_equal(mon, cheb)


class EvaluateMonomialTestCase(lsst.utils.tests.TestCase):

    def _matchesChebvalAt(self, order, xShape):
        rng = np.random.default_rng(100 + order)
        H, W = 3, 4
        cheb = rng.standard_normal((order + 1, H, W)).astype(np.float32)
        x = rng.uniform(-1.0, 1.0, size=xShape + (H, W)).astype(np.float32)
        model = PolynomialModel(order=order)
        t = model.evaluate(cheb, x)
        tRef = _refChebval(cheb, x)
        np.testing.assert_allclose(t, tRef, rtol=1e-3, atol=1e-3)

    def testOrders1To5(self):
        for order in (1, 2, 3, 4, 5):
            with self.subTest(order=order):
                self._matchesChebvalAt(order, xShape=(7,))

    def testWorksOnSingleFrame(self):
        # x shape (H, W) — what apply()'s legacy applyFrame uses.
        self._matchesChebvalAt(order=4, xShape=())

    def testEvaluateMonomialMatchesEvaluate(self):
        """evaluate() == evaluateMonomial(chebToMonomial(c), x)."""
        rng = np.random.default_rng(7)
        order = 4
        H, W = 3, 4
        cheb = rng.standard_normal((order + 1, H, W)).astype(np.float32)
        x = rng.uniform(-1.0, 1.0, size=(5, H, W)).astype(np.float32)
        model = PolynomialModel(order=order)
        mon = model.chebToMonomial(cheb).astype(np.float32, copy=False)
        t1 = model.evaluate(cheb, x)
        t2 = model.evaluateMonomial(mon, x)
        np.testing.assert_array_equal(t1, t2)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
