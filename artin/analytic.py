"""The approximate functional equation and the subfield zeta identities.

Builds the Dirichlet coefficients from the Euler factors, the test function as
a numerical Mellin inverse of the Gamma factor together with rigorous decay
bounds, the truncated theta sums and their tails, and the coefficients of the
Dedekind zeta function of a subfield from the residue degrees of its
resolvent.
"""
from __future__ import annotations
import cmath, math
from math import gcd

from .cyclo import Cyc
from .euler import cpoly_mul

# ------------------------------------------------------------------ coefficients

def euler_factor_coeffs(P_cyc):
    """Complex coefficients of P(T)."""
    return [c.to_complex() for c in P_cyc]

def local_series(P, kmax):
    """1/P(T) = sum b_k T^k up to T^kmax (P[0] = 1)."""
    b = [1.0 + 0j]
    for k in range(1, kmax + 1):
        s = 0j
        for i in range(1, min(k, len(P) - 1) + 1):
            s -= P[i] * b[k - i]
        b.append(s)
    return b

def dirichlet_coefficients(euler, X):
    """a_m for m <= X from the Euler factors euler[ell] = P_ell(T) (complex lists); ell not in euler is
    an error unless all primes <= X are present."""
    a = [0j] * (X + 1)
    a[1] = 1.0 + 0j
    primes = sorted(euler)
    # multiplicative build
    for ell in primes:
        if ell > X:
            continue
        kmax = int(math.log(X, ell)) + 1
        b = local_series(euler[ell], kmax)
        new = a[:]
        for m in range(1, X + 1):
            if a[m] == 0:
                continue
            pk, k = ell, 1
            while m * pk <= X:
                new[m * pk] += a[m] * b[k]
                pk *= ell; k += 1
        a = new
    return a

def divisor_bound(d, X):
    """d_d(m) = number of ways to write m as a product of d factors (bound |a_m| <= d_{chi(1)}(m))."""
    D = [1] * (X + 1)
    for _ in range(d - 1):
        E = [0] * (X + 1)
        for i in range(1, X + 1):
            for j in range(i, X + 1, i):
                E[j] += D[i]
        D = E
    return D

# ------------------------------------------------------------------ kernels

def gamma_factor(a, b, s):
    """gamma(s) = Gamma_R(s)^a Gamma_R(s+1)^b, Gamma_R(s) = pi^{-s/2} Gamma(s/2), via log-gamma."""
    import mpmath
    val = mpmath.mpf(1)
    lg = a * (-(s / 2) * mpmath.log(mpmath.pi) + mpmath.loggamma(s / 2)) + b * (-((s + 1) / 2) * mpmath.log(mpmath.pi) + mpmath.loggamma((s + 1) / 2))
    return mpmath.exp(lg)

class Kernel:
    """g(x) = (1/2 pi i) int_{Re s = c} gamma(s) x^{-s} ds by the trapezoid rule on a vertical line;
    the integrand is analytic in 0 < Re s < 2c so the error is O(e^{-2 pi c / h})."""
    def __init__(self, a, b, c=1.0, h=0.05, T=None):
        import mpmath
        mpmath.mp.dps = 30
        d = a + b
        self.a, self.b, self.c, self.h = a, b, c, h
        if T is None:
            T = max(30.0, 20.0 + 60.0 / d)   # |gamma(c+iu)| ~ e^{-pi d |u| / 4}
        self.T = T
        n = int(T / h)
        self.us = [k * h for k in range(-n, n + 1)]
        self.gs = [complex(gamma_factor(a, b, mpmath.mpc(c, u))) for u in self.us]
        # rigorous decay bounds |g(x)| <= B(c) x^{-c}, B(c) = (1/2 pi) int |gamma(c+iu)| du, for a grid of c
        self.bounds = []   # (c, log B(c))
        for cc in [1, 2, 3, 5, 8, 12, 18, 26, 36, 50, 70, 100]:
            hh = 0.1
            nn = int((self.T + 4 * cc) / hh)
            tot = mpmath.mpf(0)
            for k in range(-nn, nn + 1):
                tot += abs(gamma_factor(a, b, mpmath.mpc(cc, k * hh)))
            tot = tot * hh / (2 * mpmath.pi) * mpmath.mpf("1.02")
            self.bounds.append((cc, float(mpmath.log(tot))))
    def bound(self, x):
        lx = math.log(x)
        ex = min(lB - c * lx for c, lB in self.bounds)
        return float("inf") if ex > 700 else math.exp(ex)

    def __call__(self, x):
        if x <= 0:
            raise ValueError
        if not hasattr(self, "_cache"):
            self._cache = {}
        v = self._cache.get(x)
        if v is not None:
            return v
        v = self._eval(x)
        if len(self._cache) < 400000:
            self._cache[x] = v
        return v

    def _eval(self, x):
        lx = math.log(x)
        s = 0j
        for u, gval in zip(self.us, self.gs):
            s += gval * cmath.exp(-(self.c + 1j * u) * lx)
        return (s * self.h / (2 * math.pi)).real

def closed_form(a, b):
    """g for degree <= 2 (cross-check of the numerical kernel)."""
    import mpmath
    if (a, b) == (1, 0): return lambda x: 2 * math.exp(-math.pi * x * x)
    if (a, b) == (0, 1): return lambda x: 2 * x * math.exp(-math.pi * x * x)
    if (a, b) == (1, 1): return lambda x: 2 * math.exp(-2 * math.pi * x)
    if (a, b) == (2, 0): return lambda x: 4 * float(mpmath.besselk(0, 2 * math.pi * x))
    if (a, b) == (0, 2): return lambda x: 4 * x * float(mpmath.besselk(0, 2 * math.pi * x))
    return None

# ------------------------------------------------------------------ the test

def theta(coeffs, g, sqrt_f, t, X):
    """sum_{m <= X} a_m g(m t / sqrt f); terms beyond the point where the rigorous kernel bound
    times m^6 is below 1e-24 are dropped (their total is far below double precision)."""
    s = 0j
    has_bound = hasattr(g, "bound")
    for m in range(1, X + 1):
        x = m * t / sqrt_f
        if has_bound and x > 3 and m % 16 == 0 and g.bound(x) * m ** 6 < 1e-24:
            break
        s += coeffs[m] * g(x)
    return s

def tail_estimate(d, g, sqrt_f, t, X, factor=6):
    """Upper bound for sum_{m > X} |a_m| |g(m t / sqrt f)| with |a_m| <= d_d(m) and the rigorous
    kernel bounds of Kernel.bound: exact divisor counts up to factor*X, then d_d(m) <= m^{d-1}
    and the largest c of the bound grid for the remainder."""
    D = divisor_bound(d, factor * X)
    s = sum(D[m] * g.bound(m * t / sqrt_f) for m in range(X + 1, factor * X + 1))
    c, lB = g.bounds[-1]
    M0 = factor * X
    if c > d:
        ex = lB - c * math.log(t / sqrt_f) + (d - c) * math.log(M0) - math.log(c - d)
        s += float("inf") if ex > 700 else math.exp(ex)
    return s

def functional_equation_defect(coeffs, coeffs_bar, g, sqrt_f, W, t, X):
    lhs = theta(coeffs, g, sqrt_f, t, X)
    rhs = W * theta(coeffs_bar, g, sqrt_f, 1 / t, X) / t
    return lhs, rhs, abs(lhs - rhs)

# ------------------------------------------------------------------ the subfield zeta identities

def zeta_subfield_coefficients(residue_degrees_by_prime, X):
    """b_m = number of ideals of norm m in N^H, from the residue degrees f(p/ell) of the primes
    above each ell <= X (given as lists)."""
    b = [0] * (X + 1)
    b[1] = 1
    for ell, fs in sorted(residue_degrees_by_prime.items()):
        if ell > X:
            continue
        # local factor prod (1 - T^f)^{-1}: number of ideals of norm ell^k
        kmax = int(math.log(X, ell)) + 1
        loc = [0] * (kmax + 1); loc[0] = 1
        for f in fs:
            new = loc[:]
            for k in range(f, kmax + 1):
                new[k] += new[k - f]
            loc = new
        new = b[:]
        for m in range(1, X + 1):
            if b[m] == 0:
                continue
            pk, k = ell, 1
            while m * pk <= X:
                new[m * pk] += b[m] * loc[k]
                pk *= ell; k += 1
        b = new
    return b