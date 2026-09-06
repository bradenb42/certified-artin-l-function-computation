"""Exact arithmetic in the cyclotomic fields Q(zeta_e).

An element is stored in normal form: coefficients on the power basis
1, zeta, ..., zeta^{phi(e)-1}, reduced modulo the cyclotomic polynomial.
Coefficients are integers or fractions, so equality of normal forms is
equality in the field and division by an integer is exact when it should be.
"""
from __future__ import annotations
from fractions import Fraction
from functools import lru_cache
from math import gcd
import cmath

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)

def _poly_divmod_exact(num, den):
    """Exact division of integer polynomials (lists, ascending), den monic."""
    num = list(num)
    q = [0] * (len(num) - len(den) + 1)
    for i in range(len(num) - len(den), -1, -1):
        c = num[i + len(den) - 1]
        q[i] = c
        if c:
            for j, d in enumerate(den):
                num[i + j] -= c * d
    assert all(x == 0 for x in num), "non-exact division"
    return q

@lru_cache(maxsize=None)
def cyclotomic_poly(e):
    """Phi_e as an ascending integer coefficient list."""
    p = [-1] + [0] * (e - 1) + [1]  # x^e - 1
    for d in divisors(e):
        if d < e:
            p = _poly_divmod_exact(p, cyclotomic_poly(d))
    return tuple(p)

def _reduce(coeffs, e):
    """Reduce an ascending coefficient list modulo Phi_e; return length phi(e)."""
    Phi = cyclotomic_poly(e)
    deg = len(Phi) - 1
    c = list(coeffs)
    if len(c) < deg:
        c += [0] * (deg - len(c))
    for i in range(len(c) - 1, deg - 1, -1):
        t = c[i]
        if t:
            for j in range(deg + 1):
                c[i - deg + j] -= t * Phi[j]
    return c[:deg]

_ZETA_CACHE = {}

def _zeta_power(e, k):
    """Normal form of zeta_e^k (cached: all powers computed by one-step reductions)."""
    if e not in _ZETA_CACHE:
        Phi = cyclotomic_poly(e)
        deg = len(Phi) - 1
        cur = [1] + [0] * (deg - 1)
        pw = [tuple(cur)]
        for _ in range(e - 1):
            top = cur[-1]
            nxt = [0] + cur[:-1]
            if top:
                for j in range(deg):
                    nxt[j] -= top * Phi[j]
            cur = nxt
            pw.append(tuple(cur))
        _ZETA_CACHE[e] = pw
    return list(_ZETA_CACHE[e][k])

def _norm(x):
    if isinstance(x, Fraction) and x.denominator == 1:
        return int(x.numerator)
    return x

class Cyc:
    __slots__ = ("e", "c")

    def __init__(self, e, coeffs=None, reduce=True):
        self.e = e
        if coeffs is None:
            coeffs = []
        if reduce:
            coeffs = _reduce(coeffs, e)
        self.c = tuple(_norm(x) for x in coeffs)

    # constructors
    @staticmethod
    def zero(e):
        return Cyc(e, [])

    @staticmethod
    def one(e):
        return Cyc(e, [1])

    @staticmethod
    def integer(e, a):
        return Cyc(e, [a])

    @staticmethod
    def zeta(e, k=1):
        k %= e
        if e > 3000:      # no power cache for large fields (memory); direct reduction
            return Cyc(e, [0] * k + [1])
        return Cyc(e, _zeta_power(e, k), reduce=False)

    @staticmethod
    def root_of_unity(e, o, j):
        """zeta_o^j as an element of Q(zeta_e); requires o | e."""
        assert e % o == 0
        return Cyc.zeta(e, (e // o) * j)

    def _coerce(self, other):
        if isinstance(other, Cyc):
            if other.e == self.e:
                return other
            E = self.e * other.e // gcd(self.e, other.e)
            return self.embed(E), other.embed(E)
        return Cyc(self.e, [other])

    def embed(self, E):
        """Image in Q(zeta_E), e | E, zeta_e -> zeta_E^{E/e}."""
        assert E % self.e == 0
        m = E // self.e
        coeffs = [0] * (m * (len(self.c) - 1) + 1) if self.c else []
        for j, a in enumerate(self.c):
            if a:
                coeffs[m * j] = a
        return Cyc(E, coeffs)

    def _pair(self, other):
        o = self._coerce(other)
        if isinstance(o, tuple):
            return o
        return self, o

    def __add__(self, other):
        a, b = self._pair(other)
        n = max(len(a.c), len(b.c))
        return Cyc(a.e, [(a.c[i] if i < len(a.c) else 0) + (b.c[i] if i < len(b.c) else 0) for i in range(n)], reduce=False)

    __radd__ = __add__

    def __neg__(self):
        return Cyc(self.e, [-x for x in self.c], reduce=False)

    def __sub__(self, other):
        return self + (-Cyc(self.e, [other]) if not isinstance(other, Cyc) else -other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if not isinstance(other, Cyc):
            return Cyc(self.e, [x * other for x in self.c], reduce=False)
        a, b = self._pair(other)
        if not a.c or not b.c:
            return Cyc.zero(a.e)
        prod = [0] * (len(a.c) + len(b.c) - 1)
        for i, x in enumerate(a.c):
            if x:
                for j, y in enumerate(b.c):
                    if y:
                        prod[i + j] += x * y
        return Cyc(a.e, prod)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Cyc):
            return self * other.inverse()
        return Cyc(self.e, [Fraction(x, 1) / other for x in self.c], reduce=False)

    def divexact(self, a):
        """Coordinatewise exact division by the integer a (error if inexact)."""
        out = []
        for x in self.c:
            q = Fraction(x) / a
            if q.denominator != 1:
                raise ArithmeticError("inexact division in Z[zeta]")
            out.append(int(q))
        return Cyc(self.e, out, reduce=False)

    def __pow__(self, k):
        r = Cyc.one(self.e)
        b = self
        while k:
            if k & 1:
                r = r * b
            b = b * b
            k >>= 1
        return r

    def __eq__(self, other):
        if not isinstance(other, Cyc):
            other = Cyc(self.e, [other])
        a, b = self._pair(other)
        return a.c == b.c

    def __hash__(self):
        return hash((self.e, self.c))

    def __repr__(self):
        return f"Cyc({self.e}, {list(self.c)})"

    def is_zero(self):
        return all(x == 0 for x in self.c)

    def is_rational(self):
        return all(x == 0 for x in self.c[1:])

    def rational(self):
        if not self.is_rational():
            raise ValueError("not rational")
        return Fraction(self.c[0]) if self.c else Fraction(0)

    def galois(self, t):
        """sigma_t: zeta -> zeta^t (t prime to e)."""
        assert gcd(t, self.e) == 1
        coeffs = [0] * self.e
        for j, a in enumerate(self.c):
            if a:
                coeffs[(j * t) % self.e] += a
        return Cyc(self.e, coeffs)

    def conj(self):
        return self.galois(self.e - 1)

    def stabilizer(self):
        """{t in (Z/e)^x: sigma_t(x) = x}."""
        return [t for t in range(1, self.e + 1) if gcd(t, self.e) == 1 and self.galois(t) == self]

    def conductor(self):
        """Smallest c | e with x in Q(zeta_c)."""
        for c in divisors(self.e):
            if all(self.galois(t) == self for t in range(1, self.e + 1)
                   if gcd(t, self.e) == 1 and t % c == 1 % c):
                return c
        return self.e

    def restrict(self, c):
        """Rewrite in Q(zeta_c) (requires x in Q(zeta_c), c | e)."""
        assert self.e % c == 0
        m = self.e // c
        # solve: find coefficients on zeta_c^j = zeta_e^{m j}
        # Use linear algebra over Q on the basis of Q(zeta_e).
        basis = [Cyc.zeta(self.e, m * j) for j in range(phi(c))]
        sol = _solve_in_basis(basis, self)
        if sol is None:
            raise ValueError("element not in the subfield")
        return Cyc(c, sol)

    def to_complex(self):
        z = cmath.exp(2j * cmath.pi / self.e)
        return sum(float(a) * z ** j for j, a in enumerate(self.c))

    def norm_bound(self):
        """max_tau |tau(x)| (computed numerically; only used for diagnostics)."""
        return max(abs(self.galois(t).to_complex()) for t in range(1, self.e + 1) if gcd(t, self.e) == 1)

    def inverse(self):
        """Inverse in Q(zeta_e) by solving x*y = 1 over Q."""
        d = phi(self.e)
        cols = [(self * Cyc.zeta(self.e, j)).c for j in range(d)]
        # matrix M with M[i][j] = coefficient i of x*zeta^j; solve M y = e_0
        M = [[Fraction(cols[j][i] if i < len(cols[j]) else 0) for j in range(d)] for i in range(d)]
        rhs = [Fraction(1 if i == 0 else 0) for i in range(d)]
        y = _solve_linear(M, rhs)
        if y is None:
            raise ZeroDivisionError("zero element")
        return Cyc(self.e, y, reduce=False)

    def to_json(self):
        return {"e": self.e, "c": [str(x) if isinstance(x, Fraction) else x for x in self.c]}

    @staticmethod
    def from_json(d):
        return Cyc(d["e"], [Fraction(x) if isinstance(x, str) else x for x in d["c"]], reduce=False)


def _solve_linear(M, rhs):
    """Solve square system over Q by Gaussian elimination; None if singular."""
    n = len(M)
    A = [row[:] + [rhs[i]] for i, row in enumerate(M)]
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        inv = 1 / A[col][col]
        A[col] = [x * inv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [A[i][n] for i in range(n)]

def _solve_in_basis(basis, x):
    """Coefficients of x in the Q-span of basis (elements of the same Q(zeta_e)); None if not in span."""
    e = x.e
    d = phi(e)
    m = len(basis)
    rows = [[Fraction(b.c[i] if i < len(b.c) else 0) for b in basis] + [Fraction(x.c[i] if i < len(x.c) else 0)] for i in range(d)]
    # row reduce (d x (m+1))
    pivots = []
    r = 0
    for col in range(m):
        piv = next((i for i in range(r, d) if rows[i][col] != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = 1 / rows[r][col]
        rows[r] = [v * inv for v in rows[r]]
        for i in range(d):
            if i != r and rows[i][col] != 0:
                f = rows[i][col]
                rows[i] = [v - f * w for v, w in zip(rows[i], rows[r])]
        pivots.append(col)
        r += 1
    for i in range(r, d):
        if rows[i][m] != 0:
            return None
    sol = [Fraction(0)] * m
    for i, col in enumerate(pivots):
        sol[col] = rows[i][m]
    return sol