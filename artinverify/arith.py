"""Arithmetic for the checker: cyclotomic fields, F_p polynomials, Sturm.
"""
from __future__ import annotations
from fractions import Fraction
from math import gcd
import cmath, math

# ---- cyclotomic polynomials and normal forms (own implementation)
_PHI = {}
def cyclo(e):
    if e in _PHI:
        return _PHI[e]
    p = [-1] + [0] * (e - 1) + [1]
    for d in range(1, e):
        if e % d == 0:
            q = cyclo(d)
            out = [0] * (len(p) - len(q) + 1)
            for i in range(len(p) - len(q), -1, -1):
                c = p[i + len(q) - 1]
                out[i] = c
                if c:
                    for j, y in enumerate(q):
                        p[i + j] -= c * y
            p = out
    _PHI[e] = p
    return p

def phi(e):
    return sum(1 for k in range(1, e + 1) if gcd(k, e) == 1)

class Z:
    """Element of Q(zeta_e) in normal form (coefficients Fractions/ints)."""
    __slots__ = ("e", "c")
    def __init__(self, e, c, red=True):
        self.e = e
        if red:
            c = self._reduce(list(c), e)
        self.c = tuple(int(x) if isinstance(x, Fraction) and x.denominator == 1 else x for x in c)
    @staticmethod
    def _reduce(c, e):
        P = cyclo(e); d = len(P) - 1
        if len(c) < d:
            c += [0] * (d - len(c))
        for i in range(len(c) - 1, d - 1, -1):
            t = c[i]
            if t:
                for j in range(d + 1):
                    c[i - d + j] -= t * P[j]
        return c[:d]
    @staticmethod
    def zeta(e, k=1):
        return Z(e, [0] * (k % e) + [1])
    @staticmethod
    def const(e, a):
        return Z(e, [a])
    def lift(self, E):
        m = E // self.e
        c = [0] * (m * (len(self.c) - 1) + 1)
        for j, a in enumerate(self.c):
            c[m * j] = a
        return Z(E, c)
    def _pair(self, o):
        if not isinstance(o, Z):
            return self, Z.const(self.e, o)
        if o.e == self.e:
            return self, o
        E = self.e * o.e // gcd(self.e, o.e)
        return self.lift(E), o.lift(E)
    def __add__(self, o):
        a, b = self._pair(o)
        return Z(a.e, [x + y for x, y in zip(a.c, b.c)], False)
    __radd__ = __add__
    def __neg__(self):
        return Z(self.e, [-x for x in self.c], False)
    def __sub__(self, o):
        a, b = self._pair(o)
        return Z(a.e, [x - y for x, y in zip(a.c, b.c)], False)
    def __mul__(self, o):
        if not isinstance(o, Z):
            return Z(self.e, [x * o for x in self.c], False)
        a, b = self._pair(o)
        prod = [0] * (2 * len(a.c))
        for i, x in enumerate(a.c):
            if x:
                for j, y in enumerate(b.c):
                    if y:
                        prod[i + j] += x * y
        return Z(a.e, prod)
    __rmul__ = __mul__
    def __truediv__(self, k):
        return Z(self.e, [Fraction(x) / k for x in self.c], False)
    def __eq__(self, o):
        a, b = self._pair(o)
        return a.c == b.c
    def galois(self, t):
        c = [0] * self.e
        for j, a in enumerate(self.c):
            if a:
                c[(j * t) % self.e] += a
        return Z(self.e, c)
    def conj(self):
        return self.galois(self.e - 1) if self.e > 1 else self
    def is_zero(self):
        return all(x == 0 for x in self.c)
    def rational(self):
        assert all(x == 0 for x in self.c[1:]), "not rational"
        return Fraction(self.c[0])
    def to_complex(self):
        w = cmath.exp(2j * cmath.pi / self.e)
        return sum(float(a) * w ** j for j, a in enumerate(self.c))
    @staticmethod
    def from_json(d):
        return Z(d["e"], [Fraction(x) if isinstance(x, str) else x for x in d["c"]], False)

def root_of_unity(E, o, j):
    assert E % o == 0
    return Z.zeta(E, (E // o) * j)

# ---- polynomials over F_p
def ptrim(a):
    while a and a[-1] == 0:
        a.pop()
    return a
def pmul(a, b, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % p
    return ptrim(r)
def pdivmod(a, b, p):
    a = list(a); b = ptrim(list(b))
    inv = pow(b[-1], -1, p)
    q = [0] * max(0, len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] * inv % p
        q[i] = c
        if c:
            for j, y in enumerate(b):
                a[i + j] = (a[i + j] - c * y) % p
    return ptrim(q), ptrim(a[:len(b) - 1])
def pmod(a, b, p):
    return pdivmod(a, b, p)[1]
def pgcd(a, b, p):
    a, b = ptrim(list(a)), ptrim(list(b))
    while b:
        a, b = b, pmod(a, b, p)
    if not a:
        return []
    inv = pow(a[-1], -1, p)
    return [x * inv % p for x in a]
def ppowmod(base, k, mod, p):
    r = [1]; base = pmod(base, mod, p)
    while k:
        if k & 1:
            r = pmod(pmul(r, base, p), mod, p)
        base = pmod(pmul(base, base, p), mod, p)
        k >>= 1
    return r
def factorization_type(f, p):
    """Degrees (with multiplicity) of the irreducible factors of f mod p, by distinct-degree
    factorization of the squarefree part and multiplicities by division."""
    f = ptrim([c % p for c in f])
    inv = pow(f[-1], -1, p); f = [c * inv % p for c in f]
    degs = []
    # multiplicities: repeatedly remove gcd with derivative
    def sqfree_part(g):
        dg = ptrim([(i * c) % p for i, c in enumerate(g)][1:])
        if not dg:
            return None
        return pdivmod(g, pgcd(g, dg, p), p)[0]
    rem = f
    while len(rem) > 1:
        s = sqfree_part(rem)
        if s is None:
            # rem is a p-th power: rem(x) = h(x^p)
            h = [rem[i] for i in range(0, len(rem), p)]
            sub = factorization_type(h, p)
            degs += [d for d in sub for _ in range(p)]
            return sorted(degs, reverse=True)
        # distinct-degree factorization of s
        g = s; x = [0, 1]; xp = x; d = 1
        parts = []
        while len(g) - 1 >= 2 * d:
            xp = ppowmod(xp, p, g, p)
            gd = pgcd(g, ptrim([(a - b) % p for a, b in zip(xp + [0] * 2, x + [0] * len(xp))]) or [], p) if True else []
            if len(gd) > 1:
                parts += [d] * ((len(gd) - 1) // d)
                g = pdivmod(g, gd, p)[0]
                xp = pmod(xp, g, p)
            d += 1
        if len(g) > 1:
            parts.append(len(g) - 1)
        degs += parts
        rem = pdivmod(rem, s, p)[0]
    return sorted(degs, reverse=True)

def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

# ---- Sturm
def sturm_count(f):
    def trim(a):
        while a and a[-1] == 0: a.pop()
        return a
    def rem(a, b):
        a = [Fraction(x) for x in a]; b = trim([Fraction(x) for x in b])
        while len(a) >= len(b) and trim(a):
            c = a[-1] / b[-1]; d = len(a) - len(b)
            for i, y in enumerate(b): a[i + d] -= c * y
            trim(a)
        return a
    seq = [trim([Fraction(c) for c in f])]
    seq.append(trim([i * c for i, c in enumerate(seq[0])][1:]))
    while trim(list(seq[-1])):
        r = [-x for x in rem(seq[-2], seq[-1])]
        if not trim(list(r)): break
        seq.append(r)
    def changes(sign):
        s = [1 if p[-1] * sign ** (len(p) - 1) > 0 else -1 for p in seq]
        return sum(1 for a, b in zip(s, s[1:]) if a != b)
    return changes(-1) - changes(1)

# ---- integer polynomial discriminant valuation helpers
def disc_mod_p_zero(f, p):
    """p | disc f  <=>  gcd(f, f') mod p nontrivial (f monic)."""
    fb = [c % p for c in f]
    dfb = ptrim([(i * c) % p for i, c in enumerate(fb)][1:])
    return (not dfb) or len(pgcd(fb, dfb, p)) > 1
