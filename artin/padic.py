"""Arithmetic modulo powers of a prime: finite fields, Galois rings, roots.

Provides F_{ell^r}, the Galois rings Z_{ell^r}/ell^k, Hensel lifting of simple
roots, the Frobenius permutation of a set of roots, and Panayi's exact count
of the roots of an integer polynomial in Z_ell.
"""
from __future__ import annotations
import random

from .fpoly import (poly_trim, poly_mul, poly_divmod, poly_mod, poly_gcd, poly_monic, poly_sub,
                    poly_add, poly_powmod, poly_eval, roots_mod_p)
from .ramified import factor_mod

# ------------------------------------------------------------------ F_{ell^r}

class GF:
    """F_{ell^r} = F_ell[x]/(g), elements as lists of length r."""
    def __init__(self, ell, g):
        self.p = ell
        self.g = list(g)
        self.r = len(g) - 1
        self.q = ell ** self.r

    def zero(self): return [0] * self.r
    def one(self): return [1] + [0] * (self.r - 1)
    def from_int(self, a): return [a % self.p] + [0] * (self.r - 1)
    def add(self, a, b): return [(x + y) % self.p for x, y in zip(a, b)]
    def sub(self, a, b): return [(x - y) % self.p for x, y in zip(a, b)]
    def neg(self, a): return [(-x) % self.p for x in a]
    def mul(self, a, b):
        m = poly_mod(poly_mul(a, b, self.p), self.g, self.p)
        return m + [0] * (self.r - len(m))
    def power(self, a, k):
        r = self.one()
        while k:
            if k & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            k >>= 1
        return r
    def inv(self, a):
        return self.power(a, self.q - 2)
    def is_zero(self, a): return not any(a)
    def eq(self, a, b): return all(x == y for x, y in zip(a, b))
    def frob(self, a): return self.power(a, self.p)
    def elements(self):
        import itertools
        for t in itertools.product(range(self.p), repeat=self.r):
            yield list(t)

def irreducible_poly(ell, r, rng=random):
    """A monic irreducible polynomial of degree r over F_ell (ascending)."""
    if r == 1:
        return [0, 1]
    while True:
        g = [rng.randrange(ell) for _ in range(r)] + [1]
        if _is_irreducible(g, ell):
            return g

def _is_irreducible(g, p):
    r = len(g) - 1
    x = [0, 1]
    xp = x
    for i in range(1, r // 2 + 1):
        xp = poly_powmod(xp, p, g, p)
        if len(poly_gcd(g, poly_sub(xp, x, p), p)) > 1:
            return False
    return True

def roots_in_GF(h, K: GF, rng=random):
    """Roots in K = F_{q} of h in F_ell[y] (coefficients ints) of degree d | r."""
    d = len(h) - 1
    if d == 1:
        return [K.from_int((-h[0]) % K.p)]
    if K.q <= 512:
        return [a for a in K.elements() if K.is_zero(_eval_over(h, a, K))]
    # Cantor-Zassenhaus over K for odd q
    H = [K.from_int(c) for c in h]
    out = []
    _cz(H, K, out, rng)
    return out

def _eval_over(h, a, K):
    r = K.zero()
    for c in reversed(h):
        r = K.add(K.mul(r, a), K.from_int(c) if isinstance(c, int) else c)
    return r

# polynomials over K: lists of K-elements
def _ptrim(a, K):
    while a and K.is_zero(a[-1]):
        a.pop()
    return a
def _pmul(a, b, K):
    if not a or not b:
        return []
    out = [K.zero() for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = K.add(out[i + j], K.mul(x, y))
    return _ptrim(out, K)
def _pdivmod(a, b, K):
    a = list(a); b = _ptrim(list(b), K)
    inv = K.inv(b[-1])
    q = [K.zero() for _ in range(max(0, len(a) - len(b) + 1))]
    for i in range(len(a) - len(b), -1, -1):
        c = K.mul(a[i + len(b) - 1], inv)
        q[i] = c
        for j, y in enumerate(b):
            a[i + j] = K.sub(a[i + j], K.mul(c, y))
    return _ptrim(q, K), _ptrim(a[:len(b) - 1], K)
def _pmod(a, b, K): return _pdivmod(a, b, K)[1]
def _pgcd(a, b, K):
    a, b = _ptrim(list(a), K), _ptrim(list(b), K)
    while b:
        a, b = b, _pmod(a, b, K)
    inv = K.inv(a[-1])
    return [K.mul(x, inv) for x in a]
def _ppowmod(base, k, mod, K):
    r = [K.one()]
    base = _pmod(base, mod, K)
    while k:
        if k & 1:
            r = _pmod(_pmul(r, base, K), mod, K)
        base = _pmod(_pmul(base, base, K), mod, K)
        k >>= 1
    return r
def _psub(a, b, K):
    n = max(len(a), len(b))
    return _ptrim([K.sub(a[i] if i < len(a) else K.zero(), b[i] if i < len(b) else K.zero()) for i in range(n)], K)

def _cz(g, K, out, rng):
    """Roots in K of g in K[y]: split g | y^q - y by Cantor-Zassenhaus (odd q) or by
    traces (q even)."""
    g = _ptrim(list(g), K)
    inv = K.inv(g[-1]); g = [K.mul(x, inv) for x in g]
    d = len(g) - 1
    if d <= 0:
        return
    if d == 1:
        out.append(K.neg(g[0]))
        return
    y = [K.zero(), K.one()]
    yq = _ppowmod(y, K.q, g, K)
    g = _pgcd(g, _psub(yq, y, K), K)
    d = len(g) - 1
    if d <= 1:
        return _cz(g, K, out, rng)
    while True:
        a = [rng.randrange(K.p) for _ in range(K.r)]
        if K.p % 2:
            h = _ppowmod([a, K.one()], (K.q - 1) // 2, g, K)
            h = _psub(h, [K.one()], K)
        else:
            # trace map Tr(a y) = sum_{i<r} (a y)^{2^i} splits y^q - y over F_{2^r}
            t = _pmod([K.zero(), a], g, K)
            h, cur = t, t
            for _ in range(K.r - 1):
                cur = _pmod(_pmul(cur, cur, K), g, K)
                h = _padd(h, cur, K)
        w = _pgcd(g, h, K) if h else g
        if 0 < len(w) - 1 < d:
            _cz(w, K, out, rng)
            _cz(_pdivmod(g, w, K)[0], K, out, rng)
            return

def _padd(a, b, K):
    n = max(len(a), len(b))
    return _ptrim([K.add(a[i] if i < len(a) else K.zero(), b[i] if i < len(b) else K.zero()) for i in range(n)], K)

# ------------------------------------------------------------------ Galois rings

class GaloisRing:
    """Z_{ell^r}/ell^k = (Z/ell^k)[x]/(g) with g a monic lift of an irreducible of degree r."""
    def __init__(self, ell, k, g):
        self.p, self.k, self.mod = ell, k, ell ** k
        self.g = list(g)
        self.r = len(g) - 1
        self.res = GF(ell, [c % ell for c in g])

    def zero(self): return [0] * self.r
    def one(self): return [1] + [0] * (self.r - 1)
    def from_int(self, a): return [a % self.mod] + [0] * (self.r - 1)
    def add(self, a, b): return [(x + y) % self.mod for x, y in zip(a, b)]
    def sub(self, a, b): return [(x - y) % self.mod for x, y in zip(a, b)]
    def mul(self, a, b):
        m = self._mod_g(_imul(a, b))
        return m
    def _mod_g(self, a):
        a = list(a)
        r = self.r
        for i in range(len(a) - 1, r - 1, -1):
            c = a[i] % self.mod
            if c:
                for j in range(r + 1):
                    a[i - r + j] = (a[i - r + j] - c * self.g[j]) % self.mod
        a = a[:r] + [0] * (r - len(a[:r]))
        return [x % self.mod for x in a]
    def power(self, a, k):
        r = self.one()
        while k:
            if k & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            k >>= 1
        return r
    def residue(self, a): return [x % self.p for x in a]
    def truncate(self, a, k2):
        m = self.p ** k2
        return [x % m for x in a]
    def is_unit(self, a): return not self.res.is_zero(self.residue(a))
    def inv(self, a):
        v = self.res.inv(self.residue(a)) + []
        v = [x % self.mod for x in v]
        prec = 1
        two = self.from_int(2)
        while prec < self.k:
            v = self.mul(v, self.sub(two, self.mul(a, v)))
            prec *= 2
        return v
    def eval_poly(self, coeffs, a):
        r = self.zero()
        for c in reversed(coeffs):
            r = self.add(self.mul(r, a), self.from_int(c))
        return r
    def valuation(self, a):
        """min v_ell over coordinates; k if zero."""
        v = self.k
        for x in a:
            x %= self.mod
            if x:
                t = 0
                while x % self.p == 0:
                    x //= self.p
                    t += 1
                v = min(v, t)
        return v
    def centered_int(self, a):
        """The centered integer representative if a is in Z/ell^k (constant), else None."""
        if any(x % self.mod for x in a[1:]):
            return None
        x = a[0] % self.mod
        return x - self.mod if x > self.mod // 2 else x

def _imul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out

def hensel_lift(ring: GaloisRing, f, a):
    """Lift a simple residue root a (in ring.res) of the integer polynomial f to precision ring.k."""
    x = [c % ring.mod for c in a]
    fp = [i * c for i, c in enumerate(f)][1:]
    prec = 1
    while prec < ring.k:
        d = ring.inv(ring.eval_poly(fp, x))
        x = ring.sub(x, ring.mul(ring.eval_poly(f, x), d))
        prec *= 2
    assert ring.valuation(ring.eval_poly(f, x)) >= ring.k
    return x

def lifted_roots(f, ell, k, r=None, rng=random):
    """All roots of f in Z_{ell^r}/ell^k that reduce to simple roots of f mod ell
    (all n roots when ell does not divide disc f).  Returns (ring, roots, factorization mod ell)."""
    facs = factor_mod(f, ell)
    if r is None:
        r = 1
        for g, e in facs:
            d = len(g) - 1
            r = r * d // __import__("math").gcd(r, d)
    g = irreducible_poly(ell, r, rng) if r > 1 else [0, 1]
    ring = GaloisRing(ell, k, g)
    roots = []
    for h, e in facs:
        if e != 1:
            continue
        d = len(h) - 1
        if r % d:
            continue
        for a in roots_in_GF(h, ring.res, rng):
            roots.append(hensel_lift(ring, f, a))
    return ring, roots, facs

def frobenius_perm(ring: GaloisRing, roots):
    """Permutation phi with residue(root_i)^ell = residue(root_phi(i))."""
    K = ring.res
    res = [ring.residue(x) for x in roots]
    perm = []
    for a in res:
        b = K.frob(a)
        j = next(i for i, c in enumerate(res) if K.eq(b, c))
        perm.append(j)
    return tuple(perm)

# ------------------------------------------------------------------ Panayi (T5)

class PrecisionExhausted(Exception):
    pass

def count_Zl_roots(h, ell, k):
    """Number of roots in Z_ell of the monic squarefree integer polynomial h,
    computed from h modulo ell^k by Panayi's recursion.  Raises PrecisionExhausted
    if k digits do not suffice (never at the policy precision, the precision policy"""
    h = [c % ell ** k for c in h]
    return _panayi(h, ell, k)

def _panayi(h, ell, k):
    if k <= 0:
        raise PrecisionExhausted
    mod = ell ** k
    hbar = poly_trim([c % ell for c in h])
    if not hbar:
        raise PrecisionExhausted
    dh = poly_trim([(i * c) % ell for i, c in enumerate(hbar)][1:])
    count = 0
    for a in roots_mod_p(hbar, ell):
        if poly_eval(dh, a, ell) % ell:
            count += 1          # simple residue root: Hensel gives exactly one Z_ell root
            continue
        # h(a + ell y) / ell^s
        n = len(h) - 1
        # Taylor coefficients of h at a by repeated synthetic division by (x - a)
        coeffs = list(h)
        taylor = []
        for j in range(n + 1):
            q, acc = [], 0
            for c in reversed(coeffs):
                acc = (acc * a + c) % mod
                q.append(acc)
            rem = q.pop()
            coeffs = list(reversed(q))
            taylor.append(rem)
            if not coeffs:
                break
        # h(a + ell y) = sum_j taylor[j] (ell y)^j
        g = [(taylor[j] * pow(ell, j, mod)) % mod for j in range(len(taylor))]
        s = min((_val(c, ell, k) for c in g), default=k)
        if s >= k:
            raise PrecisionExhausted
        g = [c // ell ** s for c in g]
        count += _panayi(g, ell, k - s)
    return count

def _val(c, ell, k):
    if c == 0:
        return k
    v = 0
    while c % ell == 0:
        c //= ell
        v += 1
    return v


def count_Zl_roots_policy(h, ell, policy, label="R"):
    """T5 under the policy: start at policy.k(ell); Panayi returns only exact counts, so on
    PrecisionExhausted (k <= v_ell(disc h) is possible when disc h was not registered) the
    precision is doubled and the consultation recorded (policy.extra)."""
    k = policy.k(ell)
    while True:
        try:
            c = count_Zl_roots(h, ell, k)
            policy.consultations[ell] = max(policy.consultations.get(ell, 0), k)
            return c, k
        except PrecisionExhausted:
            k *= 2
            policy.extra.setdefault(str(ell), []).append({"label": label, "k": k})
            if k > 10 ** 6:
                raise