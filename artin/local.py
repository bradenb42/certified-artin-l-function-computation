"""Tamely ramified local fields and their Galois action.

At a prime whose ramification indices are prime to it, the splitting field of
f over the p-adic numbers embeds in an explicit Galois extension
K = Q_{ell^F}(varpi) with varpi^E = ell.  The module builds K, finds the roots
of f in it by Panayi's method, and reads the decomposition and inertia groups
and the Frobenius coset off the action of the two generating automorphisms on
those roots.
"""
from __future__ import annotations
import random
from math import gcd

from .padic import GaloisRing, GF, irreducible_poly, roots_in_GF, hensel_lift, _cz, _ptrim, PrecisionExhausted
from .cyclo import cyclotomic_poly

def lcm(a, b):
    return a * b // gcd(a, b)

def mult_order(a, m):
    if m == 1:
        return 1
    assert gcd(a, m) == 1
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k

class TameField:
    """K = Q_{ell^F}(varpi), varpi^E = ell, at ell-adic precision k (varpi-adic precision E k).
    Elements: lists of E Galois-ring elements (coefficients of varpi^0..varpi^{E-1})."""
    def __init__(self, ell, k, E, F, seed=0):
        self.p, self.k, self.E, self.F = ell, k, E, F
        rng = random.Random(seed)
        g = irreducible_poly(ell, F, rng) if F > 1 else [0, 1]
        self.ring = GaloisRing(ell, k, g)
        self.res = self.ring.res
        self.prec = E * k
        # zeta_E in the Galois ring (E | ell^F - 1)
        if E > 1:
            Phi = [c % ell for c in cyclotomic_poly(E)]
            w = roots_in_GF(Phi, self.res, rng)[0]
            self.zeta = hensel_lift(self.ring, list(cyclotomic_poly(E)), w)
        else:
            self.zeta = self.ring.one()
        # Frobenius of the Galois ring: image of the generator
        if F > 1:
            xbar = [0, 1] + [0] * (F - 2)
            self.frob_x = hensel_lift(self.ring, [c % self.ring.mod for c in g], self.res.frob(xbar))
        else:
            self.frob_x = self.ring.one()

    # -- elements
    def zero(self): return [self.ring.zero() for _ in range(self.E)]
    def one(self):
        z = self.zero(); z[0] = self.ring.one(); return z
    def from_int(self, a):
        z = self.zero(); z[0] = self.ring.from_int(a); return z
    def from_ring(self, a):
        z = self.zero(); z[0] = list(a); return z
    def pi(self):
        z = self.zero()
        if self.E > 1:
            z[1] = self.ring.one()
        else:
            z[0] = self.ring.from_int(self.p)
        return z
    def add(self, a, b): return [self.ring.add(x, y) for x, y in zip(a, b)]
    def sub(self, a, b): return [self.ring.sub(x, y) for x, y in zip(a, b)]
    def neg(self, a): return [self.ring.sub(self.ring.zero(), x) for x in a]
    def mul(self, a, b):
        E, R = self.E, self.ring
        out = [R.zero() for _ in range(E)]
        for i, x in enumerate(a):
            if not any(x):
                continue
            for j, y in enumerate(b):
                if not any(y):
                    continue
                t = R.mul(x, y)
                if i + j < E:
                    out[i + j] = R.add(out[i + j], t)
                else:
                    out[i + j - E] = R.add(out[i + j - E], [(c * self.p) % R.mod for c in t])
        return out
    def scalar(self, a, c):
        return [[(x * c) % self.ring.mod for x in coeff] for coeff in a]
    def power(self, a, n):
        r = self.one()
        while n:
            if n & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            n >>= 1
        return r
    def val(self, a):
        """varpi-adic valuation (E * v_ell(a_i) + i minimized), self.prec if zero."""
        v = self.prec
        for i, x in enumerate(a):
            vx = self.ring.valuation(x)
            if vx < self.k:
                v = min(v, self.E * vx + i)
        return v
    def is_zero(self, a, prec=None):
        return self.val(a) >= (self.prec if prec is None else prec)
    def eq(self, a, b, prec=None):
        return self.is_zero(self.sub(a, b), prec)
    def residue(self, a):
        return self.ring.residue(a[0])
    def is_unit(self, a):
        return not self.res.is_zero(self.residue(a))
    def inv(self, a):
        assert self.is_unit(a)
        x = self.from_ring(self.ring.inv(a[0]))
        two = self.from_int(2)
        steps = 1
        while (1 << steps) < self.prec + 1:
            steps += 1
        for _ in range(steps + 1):
            x = self.mul(x, self.sub(two, self.mul(a, x)))
        return x
    def div_pi(self, a, s):
        """a / varpi^s, assuming varpi^s | a (loses precision; caller tracks it)."""
        E, R = self.E, self.ring
        out = [R.zero() for _ in range(E)]
        for i, x in enumerate(a):
            j = i - s
            q = 0
            while j < 0:
                j += E
                q += 1
            # divide x by ell^q exactly (up to the lost top digits)
            y = list(x)
            for _ in range(q):
                assert all(c % self.p == 0 for c in y), "not divisible by varpi^s"
                y = [c // self.p for c in y]
            out[j] = R.add(out[j], [c % R.mod for c in y])
        return out
    def eval_poly(self, coeffs, x):
        r = self.zero()
        for c in reversed(coeffs):
            r = self.add(self.mul(r, x), c)
        return r
    # -- Galois action
    def frob(self, a):
        """phi: Frobenius on the coefficients, varpi fixed."""
        return [self.ring.eval_poly(list(x), self.frob_x) for x in a]
    def tau(self, a):
        """tau: varpi -> zeta_E varpi."""
        out, z = [], self.ring.one()
        for x in a:
            out.append(self.ring.mul(x, z))
            z = self.ring.mul(z, self.zeta)
        return out

# ------------------------------------------------------------------ roots of f in K (Panayi)

def roots_in_field(f, K: TameField, prec=None, rng=random):
    """All roots of the integer polynomial f in O_K, to precision prec (varpi-adic).  Roots found
    below a shift (whose precision Panayi's recursion reduces) are refined by the general Newton
    iteration (Hensel with non-unit derivative) to the full precision."""
    if prec is None:
        prec = K.prec
    coeffs = [K.from_int(c) for c in f]
    roots = _panayi_K(coeffs, K, prec, rng)
    return [_newton_general(coeffs, x, K, prec) for x in roots]

def _newton_general(h, x, K, prec):
    """x <- x - h(x)/h'(x) with h'(x) = varpi^v u, u a unit: converges once v(h(x)) > 2v."""
    dh = [K.scalar(c, i) for i, c in enumerate(h)][1:]
    for _ in range(64):
        fx = K.eval_poly(h, x)
        vf = K.val(fx)
        if vf >= prec:
            return x
        d = K.eval_poly(dh, x)
        v = K.val(d)
        if vf <= 2 * v:
            return x     # not in the quadratic-convergence region: keep the Panayi root
        u = K.inv(K.div_pi(d, v))
        x = K.sub(x, K.mul(K.div_pi(fx, v), u))
    return x

def _residue_poly(coeffs, K):
    return _ptrim([K.residue(c) for c in coeffs], K.res)

def _panayi_K(h, K, prec, rng):
    if prec <= 0:
        raise PrecisionExhausted
    hbar = _residue_poly(h, K)
    if not hbar:
        raise PrecisionExhausted
    Kr = K.res
    out = []
    dbar = _ptrim([Kr.mul(Kr.from_int(i), c) for i, c in enumerate(hbar)][1:], Kr)
    _cz(hbar, Kr, out, rng)
    roots = []
    for a in out:
        da = _eval(dbar, a, Kr)
        A = K.from_ring([c % K.ring.mod for c in a])
        if not Kr.is_zero(da):
            roots.append(_newton(h, A, K, prec))
            continue
        # shift: h(a + varpi y)
        taylor = _taylor(h, A, K)
        g = []
        for j, t in enumerate(taylor):
            g.append(K.mul(t, K.power(K.pi(), j)))
        s = min(K.val(c) for c in g)
        if s >= prec:
            raise PrecisionExhausted
        g = [K.div_pi(c, s) for c in g]
        for y in _panayi_K(g, K, prec - s, rng):
            roots.append(K.add(A, K.mul(K.pi(), y)))
    return roots

def _eval(poly, a, Kr):
    r = Kr.zero()
    for c in reversed(poly):
        r = Kr.add(Kr.mul(r, a), c)
    return r

def _taylor(h, a, K):
    coeffs = list(h)
    out = []
    while coeffs:
        q, acc = [], K.zero()
        for c in reversed(coeffs):
            acc = K.add(K.mul(acc, a), c)
            q.append(acc)
        rem = q.pop()
        out.append(rem)
        coeffs = list(reversed(q))
    return out

def _newton(h, x, K, prec):
    dh = [K.scalar(c, i) for i, c in enumerate(h)][1:]
    steps = 1
    while (1 << steps) < prec + 1:
        steps += 1
    for _ in range(steps + 1):
        x = K.sub(x, K.mul(K.eval_poly(h, x), K.inv(K.eval_poly(dh, x))))
    assert K.val(K.eval_poly(h, x)) >= prec - K.E, "Newton did not converge"
    return x

# ------------------------------------------------------------------ the local Galois group

class LocalGalois:
    """Decomposition and inertia groups of f at a tame prime ell as permutations of the
    local roots beta_1..beta_n (numbered as found), with the Frobenius coset."""
    def __init__(self, f, ell, k, ram_data, seed=0, log=print):
        n = len(f) - 1
        E, F0 = 1, 1
        for e, fdeg in ram_data:
            E, F0 = lcm(E, e), lcm(F0, fdeg)
        if E % ell == 0:
            raise ValueError(f"wild ramification at {ell} (e = {E}): not handled by the tame construction")
        # F also large enough that mu_{2E} lies in Q_{ell^F} (needed by the tame root-number construction: (-ell)^{1/E} in K)
        F = lcm(F0, mult_order(ell, (2 * E if E % 2 == 0 else E) * (ell ** F0 - 1)))
        self.ell, self.E, self.F, self.k = ell, E, F, k
        K = TameField(ell, k, E, F, seed)
        self.K = K
        roots = roots_in_field(f, K, rng=random.Random(seed))
        if len(roots) != n:
            raise RuntimeError(f"found {len(roots)} roots of f in K at {ell}, expected {n}")
        # distinctness
        for i in range(n):
            for j in range(i + 1, n):
                assert not K.eq(roots[i], roots[j], K.prec - 2 * K.E), "roots not separated at this precision"
        self.roots = roots
        self.phi = self._perm(K.frob)
        self.tau = self._perm(K.tau)
        from .perm import PermGroup
        self.D = PermGroup([self.phi, self.tau], n=n)
        self.I = PermGroup([self.tau], n=n)
        self.e, self.f = self.I.order(), self.D.order() // self.I.order()
        log(f"local field at {ell}: K = Q_{{{ell}^{F}}}(varpi), varpi^{E} = {ell}; |D| = {self.D.order()}, e = {self.e}, f = {self.f}")

    def _perm(self, auto):
        K, n = self.K, len(self.roots)
        imgs = [auto(r) for r in self.roots]
        perm = []
        for im in imgs:
            j = next((i for i, r in enumerate(self.roots) if K.eq(im, r, K.prec - 2 * K.E)), None)
            if j is None:
                raise RuntimeError("automorphism image is not a root")
            perm.append(j)
        return tuple(perm)

    def evaluate(self, F, perm):
        """F(beta^perm) as an element of K (for the matching test)."""
        vals = [self.roots[perm[i]] for i in range(len(self.roots))]
        return F.evaluate(vals, self.K)

    def inertia_orbit_sizes(self):
        from .perm import cycles
        seen, sizes = set(), []
        for orb in _orbits(self.I, len(self.roots)):
            sizes.append(len(orb))
        return sorted(sizes)

def _orbits(G, n):
    seen, out = set(), []
    for i in range(n):
        if i in seen:
            continue
        orb, queue = {i}, [i]
        while queue:
            x = queue.pop()
            for g in G.generators:
                y = g[x]
                if y not in orb:
                    orb.add(y); queue.append(y)
        seen |= orb
        out.append(sorted(orb))
    return out