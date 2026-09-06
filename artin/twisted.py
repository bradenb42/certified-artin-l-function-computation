"""Twisted resolvents: refinement inside a rational class.

For a cyclic subgroup C and a faithful character chi of C, the twisted
resolvent has coefficients in Z[zeta_o].  Counting the roots of its reduction
that are given powers of a root of unity yields the value of the induced
character on a power of the Frobenius, which separates classes that the
rational-class tests leave together.
"""
from __future__ import annotations
from math import gcd

from .cyclo import Cyc, cyclotomic_poly, phi
from .invariants import Invariant
from .perm import mul, inverse, identity, PermGroup
from .precision import Policy
from .resolvent import coset_transversal, roots_at
from .padic import GaloisRing, GF, hensel_lift, roots_in_GF, _pgcd, _ppowmod, _psub, _pmod, _pmul, _ptrim
from .ramified import factor_mod
from .fpoly import poly_mod, poly_mul

def multiplicative_order(a, m):
    a %= m
    k, x = 1, a
    while x != 1:
        x = x * a % m
        k += 1
    return k

def ring_with_root_of_unity(f, p, k, o, seed=0):
    """Galois ring at p of degree lcm(r_f, ord_o(p)) containing an element omega of exact
    order o (Hensel lift of a residue root of Phi_o), and the roots of f in it."""
    facs = factor_mod(f, p)
    r = 1
    for g, e in facs:
        d = len(g) - 1
        r = r * d // gcd(r, d)
    ro = multiplicative_order(p, o)
    R = r * ro // gcd(r, ro)
    ring, roots, _ = roots_at(f, p, k, seed=seed, r=R)
    Phi = [c % p for c in cyclotomic_poly(o)]
    w = roots_in_GF(Phi, ring.res)[0]
    omega = hensel_lift(ring, list(cyclotomic_poly(o)), w)
    return ring, roots, omega

def recover_cyclotomic(ring: GaloisRing, x, omega, o):
    """c in Z[zeta_o] with c(omega) = x mod p^k: solve over Z/p^k on the basis 1..omega^{phi-1}
    (the power basis of omega is independent modulo p), then centre; the size bound is covered by the policy."""
    d = phi(o)
    cols = [ring.one()]
    for _ in range(d - 1):
        cols.append(ring.mul(cols[-1], omega))
    mod, p = ring.mod, ring.p
    # rows: ring coordinates (r), columns: d unknowns, rhs x
    A = [[cols[j][i] for j in range(d)] + [x[i]] for i in range(ring.r)]
    piv_rows = []
    for c in range(d):
        piv = next((i for i in range(len(A)) if i not in piv_rows and A[i][c] % p), None)
        if piv is None:
            raise ValueError("power basis of omega not independent mod p")
        inv = pow(A[piv][c], -1, mod)
        A[piv] = [(v * inv) % mod for v in A[piv]]
        for i in range(len(A)):
            if i != piv and A[i][c]:
                fct = A[i][c]
                A[i] = [(v - fct * w) % mod for v, w in zip(A[i], A[piv])]
        piv_rows.append(piv)
    for i in range(len(A)):
        if i not in piv_rows and A[i][d] % mod:
            raise ValueError("element not in Z/p^k[omega]")
    coeffs = []
    for c in range(d):
        v = A[piv_rows[c]][d] % mod
        coeffs.append(v - mod if v > mod // 2 else v)
    return Cyc(o, coeffs)

class TwistedResolvent:
    """R_{C,chi_j}(x) = prod_{g in G/C}(x - (g Lambda)(alpha)^o), C = <sigma>, chi_j(sigma) = zeta_o^j."""
    def __init__(self, G: PermGroup, cl, k_class, j, f, p, policy: Policy, seed=0):
        self.G, self.cl, self.k_class, self.j = G, cl, k_class, j
        sigma = cl.reps[k_class]
        self.o = cl.orders[k_class]
        o = self.o
        C = PermGroup([sigma], n=G.n)
        self.C = C
        self.T = coset_transversal(G, C)
        self.base = G.base()
        # policy pair: a synthetic invariant with the height data of (g Lambda)^o
        # height data of the twisted invariant: ||.||_1 = phi(o)^{phi(o)/2} o^o, degree = o deg(m_t)
        deg_mt = sum(range(1, len(self.base) + 1))
        synth = Invariant(G.n, {tuple([deg_mt * o] + [0] * (G.n - 1)): int(phi(o) ** (phi(o) / 2 + 0.5)) * o ** o},
                          label=f"twisted R_(C=<g_{k_class+1}>, chi_{j})")
        policy.register(synth, len(self.T), label=synth.label + f", m={len(self.T)}")
        k = policy.k(p)
        ring, roots, omega = ring_with_root_of_unity(f, p, k, o, seed)
        self.ring_p = ring
        powers = [identity(G.n)]
        for _ in range(o - 1):
            powers.append(mul(powers[-1], sigma))
        def lam(g):
            # (g Lambda)(alpha) = sum_i omega^{-ij} (g sigma^i m_t)(alpha)
            tot = ring.zero()
            for i, s in enumerate(powers):
                h = mul(g, s)
                term = ring.one()
                for jj, t in enumerate(self.base):
                    term = ring.mul(term, ring.power(roots[h[t]], jj + 1))
                tot = ring.add(tot, ring.mul(term, ring.power(omega, (-i * j) % o)))
            return tot
        vals = [ring.power(lam(g), o) for g in self.T]
        # product, then recover coefficients in Z[zeta_o]
        poly = [ring.one()]
        for v in vals:
            new = [ring.zero() for _ in range(len(poly) + 1)]
            for i, c in enumerate(poly):
                new[i + 1] = ring.add(new[i + 1], c)
                new[i] = ring.sub(new[i], ring.mul(c, v))
            poly = new
        self.R = [recover_cyclotomic(ring, c, omega, o) for c in poly]
        # values of Ind_C^G chi_j on all classes, from the table-free formula
        self.ind = {}
        for kk in range(cl.r):
            self.ind[kk] = self._ind_value(kk)

    def _ind_value(self, kk):
        """Ind_C^G chi(y) = |C_G(y)|/|C| sum_{c in C cap y^G} chi(c)."""
        cl, o, j = self.cl, self.o, self.j
        s = Cyc.zero(o)
        c = identity(self.G.n)
        for i in range(o):
            if cl.class_of(c) == kk:
                s = s + Cyc.root_of_unity(o, o, (i * j) % o)
            c = mul(c, self.cl.reps[self.k_class])
        return s * (cl.centralizer_order(kk) // o)

    def counts(self, ell):
        """n_j' for all j' at a prime lambda | ell of Z[zeta_o]; None if ell is
        excluded (ell | o, or R_{C,chi} not squarefree mod lambda)."""
        o = self.o
        if o % ell == 0 or any(c.e != o for c in self.R):
            return None
        fo = multiplicative_order(ell, o)
        Phi = [c % ell for c in cyclotomic_poly(o)]
        facs = factor_mod(Phi, ell)
        g = next(gg for gg, e in facs if len(gg) - 1 == fo)
        K = GF(ell, g)
        w = K.from_int((-g[0]) % ell) if K.r == 1 else [0, 1] + [0] * (K.r - 2)
        def red(c):
            out = K.zero()
            pw = K.one()
            for a in c.c:
                if a:
                    out = K.add(out, [(x * int(a)) % ell for x in pw])
                pw = K.mul(pw, w)
            return out
        Rbar = _ptrim([red(c) for c in self.R], K)
        if len(Rbar) - 1 != len(self.T):
            return None
        # squarefree mod lambda?
        dR = _ptrim([K.mul(K.from_int(i), Rbar[i]) for i in range(1, len(Rbar))], K)
        if len(_pgcd(Rbar, dR, K)) > 1:
            return None
        N = K.q
        y = [K.zero(), K.one()]
        ypow = _ppowmod(y, (N - 1) // o, Rbar, K)
        counts = []
        for jj in range(o):
            wj = K.power(w, jj)
            h = _psub(ypow, [wj], K)
            gcdp = _pgcd(Rbar, h, K) if h else Rbar
            counts.append(len(gcdp) - 1)
        return fo, counts

    def ind_value_at(self, ell):
        """Ind_C^G chi(sigma_ell^{f_o}) as a Cyc, or None."""
        res = self.counts(ell)
        if res is None:
            return None
        fo, counts = res
        s = Cyc.zero(self.o)
        for jj, n in enumerate(counts):
            s = s + Cyc.root_of_unity(self.o, self.o, jj) * n
        return fo, s, counts

    def consistent_classes(self, ell, cands):
        """Classes x among cands with Ind_C^G chi(x^{f_o}) equal to the value read at ell."""
        r = self.ind_value_at(ell)
        if r is None:
            return None, None
        fo, val, counts = r
        keep = [x for x in cands if self.ind[self.cl.power_map(x, fo)] == val]
        return keep, {"f_o": fo, "counts": counts}