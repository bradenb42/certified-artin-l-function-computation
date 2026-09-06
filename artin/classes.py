"""Rational conjugacy classes and the subgroups that separate them.

Groups the classes into rational classes (orbits under raising to powers
coprime to the order) and into blocks (equal cycle type), chooses a family of
subgroups separating the rational classes within each block by weighted
greedy set cover, and builds their resolvents.  At a prime, the factorisation
type of a resolvent gives the cycle type of the Frobenius on the cosets, and
the family's answers pin down the rational class of the Frobenius.
"""
from __future__ import annotations
from math import gcd

from .invariants import orbit_sum, set_product, stabilizer_in, Invariant
from .perm import PermGroup, mul, inverse, identity, cycle_type, to_json
from .precision import Policy
from .ramified import factor_mod, discriminant
from .resolvent import coset_transversal, roots_at, resolvent, evaluate_permuted, squarefree_certificate, divides_disc
from .padic import count_Zl_roots_policy

TSCHIRNHAUS = [None, [0, 1, 1], [0, 1, 2], [0, 2, 1], [1, 1, 1, 1], [0, 1, -1, 1], [0, 3, 1, 1], [0, 1, 0, 1]]

# ------------------------------------------------------------------ rational classes and blocks

def rational_classes(cl):
    """Partition of the classes into orbits of (Z/e)^x under k -> class of g_k^t."""
    e = cl.exponent
    parent = list(range(cl.r))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for k in range(cl.r):
        for t in range(1, e):
            if gcd(t, e) == 1:
                a, b = find(k), find(cl.power_map(k, t))
                if a != b:
                    parent[b] = a
    groups = {}
    for k in range(cl.r):
        groups.setdefault(find(k), []).append(k)
    return sorted(groups.values())

def blocks(cl):
    """Classes grouped by cycle type (the block located by factoring f mod ell)."""
    out = {}
    for k in range(cl.r):
        out.setdefault(cl.cycle_types[k], []).append(k)
    return out

def coset_action(G: PermGroup, H: PermGroup):
    """Transversal T of G/H and a function sending g in G to its permutation of the cosets."""
    T = coset_transversal(G, H)
    index = {}
    for c, x in enumerate(T):
        for h in H.elements():
            index[mul(x, h)] = c
    def act(g):
        return tuple(index[mul(g, x)] for x in T)
    return T, act

def signature(cl, act, k):
    """Cycle type of g_k on G/H (the H-signature)"""
    return cycle_type(act(cl.reps[k]))

# ------------------------------------------------------------------  family and cover

class Subgroup:
    def __init__(self, H: PermGroup, F: Invariant, label):
        self.H, self.F, self.label = H, F, label

def candidate_subgroups(G: PermGroup, cl):
    """Cyclic subgroups of the class representatives (one per rational class of cyclic
    subgroups), the base stabilizers, and the point and 2-set stabilizers, each with an
    invariant F having Stab_G(F) = H."""
    n = G.n
    base = G.base()
    out = []
    seen = set()
    for k in range(1, cl.r):
        g = cl.reps[k]
        H = PermGroup([g], n=n)
        key = frozenset(H.elements())
        if key in seen:
            continue
        seen.add(key)
        F = orbit_sum(n, H.elements(), base, label=f"orbit sum, H=<g_{k+1}>")
        out.append(Subgroup(H, F, {"kind": "cyclic", "class": k + 1}))
    for i in range(n):
        F = set_product(n, [i], label=f"x_{i+1}")
        st = stabilizer_in(F, G.elements())
        H = PermGroup(st, n=n) if len(st) > 1 else PermGroup([], n=n)
        if H.order() in (1, G.order()):
            continue
        key = frozenset(H.elements())
        if key in seen:
            continue
        seen.add(key)
        out.append(Subgroup(H, F, {"kind": "point stabilizer", "point": i + 1}))
        break   # one point suffices for a transitive G; others are conjugate
    for i in range(n):
        for j in range(i + 1, n):
            F = set_product(n, [i, j], label=f"x_{i+1}x_{j+1}")
            st = stabilizer_in(F, G.elements())
            H = PermGroup(st, n=n) if len(st) > 1 else PermGroup([], n=n)
            if H.order() in (1, G.order()):
                continue
            key = frozenset(H.elements())
            if key in seen:
                continue
            seen.add(key)
            out.append(Subgroup(H, F, {"kind": "2-set stabilizer", "set": [i + 1, j + 1]}))
    for lvl in range(1, len(base)):
        gens = G._level_gens(lvl)
        if not gens:
            continue
        H = PermGroup(gens, n=n)
        if H.order() in (1, G.order()):
            continue
        key = frozenset(H.elements())
        if key in seen:
            continue
        seen.add(key)
        F = orbit_sum(n, H.elements(), base, label=f"orbit sum, base stabilizer {lvl}")
        out.append(Subgroup(H, F, {"kind": "base stabilizer", "prefix": [b + 1 for b in base[:lvl]]}))
    return out

def greedy_cover(G, cl, subgroups):
    """Weighted set cover of the pairs of distinct rational classes in the same block,
    weight [G:H].  Returns the chosen subgroups (with their coset actions)
    and the signature tables."""
    rc = rational_classes(cl)
    rc_of = {k: i for i, orb in enumerate(rc) for k in orb}
    bl = blocks(cl)
    U = set()
    for lam, ks in bl.items():
        rcs = sorted({rc_of[k] for k in ks})
        for a in range(len(rcs)):
            for b in range(a + 1, len(rcs)):
                U.add((rcs[a], rcs[b]))
    data = []
    for S in subgroups:
        T, act = coset_action(G, S.H)
        sig = {k: signature(cl, act, k) for k in range(cl.r)}
        covered = {(a, b) for (a, b) in U if sig[rc[a][0]] != sig[rc[b][0]]}
        data.append((S, T, act, sig, covered))
    chosen, D = [], set()
    while D != U:
        best = None
        for d in data:
            new = d[4] - D
            if not new:
                continue
            ratio = len(d[1]) / len(new)
            if best is None or ratio < best[0]:
                best = (ratio, d)
        if best is None:
            raise RuntimeError("family does not separate all rational classes (cyclic subgroups missing?)")
        chosen.append(best[1])
        D |= best[1][4]
    return rc, rc_of, bl, U, chosen

# ------------------------------------------------------------------  resolvents

class ResolventTest:
    """R_{G,H,F} with its exact discriminant; per-prime evaluation."""
    def __init__(self, S: Subgroup, T, act, sig, f, p, policy: Policy, seed=0):
        self.S, self.T, self.act, self.sig = S, T, act, sig
        self.m = len(T)
        for Tr in TSCHIRNHAUS:
            Ft = Invariant(S.F.n, S.F.terms, Tr, S.F.label + (f" o T{Tr}" if Tr else ""))
            trial = Policy(f, policy.Delta, mode=policy.mode)
            trial.register(Ft, self.m)
            ring, roots, _ = roots_at(f, p, trial.k(p), seed=seed)
            R, _ = resolvent(Ft, T, roots, ring)
            cert = squarefree_certificate(R)
            if cert is not None:
                self.F, self.R, self.squarefree_prime = Ft, R, cert
                self.policy = policy
                policy.register(Ft, self.m, label=f"separating resolvent {Ft.label}, m={self.m}")
                return
        raise RuntimeError("no squarefree resolvent for " + str(S.label))

    def excluded(self, ell):
        return divides_disc(self.R, ell)

    def test(self, ell, k=None):
        """Data read at ell: (mode, value) with the factorization type mod ell (ell not excluded)
        or the number of Q_ell-roots (excluded ell, at the policy precision, T5)."""
        if not self.excluded(ell):
            facs = factor_mod(self.R, ell)
            ct = tuple(sorted((len(g) - 1 for g, e in facs for _ in range(e)), reverse=True))
            return "cycle type", ct
        cnt, _ = count_Zl_roots_policy(self.R, ell, self.policy, label=self.F.label)
        return "Q_ell roots", cnt

    def consistent(self, k, mode, value):
        if mode == "cycle type":
            return self.sig[k] == value
        return sum(1 for c in self.sig[k] if c == 1) == value

# ------------------------------------------------------------------ the rational-class stage

class RationalClassStage:
    def __init__(self, G, cl, f, p, policy, seed=0, log=print):
        self.G, self.cl, self.f, self.p, self.policy = G, cl, f, p, policy
        subs = candidate_subgroups(G, cl)
        self.rc, self.rc_of, self.blocks, self.U, chosen = greedy_cover(G, cl, subs)
        self.tests = [ResolventTest(S, T, act, sig, f, p, policy, seed) for (S, T, act, sig, cov) in chosen]
        log(f"separating family: {len(self.tests)} subgroups, weights {[t.m for t in self.tests]}, {len(self.U)} pairs to separate, rational classes {self.rc}")

    def block_of(self, ell):
        facs = factor_mod(self.f, ell)
        ct = tuple(sorted((len(g) - 1 for g, e in facs for _ in range(e)), reverse=True))
        return ct

    def candidates(self, ell, v_disc=None):
        """Classes consistent with the block and all resolvent tests at ell.  Returns
        (candidate classes, records)."""
        ct = self.block_of(ell)
        cands = list(self.blocks.get(ct, []))
        recs = []
        k = None
        for t in self.tests:
            # only tests whose signatures vary on the block are informative
            if len({t.sig[c] for c in cands}) <= 1:
                continue
            mode, val = t.test(ell)
            recs.append({"subgroup": t.S.label, "mode": mode, "value": list(val) if isinstance(val, tuple) else val})
            cands = [c for c in cands if t.consistent(c, mode, val)]
        return cands, ct, recs