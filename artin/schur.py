"""Schur indices and matrix models of the irreducible characters.

Computes bounds on the Schur index (the Frobenius-Schur indicator, the
divisibility bounds, and the bound from every pair (H, lambda) of a subgroup
and a linear character), reporting an exact value when the bounds meet.  Also
builds a matrix model of a multiple of chi over a cyclotomic field, with a
common denominator cleared, by projecting an induced representation.

None of the downstream L-function data depends on these: they are recorded
for reference and used as an independent route to the Euler factors.
"""
from __future__ import annotations
from fractions import Fraction
from math import gcd

from .cyclo import Cyc, phi, divisors
from .perm import mul, inverse, power, identity, cycles


def lcm(a, b):
    return a * b // gcd(a, b)


def roots_of_unity_in_field(e, stab):
    """Order of the group of roots of unity of the fixed field of stab <= (Z/e)^x in Q(zeta_e)."""
    J = [j for j in range(e) if all(((t - 1) * j) % e == 0 for t in stab)]
    n = len(J)  # zeta_e^J is cyclic of this order
    return n if n % 2 == 0 else 2 * n


def eigen_multiplicities(table, nu, k):
    """m_j = multiplicity of eigenvalue zeta_o^j of rho_chi(g_k), j in Z/o (exact, from the table)."""
    cl, e = table.cl, table.e
    o = cl.orders[k]
    out = []
    for j in range(o):
        s = Cyc.zero(e)
        for t in range(o):
            s = s + table.values[nu][cl.power_map(k, t)] * Cyc.root_of_unity(e, o, (-j * t) % o)
        q = s.rational() / o
        assert q.denominator == 1
        out.append(int(q))
    return out


class Pair:
    """A pair (H, lambda): H a subgroup given by its elements, lambda a linear character
    given as a function on H with values in Q(zeta_{o_lam})."""
    def __init__(self, kind, H_elements, lam, o_lam, label):
        self.kind = kind            # "cyclic" or "stabilizer"
        self.H = H_elements         # list of permutations
        self.lam = lam              # dict elt -> exponent t (lambda(h) = zeta_{o_lam}^t)
        self.o_lam = o_lam
        self.label = label          # json-able description
        self._counts = None

def candidate_pairs(G, cl, include_stabilizers=True):
    """Cyclic subgroups <g_k> with all linear characters, and the pointwise
    stabilizers G_{(b_1..b_i)} of the base with the trivial character."""
    pairs = []
    n = G.n
    for k in range(cl.r):
        g = cl.reps[k]
        o = cl.orders[k]
        elts, h = [], identity(n)
        powers = {}
        for t in range(o):
            powers[h] = t
            elts.append(h)
            h = mul(h, g)
        for j in range(o):
            g0 = gcd(o, j)
            o_lam, j0 = o // g0, j // g0
            lam = {x: (j0 * t) % o_lam for x, t in powers.items()}
            pairs.append(Pair("cyclic", elts, lam, o_lam, {"kind": "cyclic", "class": k + 1, "lambda_exponent": j}))
    if include_stabilizers:
        from .perm import PermGroup
        for lvl in range(1, len(G._base)):
            gens = G._level_gens(lvl)
            H = PermGroup(gens, n=n) if gens else None
            if H is None or H.order() == 1 or H.order() == G.order():
                continue
            elts = H.elements()
            label = {"kind": "stabilizer", "base_prefix": [b + 1 for b in G._base[:lvl]], "order": H.order()}
            pairs.append(Pair("stabilizer", elts, {x: 0 for x in elts}, 1, dict(label, character="trivial")))
            # the restriction of the sign character of S_n, when nontrivial on H
            sgn = {x: (sum(len(c) - 1 for c in cycles(x)) % 2) for x in elts}
            if any(sgn.values()):
                pairs.append(Pair("stabilizer", elts, sgn, 2, dict(label, character="sign")))
    return pairs

def multiplicity(table, nu, pair):
    """<chi_H, lambda> = |H|^{-1} sum_h chi(h) conj(lambda(h)), by class sums: the counts
    #{h in H: class(h) = k, lambda(h) = zeta^t} are computed once per pair and reused for
    every character (the element enumeration is done a single time)."""
    cl, e = table.cl, table.e
    E = lcm(e, pair.o_lam)
    if not hasattr(pair, "_counts") or pair._counts is None:
        counts = {}
        for h in pair.H:
            key = (cl.class_of(h), pair.lam[h] % pair.o_lam)
            counts[key] = counts.get(key, 0) + 1
        pair._counts = counts
    s = Cyc.zero(E)
    for (k, t), cnt in pair._counts.items():
        s = s + table.values[nu][k].embed(E) * Cyc.root_of_unity(E, pair.o_lam, (-t) % pair.o_lam) * cnt
    q = s.rational() / len(pair.H)
    assert q.denominator == 1
    return int(q)

def schur_data(table, nu, pairs=None):
    cl, e = table.cl, table.e
    d = table.degrees[nu]
    stab = table.stabilizers[nu]
    units = [t for t in range(1, e + 1) if gcd(t, e) == 1]
    K_deg = len(units) // len(stab)
    ext_deg = len(stab)  # [Q(zeta_e):K]
    w = roots_of_unity_in_field(e, stab)
    m_inf = 2 if (table.is_real[nu] and table.indicator[nu] == -1) else 1
    bound = gcd(gcd(d, ext_deg), w)
    if pairs is None:
        pairs = candidate_pairs(table.cl.G, cl)
    best = None  # (bound_value, a, ext, index, pair)
    for pr in pairs:
        a = multiplicity(table, nu, pr)
        if a == 0:
            continue
        o_lam = pr.o_lam
        # [K(lambda):K]: Gal(Q(zeta_e)/K(lambda)) = Stab(chi) meet Stab(lambda); o_lam | e
        assert e % o_lam == 0
        stab_lam = [t for t in stab if t % o_lam == 1 % o_lam]
        ext = len(stab) // len(stab_lam)
        bound = gcd(bound, a * ext)          #  for every pair, not only the best
        cand = (a * ext, a, ext, table.order // len(pr.H), pr)
        if best is None or cand[:4] < best[:4]:
            best = cand
    if bound % m_inf:
        raise RuntimeError("inconsistent Schur index bounds")
    candidates = [m for m in divisors(bound) if m % m_inf == 0]
    exact = len(candidates) == 1
    return {"degree": d, "field_degree": K_deg, "field_conductor": table.conductors[nu],
            "roots_of_unity_in_K": w, "m_infinity": m_inf,
            "upper_bound": bound, "candidates": candidates,
            "value": candidates[0] if exact else None,
            "status": "exact" if exact else "bounded (local indices m_ell not computed)",
            "best_pair": {"pair": best[4].label, "multiplicity": best[1],
                          "field_extension_degree": best[2]}}


# ------------------------------------------------------------------- models

def _coset_data(G, pair):
    """Cosets xH: returns (reps, lookup elt -> (coset index, lambda-exponent of h) with elt = x h)."""
    lookup = {}
    reps = []
    for x in G.elements():
        if x in lookup:
            continue
        c = len(reps)
        reps.append(x)
        for h in pair.H:
            lookup[mul(x, h)] = (c, pair.lam[h])
    return reps, lookup


def induced_action(G, pair, E):
    """Monomial matrices of Ind_H^G lambda as sparse maps col -> (row, Cyc in Q(zeta_E))."""
    reps, lookup = _coset_data(G, pair)
    cache = [Cyc.root_of_unity(E, pair.o_lam, t) for t in range(pair.o_lam)]
    def act(g):
        out = {}
        for c, x in enumerate(reps):
            c2, t = lookup[mul(g, x)]
            out[c] = (c2, cache[t])
        return out
    return reps, act


def build_model(G, table, nu, pairs=None, max_dim=80, max_work=3_000_000):
    """Cleared-denominator model of a*chi (see module docstring). Returns dict or None."""
    cl, e = table.cl, table.e
    d = table.degrees[nu]
    c = table.conductors[nu]
    if pairs is None:
        pairs = candidate_pairs(G, cl)
    # choose pair: minimal multiplicity a, then minimal index [G:H], then smallest field
    best = None
    for pr in pairs:
        a = multiplicity(table, nu, pr)
        if a == 0:
            continue
        idx = G.order() // len(pr.H)
        cand = (a, idx, lcm(c, pr.o_lam), pr)
        if best is None or cand[:3] < best[:3]:
            best = cand
    a, idx, F, pair = best
    if idx > max_dim or G.order() * idx > max_work:
        return {"status": "skipped", "reason": f"index {idx} exceeds model_dim_limit or work bound",
                "multiplicity": a, "pair": pair.label}
    reps, act = induced_action(G, pair, F)
    m = len(reps)
    # projector P = sum_g chi(g^{-1}) rho(g), dense m x m over Z[zeta_F]
    P = [[Cyc.zero(F) for _ in range(m)] for _ in range(m)]
    elts = G.elements()
    chi_inv = {}
    for g in elts:
        kk = cl.class_of(inverse(g))
        if kk not in chi_inv:
            chi_inv[kk] = table.values[nu][kk].restrict(c).embed(F)
        coef = chi_inv[kk]
        if coef.is_zero():
            continue
        A = act(g)
        for col, (row, val) in A.items():
            P[row][col] = P[row][col] + val * coef
    # column echelon basis of the image (rank should be a*d)
    target_rank = a * d
    basis = []   # list of (pivot, vector) with vector reduced
    for col in range(m):
        v = [P[r][col] for r in range(m)]
        for piv, b in basis:
            if not v[piv].is_zero():
                f = v[piv]
                v = [x - f * y for x, y in zip(v, b)]
        piv = next((r for r in range(m) if not v[r].is_zero()), None)
        if piv is None:
            continue
        inv = v[piv].inverse()
        v = [x * inv for x in v]
        # reduce earlier basis vectors at the new pivot
        basis = [(pv, [x - b[piv] * y for x, y in zip(b, v)]) for pv, b in basis]
        basis.append((piv, v))
        if len(basis) == target_rank:
            break
    if len(basis) != target_rank:
        return {"status": "failed", "reason": f"rank {len(basis)} != a*chi(1) = {target_rank}"}
    pivots = [pv for pv, _ in basis]
    vecs = [b for _, b in basis]
    def coords(v):
        """Coordinates of v in the basis; verifies v lies in the span."""
        co = [v[pv] for pv in pivots]
        w = [Cyc.zero(F) for _ in range(m)]
        for ci, b in zip(co, vecs):
            if not ci.is_zero():
                w = [x + ci * y for x, y in zip(w, b)]
        if any(not (x == y) for x, y in zip(v, w)):
            raise RuntimeError("image vector not in the span: W is not G-stable?")
        return co
    def rho(g):
        A = act(g)
        cols = []
        for b in vecs:
            img = [Cyc.zero(F) for _ in range(m)]
            for col, x in enumerate(b):
                if not x.is_zero():
                    row, val = A[col]
                    img[row] = img[row] + val * x
            cols.append(coords(img))
        # matrix with columns = images of basis vectors
        return [[cols[jj][ii] for jj in range(target_rank)] for ii in range(target_rank)]
    gens = {("gen", i): g for i, g in enumerate(G.generators)}
    mats = {key: rho(g) for key, g in gens.items()}
    # trace check on class representatives
    trace_ok = True
    for kk, g in enumerate(cl.reps):
        M = rho(g)
        tr = sum((M[i][i] for i in range(target_rank)), Cyc.zero(F))
        want = (table.values[nu][kk] * a)
        if not (tr == want):
            trace_ok = False
            break
    # clear denominators
    den = 1
    for M in mats.values():
        for row in M:
            for x in row:
                for co in x.c:
                    if isinstance(co, Fraction):
                        den = lcm(den, co.denominator)
    def clear(M):
        return [[(x * den) for x in row] for row in M]
    den_basis = 1
    for b in vecs:
        for x in b:
            for co in x.c:
                if isinstance(co, Fraction):
                    den_basis = lcm(den_basis, co.denominator)
    return {"status": "ok", "multiplicity": a, "dimension": target_rank, "field": F, "_rho": rho,
            "pair": pair.label, "H_order": len(pair.H), "index": m,
            "denominator": den,
            "generator_matrices": [[[x.to_json() for x in row] for row in clear(mats[("gen", i)])] for i in range(len(G.generators))],
            "basis_denominator": den_basis,
            "basis": [[(x * den_basis).to_json() for x in b] for b in vecs],
            "trace_check": trace_ok}