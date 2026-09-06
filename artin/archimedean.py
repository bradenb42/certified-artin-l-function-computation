"""The archimedean place: complex conjugation and Gamma factors.

Determines the class of complex conjugation from a certified count of the real
roots (Sturm sequences over the rationals) and, where needed, from the real
root counts of the separating resolvents, then computes the Gamma factor
exponents, the parity and the archimedean root number of every character.
"""
from __future__ import annotations
from fractions import Fraction
from math import gcd

from .cyclo import Cyc
from .euler import det_character
from .schur import multiplicity, Pair

# ------------------------------------------------------------------ Sturm

def _ptrim(a):
    while a and a[-1] == 0:
        a.pop()
    return a

def _prem(a, b):
    """Remainder of a by b over Q (ascending coefficient lists of Fractions)."""
    a = [Fraction(x) for x in a]
    b = _ptrim([Fraction(x) for x in b])
    while len(a) >= len(b) and _ptrim(a):
        c = a[-1] / b[-1]
        d = len(a) - len(b)
        for i, y in enumerate(b):
            a[i + d] -= c * y
        _ptrim(a)
    return a

def sturm_sequence(f):
    f0 = _ptrim([Fraction(c) for c in f])
    f1 = _ptrim([i * c for i, c in enumerate(f0)][1:])
    seq = [f0, f1]
    while _ptrim(list(seq[-1])):
        r = _prem(seq[-2], seq[-1])
        r = [-x for x in r]
        if not _ptrim(list(r)):
            break
        seq.append(r)
    return seq

def _eval(p, x):
    r = Fraction(0)
    for c in reversed(p):
        r = r * x + c
    return r

def _sign_changes(seq, x):
    signs = []
    for p in seq:
        v = _eval(p, x)
        if v != 0:
            signs.append(1 if v > 0 else -1)
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)

def _sign_changes_inf(seq, sign):
    signs = []
    for p in seq:
        lc = p[-1] * (sign ** (len(p) - 1))
        signs.append(1 if lc > 0 else -1)
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)

def real_root_count(f):
    """Number of distinct real roots of the squarefree integer polynomial f (Sturm, exact)."""
    seq = sturm_sequence(f)
    return _sign_changes_inf(seq, -1) - _sign_changes_inf(seq, 1)

def cauchy_bound(f):
    lead = abs(f[-1])
    return 1 + max(abs(Fraction(c, lead)) for c in f[:-1])

def isolate_real_roots(f):
    """Certified isolating intervals (a, b] with exactly one real root each (Sturm counts)."""
    seq = sturm_sequence(f)
    B = cauchy_bound(f)
    def count(a, b):
        return _sign_changes(seq, a) - _sign_changes(seq, b)
    out = []
    stack = [(Fraction(-B), Fraction(B))]
    while stack:
        a, b = stack.pop()
        c = count(a, b)
        if c == 0:
            continue
        if c == 1:
            out.append((a, b))
            continue
        m = (a + b) / 2
        stack.append((a, m)); stack.append((m, b))
    return sorted(out)

# ------------------------------------------------------------------ complex conjugation and the Gamma factors

def conjugation_class(G, cl, f, stage, log=print):
    """Class of complex conjugation: cycle type 1^r 2^{(n-r)/2} from the real root count r
    (Sturm), then separated inside its block by the real root counts of the separating resolvents
    (fix_{G/H}(c) = number of real roots of R_{G,H,F})."""
    n = G.n
    r = real_root_count(f)
    assert (n - r) % 2 == 0
    ct = tuple(sorted([2] * ((n - r) // 2) + [1] * r, reverse=True))
    cands = list(stage.blocks.get(ct, []))
    if not cands:
        raise RuntimeError("no class with the cycle type of complex conjugation")
    used = []
    for t in stage.tests:
        if len({t.sig[c] for c in cands}) <= 1:
            continue
        rr = real_root_count(t.R)
        used.append({"subgroup": t.S.label, "real_roots_of_R": rr})
        cands = [c for c in cands if sum(1 for x in t.sig[c] if x == 1) == rr]
    if len(cands) != 1:
        # involutions in one block not separated by fixed-point counts of the family: use all cyclic subgroups' resolvents
        raise RuntimeError(f"complex conjugation not separated: candidates {[c+1 for c in cands]}")
    log(f"complex conjugation: {r} real roots, cycle type {ct}, class {cands[0]+1}")
    return cands[0], r, used

def archimedean_data(table, cl, k_c, r, n):
    """a_chi, b_chi, Gamma exponents, epsilon_infinity, and the consistency checks."""
    out = []
    e = table.e
    for nu in range(table.r):
        d = table.degrees[nu]
        chic = table.values[nu][k_c]
        if not chic.is_rational():
            raise RuntimeError("chi(c) not rational")
        chic = chic.rational()
        if chic.denominator != 1 or (d - int(chic)) % 2:
            raise RuntimeError(f"chi(c) = {chic} not congruent to chi(1) mod 2")
        a, b = (d + int(chic)) // 2, (d - int(chic)) // 2
        if a < 0 or b < 0:
            raise RuntimeError("a or b negative")
        # parity: det rho_chi(c) = (-1)^b
        dv = det_character(table, nu)[k_c]
        detc = dv.rational()
        assert detc in (1, -1)
        if detc != (-1) ** b:
            raise RuntimeError(f"parity classification fails for chi_{nu+1}: det(c) = {detc}, b = {b}")
        out.append({"chi": nu + 1, "degree": d, "chi(c)": int(chic), "a": a, "b": b, "a+b=chi(1)": a + b == d,
                    "gamma": {"Gamma_R(s)": a, "Gamma_R(s+1)": b},
                    "epsilon_infinity": {0: "1", 1: "-i", 2: "-1", 3: "i"}[b % 4],
                    "parity": "even" if b % 2 == 0 else "odd",
                    "odd_two_dimensional": d == 2 and b == 1})
    # permutation character: sum n_chi b_chi = number of complex-conjugate pairs
    return out

def permutation_check(table, cl, data, G, r):
    n = G.n
    # permutation character on all n roots (all G-orbits): n_chi = <chi, pi>
    fix = [sum(1 for i in range(n) if cl.reps[k][i] == i) for k in range(cl.r)]
    def n_chi(nu):
        s = Cyc.zero(table.e)
        for k in range(cl.r):
            s = s + table.values[nu][k].conj() * (fix[k] * cl.sizes[k])
        q = s.rational() / G.order()
        assert q.denominator == 1
        return int(q)
    sb = sum(n_chi(nu) * data[nu]["b"] for nu in range(table.r))
    sa = sum(n_chi(nu) * data[nu]["a"] for nu in range(table.r))
    # regular representation: sum chi(1)(a - b) = sum chi(1) chi(c) = 0 unless c = 1
    reg = sum(table.degrees[nu] * (data[nu]["a"] - data[nu]["b"]) for nu in range(table.r))
    return {"sum_n_chi_b": sb, "complex_pairs": (n - r) // 2, "sum_n_chi_a": sa, "real_roots": r,
            "ok": sb == (n - r) // 2 and sa == (n + r) // 2, "regular_sum_chi1_chi(c)": reg}