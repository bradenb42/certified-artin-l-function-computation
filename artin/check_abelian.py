"""Cross-check against class field theory for abelian groups.

Fields inside a cyclotomic field are defined by Gauss periods, run through the
pipeline, and compared character by character with the data computed directly
from Dirichlet characters: the matching of characters through the Frobenius
values, the conductor, the Euler factors at ramified primes, the parity and
the root number.
"""
from __future__ import annotations
import cmath, math, os, tempfile, json, time
from math import gcd

from .cyclo import Cyc, phi
from .perm import PermGroup
from .run import run_pipeline
from .certificate import load_json

def units(c):
    return [a for a in range(1, c) if gcd(a, c) == 1]

def subgroup_generated(c, gens):
    H = {1}
    frontier = [1]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = x * g % c
            if y not in H:
                H.add(y); frontier.append(y)
    return sorted(H)

def period_polynomial(c, H):
    """Minimal polynomial of eta_H = sum_{h in H} zeta_c^h (exact in Z[zeta_c]); None if not squarefree."""
    U = units(c)
    cosets, seen = [], set()
    for a in U:
        if a in seen:
            continue
        S = sorted({a * h % c for h in H})
        seen |= set(S)
        cosets.append(S)
    e = c
    periods = [sum((Cyc.zeta(e, s) for s in S), Cyc.zero(e)) for S in cosets]
    poly = [Cyc.one(e)]
    for eta in periods:
        new = [Cyc.zero(e) for _ in range(len(poly) + 1)]
        for i, coef in enumerate(poly):
            new[i + 1] = new[i + 1] + coef
            new[i] = new[i] - coef * eta
        poly = new
    coeffs = []
    for coef in poly:
        assert coef.is_rational()
        q = coef.rational()
        assert q.denominator == 1
        coeffs.append(int(q))
    # distinct periods?
    for i in range(len(periods)):
        for j in range(i + 1, len(periods)):
            if periods[i] == periods[j]:
                return None, cosets
    return coeffs, cosets

def regular_generators(c, H, cosets):
    """Generators of (Z/c)^x / H acting on the cosets (as permutations, 0-indexed)."""
    idx = {}
    for i, S in enumerate(cosets):
        for s in S:
            idx[s] = i
    gens = []
    for a in units(c):
        perm = tuple(idx[a * S[0] % c] for S in cosets)
        if perm != tuple(range(len(cosets))) and perm not in gens:
            gens.append(perm)
    G = PermGroup(gens, n=len(cosets))
    # reduce generator set
    small = []
    for g in gens:
        if not PermGroup(small, n=len(cosets)).contains(g) if small else True:
            small.append(g)
            if PermGroup(small, n=len(cosets)).order() == G.order():
                break
    return small, idx

def dirichlet_characters(c, H, cosets, idx):
    """Characters of (Z/c)^x trivial on H, as dicts a -> exponent k with psi(a) = zeta_n^k (n = |cosets|),
    together with their primitive conductors and values."""
    # quotient group A = (Z/c)^x / H; its characters via a Smith-like enumeration: use the pipeline's
    # own machinery? no: enumerate homomorphisms to roots of unity by brute force on generators.
    n = len(cosets)
    A = list(range(n))
    reps = [S[0] for S in cosets]
    def amul(i, j):
        return idx[reps[i] * reps[j] % c]
    # find generators of A and their orders (small group)
    gens_A, span = [], {idx[1]}
    for i in A:
        if i in span:
            continue
        gens_A.append(i)
        new = set(span)
        frontier = list(span)
        while frontier:
            x = frontier.pop()
            y = amul(x, i)
            while y not in new:
                new.add(y); frontier.append(y)
                y = amul(y, i)
        span = new
        if len(span) == n:
            break
    orders = []
    for g in gens_A:
        o, x = 1, g
        while x != idx[1]:
            x = amul(x, g); o += 1
        orders.append(o)
    # express each element as a word in the generators (BFS)
    word = {idx[1]: tuple([0] * len(gens_A))}
    frontier = [idx[1]]
    while frontier:
        x = frontier.pop()
        for t, g in enumerate(gens_A):
            y = amul(x, g)
            if y not in word:
                w = list(word[x]); w[t] += 1
                word[y] = tuple(w); frontier.append(y)
    chars = []
    import itertools
    for ks in itertools.product(*[range(o) for o in orders]):
        # candidate: psi(g_t) = zeta_{o_t}^{k_t}; must be consistent (relations); check homomorphism on all pairs
        def val(i):
            w = word[i]
            s = 0
            for t in range(len(gens_A)):
                s += ks[t] * w[t] * (n // orders[t])
            return s % n   # exponent of zeta_n
        ok = all((val(i) + val(j) - val(amul(i, j))) % n == 0 for i in A for j in A)
        if ok:
            table = {a: val(idx[a]) for a in units(c)}
            chars.append(table)
    chars = [dict(t) for t in {tuple(sorted(ch.items())) for ch in chars}]
    assert len(chars) == n, (len(chars), n)
    out = []
    for ch in chars:
        # primitive conductor: smallest f | c with psi(a) depending only on a mod f
        cond = c
        for fdiv in range(1, c + 1):
            if c % fdiv:
                continue
            if all(ch[a] == ch[b] for a in units(c) for b in units(c) if (a - b) % fdiv == 0):
                cond = fdiv
                break
        out.append({"exponents": ch, "n": n, "conductor": cond})
    return out

def cft_data(c, H, cosets, idx, X):
    chars = dirichlet_characters(c, H, cosets, idx)
    n = len(cosets)
    zn = lambda k: cmath.exp(2j * cmath.pi * k / n)
    for ch in chars:
        f = ch["conductor"]
        # primitive character values mod f: psi_0(a) for a prime to f: pick b = a mod c with b prime to c
        def psi0(a, f=f, ch=ch):
            if f == 1:
                return 1 + 0j
            a %= f
            if gcd(a, f) != 1:
                return 0
            # find b = a mod f, gcd(b, c) = 1
            b = a
            while gcd(b, c) != 1:
                b += f
            return zn(ch["exponents"][b])
        ch["psi0"] = psi0
        ch["parity"] = 0 if abs(psi0(f - 1) - 1) < 1e-9 else 1
        tau = sum(psi0(a) * cmath.exp(2j * cmath.pi * a / f) for a in range(1, f + 1) if gcd(a, f) == 1)
        ch["W"] = tau / (1j ** ch["parity"] * math.sqrt(f)) if f > 1 else 1 + 0j
    return chars

def check_field(c, H, X=120, log=print, workdir=None):
    poly, cosets = period_polynomial(c, H)
    if poly is None:
        return {"c": c, "H": H, "status": "periods not distinct"}
    n = len(cosets)
    gens, idx = regular_generators(c, H, cosets)
    d = workdir or tempfile.mkdtemp()
    cfg = {"f": poly, "generators": [[i + 1 for i in g] for g in gens], "run_dir": d,
           "options": {"precision_check": False, "class_confirmation": False, "class_bound_X": X,
                       "local_unramified_check": 0, "fe_test": False, "compute_models": False, "quiet": True}}
    t0 = time.time()
    try:
        run_pipeline(cfg)
    except Exception as ex:
        return {"c": c, "H": H, "n": n, "f": poly, "status": f"pipeline failed: {ex}"}
    T = load_json(os.path.join(d, "chartable.json"))
    C = load_json(os.path.join(d, "classes.json"))
    CJ = load_json(os.path.join(d, "conductors.json"))
    EU = load_json(os.path.join(d, "euler.json"))
    RN = load_json(os.path.join(d, "rootnumbers.json"))
    AR = load_json(os.path.join(d, "archimedean.json"))
    G = load_json(os.path.join(d, "group.json"))
    r = T["r"]
    vals = [[Cyc.from_json(v).to_complex() for v in row] for row in T["values"]]
    # the isomorphism G -> (Z/c)^x/H from the Frobenius classes: class k at ell  <->  ell H
    frob = {rec["ell"]: rec["class"] - 1 for rec in C["primes"]}
    chars = cft_data(c, H, cosets, idx, X)
    # match rows of the table with Dirichlet characters by comparing chi(Frob_ell) with psi(ell) at all recorded ell
    matching, problems = {}, []
    for nu in range(r):
        found = None
        for j, ch in enumerate(chars):
            if all(abs(vals[nu][k] - cmath.exp(2j * cmath.pi * ch["exponents"][ell % c] / n)) < 1e-9 for ell, k in frob.items()):
                found = j
                break
        if found is None:
            problems.append(f"row {nu+1}: no Dirichlet character matches its Frobenius values")
        else:
            matching[nu] = found
    if len(set(matching.values())) != r:
        problems.append("matching not bijective")
    # conductors
    for nu, j in matching.items():
        pc = CJ["conductors"][nu]["partial_conductor"]
        if CJ["conductors"][nu]["unknown_at"]:
            problems.append("wild prime present")
        if pc != chars[j]["conductor"]:
            problems.append(f"chi_{nu+1}: conductor {pc} vs CFT {chars[j]['conductor']}")
        # Euler factors at ramified primes
        for ell_s, rec in EU["ramified"].items():
            ell = int(ell_s)
            P = [Cyc.from_json(x).to_complex() for x in rec["characters"][str(nu + 1)]["P"]]
            f = chars[j]["conductor"]
            want = [1 + 0j] if f % ell == 0 else [1 + 0j, -chars[j]["psi0"](ell)]
            if len(P) != len(want) or any(abs(a - b) > 1e-9 for a, b in zip(P, want)):
                problems.append(f"chi_{nu+1}: Euler factor at {ell}: {P} vs CFT {want}")
        # Gamma factor parity
        b = AR["characters"][nu]["b"]
        if b != chars[j]["parity"]:
            problems.append(f"chi_{nu+1}: parity b={b} vs psi(-1)")
        # root number
        rn = next(x for x in RN["characters"] if x["chi"] == nu + 1)
        if rn.get("W_complex") is None:
            problems.append(f"chi_{nu+1}: no W")
        else:
            W = complex(*rn["W_complex"])
            if abs(W - chars[j]["W"]) > 1e-8:
                problems.append(f"chi_{nu+1}: W {W} vs CFT {chars[j]['W']}")
    return {"c": c, "H": H, "n": n, "f": poly, "primes_checked": len(frob), "status": "ok" if not problems else "MISMATCH",
            "problems": problems, "seconds": round(time.time() - t0, 1)}

def squarefree(c):
    for p in range(2, int(c ** 0.5) + 1):
        if c % (p * p) == 0:
            return False
    return True

def fields_of_conductor(c, degrees):
    """Subgroups H of (Z/c)^x with conductor exactly c and index in degrees: enumerate subgroups
    generated by up to two elements (enough for the small quotients used), dedupe, keep those whose
    fixed field has conductor c (no character of the quotient trivial on the kernel of reduction mod c/p)."""
    U = units(c)
    subs = set()
    for a in U:
        subs.add(tuple(subgroup_generated(c, [a])))
    for a in U:
        for b in U:
            if a < b:
                subs.add(tuple(subgroup_generated(c, [a, b])))
    out = []
    for H in subs:
        n = len(U) // len(H)
        if n not in degrees:
            continue
        # conductor exactly c: H must not contain the kernel of (Z/c)^x -> (Z/(c/p))^x for any p | c
        ok = True
        for p in range(2, c + 1):
            if c % p == 0 and all(p % q for q in range(2, p)):
                ker = [a for a in U if a % (c // p) == 1 % (c // p)]
                if all(k in H for k in ker):
                    ok = False
                    break
        if ok:
            out.append(list(H))
    return out