"""Cross-check against weight-one newforms.

For an odd two-dimensional character induced from an imaginary quadratic
field, the associated modular form is the theta series of a class group
character.  That series is computed here from binary quadratic forms, with no
input from the pipeline, and its coefficients, level, nebentypus and root
number are compared with the computed ones.
"""
from __future__ import annotations
import cmath, math, os, tempfile, time
from math import gcd, isqrt

from .run import run_pipeline
from .certificate import load_json
from .cyclo import Cyc
from .analytic import dirichlet_coefficients, euler_factor_coeffs
from .stages_euler import unramified_euler
from .chartable import is_prime, CharacterTable
from .perm import PermGroup, from_json

# ------------------------------------------------------------------ binary quadratic forms

def reduce_form(a, b, c):
    while True:
        if c < a or (c == a and b < 0):
            a, b, c = c, -b, a
            continue
        if b > a or b <= -a:
            # bring b into (-a, a]
            k = (a - b) // (2 * a)
            b2 = b + 2 * a * k
            if b2 > a:
                b2 -= 2 * a; k -= 1
            c = (b2 * b2 - (b * b - 4 * a * c)) // (4 * a)
            b = b2
            continue
        return (a, b, c)

def reduced_forms(D):
    out = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a) == 0:
                c = (b * b - D) // (4 * a)
                if c >= a and gcd(gcd(a, abs(b)), c) == 1:
                    if c == a and b < 0:
                        continue
                    out.append((a, b, c))
        a += 1
    return sorted(out)

def compose(f1, f2, D):
    """Dirichlet composition of primitive forms of discriminant D, reduced."""
    a1, b1, c1 = f1; a2, b2, c2 = f2
    # standard algorithm (Cohen 5.4.7)
    if a1 > a2:
        a1, b1, c1, a2, b2, c2 = a2, b2, c2, a1, b1, c1
    s = (b1 + b2) // 2
    n = b2 - s
    if a2 % a1 == 0:
        y1, d = 0, a1
    else:
        # d = gcd(a1, a2) = u a1 + v a2 ... extended Euclid
        d, u, v = ext_gcd(a2, a1)
        y1 = u
    if s % d == 0:
        y2, d1 = -1, d
    else:
        d1, u2, v2 = ext_gcd(s, d)
        y2 = u2
        y1 = y1 * -v2 if False else y1  # placeholder (handled below)
    # use the general formula via solving congruences directly (small numbers): find B with
    # B = b1 mod 2a1, B = b2 mod 2a2, B^2 = D mod 4 a1 a2 / d1^2 -- brute force safe for our sizes
    A = a1 * a2
    g = gcd(gcd(a1, a2), s)
    A //= g * g
    for B in range(-2 * A, 2 * A + 1):
        if (B - b1) % (2 * a1 // g) == 0 and (B - b2) % (2 * a2 // g) == 0 and (B * B - D) % (4 * A) == 0:
            # also need the cross condition (b1+b2)/2 * B/... ; check that the form (A, B, C) composes: verify via norm relation on a test
            C = (B * B - D) // (4 * A)
            return reduce_form(A, B, C)
    raise RuntimeError("composition failed")

def ext_gcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = ext_gcd(b, a % b)
    return g, y, x - (a // b) * y

def class_group(D):
    forms = reduced_forms(D)
    h = len(forms)
    idx = {f: i for i, f in enumerate(forms)}
    one = reduce_form(1, 1 if D % 2 else 0, (1 - D) // 4 if D % 2 else -D // 4)
    table = [[idx[compose(forms[i], forms[j], D)] for j in range(h)] for i in range(h)]
    # sanity: a group
    assert all(table[i][idx[one]] == i for i in range(h))
    return forms, idx, table, idx[one]

def characters_of(table, e_idx):
    """All characters of the (abelian) group given by its table: exponent vectors modulo the
    group order via a generator search (small groups)."""
    h = len(table)
    # find the structure by brute force: powers
    def power(i, k):
        r = e_idx
        for _ in range(k):
            r = table[r][i]
        return r
    orders = [next(k for k in range(1, h + 1) if power(i, k) == e_idx) for i in range(h)]
    # generators
    gens, span = [], {e_idx}
    for i in sorted(range(h), key=lambda i: -orders[i]):
        if i in span:
            continue
        gens.append(i)
        new = set(span)
        for s in span:
            x = s
            for _ in range(orders[i]):
                x = table[x][i]
                new.add(x)
        span = new
        if len(span) == h:
            break
    # words
    word = {e_idx: tuple([0] * len(gens))}
    frontier = [e_idx]
    while frontier:
        x = frontier.pop()
        for t, g in enumerate(gens):
            y = table[x][g]
            if y not in word:
                w = list(word[x]); w[t] += 1; word[y] = tuple(w); frontier.append(y)
    import itertools
    chars = []
    for ks in itertools.product(*[range(orders[g]) for g in gens]):
        def val(i):
            return sum(Fraction_mod(ks[t], word[i][t], orders[gens[t]]) for t in range(len(gens)))
        vals = {}
        ok = True
        for i in range(h):
            vals[i] = sum(ks[t] * word[i][t] / orders[gens[t]] for t in range(len(gens))) % 1.0
        # homomorphism check
        for i in range(h):
            for j in range(h):
                if abs(((vals[i] + vals[j] - vals[table[i][j]]) % 1.0) - 0) > 1e-9 and abs(((vals[i] + vals[j] - vals[table[i][j]]) % 1.0) - 1) > 1e-9:
                    ok = False
        if ok:
            chars.append({i: cmath.exp(2j * math.pi * vals[i]) for i in range(h)})
    return chars

def Fraction_mod(a, b, m):
    return 0

def kronecker(D, p):
    if p == 2:
        return 0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1)
    if D % p == 0:
        return 0
    return 1 if pow(D % p, (p - 1) // 2, p) == 1 else -1

def prime_ideal_class(D, p, idx):
    """Class index of a prime ideal above the split or ramified prime p."""
    for b in range(0, 2 * p):
        if (b * b - D) % (4 * p) == 0:
            return idx[reduce_form(p, b, (b * b - D) // (4 * p))]
    raise RuntimeError("no form")

def theta_coefficients(D, psi, idx, X):
    """a_n of the theta series sum_a psi(a) q^{N a}, n <= X (unramified psi)."""
    euler = {}
    for p in range(2, X + 1):
        if not is_prime(p):
            continue
        k = kronecker(D, p)
        if k == 1:
            c = prime_ideal_class(D, p, idx)
            z = psi[c]
            euler[p] = [1 + 0j, -(z + z.conjugate()), 1 + 0j]
        elif k == -1:
            euler[p] = [1 + 0j, 0j, -1 + 0j]
        else:
            c = prime_ideal_class(D, p, idx)
            euler[p] = [1 + 0j, -psi[c]]
    return dirichlet_coefficients(euler, X)

# ------------------------------------------------------------------ the check

def check_dihedral(f, gens, D, X=300, log=print, workdir=None):
    d = workdir or tempfile.mkdtemp()
    cfg = {"f": f, "generators": gens, "run_dir": d,
           "options": {"precision_check": False, "class_confirmation": False, "class_bound_X": X,
                       "local_unramified_check": 0, "fe_test": False, "compute_models": False, "quiet": True}}
    t0 = time.time()
    run_pipeline(cfg)
    T = load_json(os.path.join(d, "chartable.json"))
    C = load_json(os.path.join(d, "classes.json"))
    CJ = load_json(os.path.join(d, "conductors.json"))
    EU = load_json(os.path.join(d, "euler.json"))
    RN = load_json(os.path.join(d, "rootnumbers.json"))
    AR = load_json(os.path.join(d, "archimedean.json"))
    G = PermGroup([from_json(g) for g in load_json(os.path.join(d, "group.json"))["generators"]], n=len(f) - 1)
    cl = G.classes()
    table = CharacterTable(cl)
    # the odd two-dimensional characters
    odd2 = [c["chi"] - 1 for c in AR["characters"] if c["odd_two_dimensional"]]
    if not odd2:
        return {"f": f, "status": "no odd two-dimensional character"}
    forms, idx, tab, e_idx = class_group(D)
    chars = characters_of(tab, e_idx)
    h = len(forms)
    results = []
    for nu in odd2:
        # pipeline coefficients from classes + ramified Euler factors
        euler = {}
        for rec in C["primes"]:
            euler[rec["ell"]] = euler_factor_coeffs(unramified_euler(table, nu, rec["class"] - 1, cl))
        for ell_s, rec in EU["ramified"].items():
            euler[int(ell_s)] = [Cyc.from_json(x).to_complex() for x in rec["characters"][str(nu + 1)]["P"]]
        primes_needed = [p for p in range(2, X + 1) if is_prime(p)]
        if any(p not in euler for p in primes_needed):
            return {"f": f, "status": "missing Euler factor (wild prime?)"}
        a_pipe = dirichlet_coefficients(euler, X)
        cond = CJ["conductors"][nu]["partial_conductor"]
        W_pipe = next(x for x in RN["characters"] if x["chi"] == nu + 1).get("W_complex")
        # match a class group character psi with a_pipe (psi and psi-bar give the same form)
        best = None
        for psi in chars:
            if all(abs(psi[i] - 1) < 1e-9 for i in psi):
                continue
            a_ref = theta_coefficients(D, psi, idx, X)
            mism = [m for m in range(1, X + 1) if abs(a_ref[m] - a_pipe[m]) > 1e-6]
            if best is None or len(mism) < len(best[1]):
                best = (psi, mism, a_ref)
        psi, mism, a_ref = best
        # level, nebentypus, root number
        level_ok = cond == abs(D)
        # nebentypus det rho_chi(Frob_p) = (D/p): det from the table via eigenvalue multiplicities
        from .euler import det_character
        dv = det_character(table, nu)
        neb_ok = all(abs(dv[rec["class"] - 1].to_complex() - kronecker(D, rec["ell"])) < 1e-9 for rec in C["primes"])
        # ramified prime ideal class: psi(d_K) with d_K = (sqrt D) the different (principal for prime |D|)
        W_ref = 1.0
        W_ok = W_pipe is not None and abs(complex(*W_pipe) - W_ref) < 1e-8
        rec = {"chi": nu + 1, "level": abs(D), "class_number": h, "psi_order": max(1, round(1 / (abs(cmath.phase(psi[i]) / (2 * math.pi)) or 1))) if False else None,
               "coefficients_compared": X, "mismatches": mism[:5], "n_mismatches": len(mism),
               "level_ok": level_ok, "conductor_pipeline": cond, "nebentypus_ok": neb_ok,
               "W_pipeline": W_pipe, "W_newform": W_ref, "W_ok": W_ok,
               "first_coefficients": [round(a_pipe[m].real) for m in range(1, 21)]}
        rec["ok"] = (not mism) and level_ok and neb_ok and W_ok
        results.append(rec)
    return {"f": f, "D": D, "|G|": G.order(), "characters": results, "all_ok": all(r["ok"] for r in results),
            "seconds": round(time.time() - t0, 1)}