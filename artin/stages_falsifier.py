"""The decision rule, its perturbations, the trivial character and the anchor.
"""
from __future__ import annotations
import math

from .analytic import Kernel, theta
from .falsifier import (falsify, trivial_character_test, zeta_anchor, lambda_at_one, E_total, T_FIXED, decide)
from .filtration import HardFailure
from .schur import Pair, multiplicity
from .chartable import is_prime
from .perm import PermGroup

def run_falsifier(G, cl, table, f, factors, AN_raw, X, CJ, RN, AR, policy, local_json, log=print):
    coeffs, euler, kernels = AN_raw["coeffs"], AN_raw["euler"], AN_raw["kernels"]
    wild = CJ["wild_primes"]
    conds = {c["chi"]: c["partial_conductor"] for c in CJ["conductors"]}
    ab = {d["chi"]: (d["a"], d["b"]) for d in AR["characters"]}
    Wc = {c["chi"]: complex(*c["W_complex"]) for c in RN["characters"] if c.get("W_complex") is not None}
    conj_row = table.galois_action[table.e - 1] if table.e > 1 else list(range(table.r))
    out = {"X": X, "t_fixed": list(T_FIXED), "characters": [], "trivial": None, "anchor": None}
    # ramified primes with nontrivial inertia
    ell_ram = sorted(int(e) for e, r in local_json["ramified"].items() if r["status"] == "ok" and r["e"] > 1)
    small = [q for q in range(2, 60) if is_prime(q) and policy.Delta % q][:2]
    L1 = {}
    for nu in range(1, table.r):
        if wild:
            out["characters"].append({"chi": nu + 1, "status": f"no test: wild primes {wild}"})
            continue
        a, b = ab[nu + 1]
        d = a + b
        g = kernels[(a, b)]
        nb = conj_row[nu]
        f_chi = conds[nu + 1]
        ram_chi = [ell for ell in ell_ram if CJ["primes"][str(ell)]["conductor_exponents"][nu]["f_ell"] > 0]
        rep = falsify(coeffs[nu], coeffs[nb], euler[nu], euler[nb], g, d, a, b, f_chi, Wc[nu + 1], X, ram_chi, small, log)
        rep["chi"] = nu + 1
        rep["status"] = "ok"
        # L(1, chi) for the anchor
        Lam1 = lambda_at_one(coeffs[nu], coeffs[nb], g, f_chi, Wc[nu + 1], X)
        L1[nu] = Lam1 / (math.sqrt(f_chi) * math.pi ** (-b))
        rep["L(1,chi)"] = [L1[nu].real, L1[nu].imag]
        out["characters"].append(rep)
    if not wild:
        summary = {}
        for r in out["characters"]:
            for key in ("wrong_W", "wrong_conductor", "wrong_euler_ramified", "wrong_frobenius"):
                for e in r[key]:
                    summary.setdefault(key, {}).setdefault(e["verdict"], 0)
                    summary[key][e["verdict"]] += 1
        log(f"falsifier at t = {T_FIXED}: true data pass ({[('%.1e' % r['true'][0]['Delta'], '%.1e' % (r['true'][0]['E'] if r['true'][0]['E'] is not None else float('nan'))) for r in out['characters']][:4]}...); perturbations: {summary}")
        out["summary"] = summary
    out["trivial"] = trivial_character_test(X)
    log(f"trivial character: residuals {['%.1e' % t['residual'] for t in out['trivial']]} within E + evaluation error")
    # zeta anchor: K = Q(alpha) for the first irreducible factor
    f0 = factors[0]
    n = G.n
    roots0 = [i for i in range(n)]  # H = stabilizer of a root of f0: pick the first root index belonging to f0
    # point stabilizer of the first point (f irreducible) or of the first root of f0
    pt = None
    from .ramified import discriminant
    # find a point whose orbit has size deg f0: stabilizer of point 0 works for irreducible f; for reducible f use orbit sizes
    orbits = []
    seen = set()
    for i in range(n):
        if i in seen:
            continue
        orb = {i}; stack = [i]
        while stack:
            x = stack.pop()
            for g_ in G.generators:
                y = g_[x]
                if y not in orb:
                    orb.add(y); stack.append(y)
        seen |= orb
        orbits.append(sorted(orb))
    orb0 = next(o for o in orbits if len(o) == len(f0) - 1)
    i0 = orb0[0]
    H = [g_ for g_ in G.elements() if g_[i0] == i0]
    pair = Pair("pt", H, {x: 0 for x in H}, 1, "point stabilizer")
    ns = [multiplicity(table, nu, pair) for nu in range(table.r)]
    def table_side():
        if wild or ns[0] != 1:
            return None
        prod = 1 + 0j
        for nu in range(1, table.r):
            prod *= L1[nu] ** ns[nu]
        if abs(prod.imag) > 1e-9 * abs(prod):
            raise HardFailure("character-side product not real")
        return prod.real
    XA = min(X, 3000)
    out["anchor"] = zeta_anchor(f0, XA, table_side, log)
    out["anchor"]["H"] = {"index": len(f0) - 1, "n_chi": ns}
    log(f"zeta anchor for K = Q(alpha), index {len(f0)-1}: rho_K = {out['anchor']['zeta_residue_rho_K']:.12g}"
        + (f", prod L(1,chi)^n = {out['anchor']['character_side_product_L1']:.12g}, relative difference {out['anchor']['relative_difference']:.2e}" if 'character_side_product_L1' in out['anchor'] else " (character side unavailable)"))
    return out