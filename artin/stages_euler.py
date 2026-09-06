"""Euler factors at the ramified primes and the identity system.
"""
from __future__ import annotations
from math import gcd

from .cyclo import Cyc
from .euler import (det_character, conductor_exponent_of_class_function, euler_factor_A, euler_factor_B,
                    euler_factor_C, cpoly_mul, cpoly_pow, cpoly_eq, cpoly_json, subgroup_zeta_factor,
                    direct_zeta_factor)
from .filtration import HardFailure
from .invariants import orbit_sum, set_product
from .perm import PermGroup, mul, inverse, identity
from .precision import Policy
from .resolvent import coset_transversal, roots_at, resolvent, squarefree_certificate
from .ramified import decide_prime_for_factor, factor_mod, discriminant, valuation
from .schur import Pair, multiplicity, eigen_multiplicities
from .classes import TSCHIRNHAUS
from .invariants import Invariant

def product_identity(polys, table, fdeg, index):
    """prod_chi P_chi^{chi(1)} == (1 - T^fdeg)^index, computed over Z via Galois orbits."""
    import math
    prod = [1]
    for orb in table.orbits:
        Q = [Cyc.one(polys[orb[0]][0].e)]
        for nu in orb:
            Q = cpoly_mul(Q, polys[nu])
        Qi = []
        for c in Q:
            if not c.is_rational():
                return False
            q = c.rational()
            if q.denominator != 1:
                return False
            Qi.append(int(q))
        deg = table.degrees[orb[0]]
        for _ in range(deg):
            prod = _imul(prod, Qi)
    want = [1]
    base = [1] + [0] * (fdeg - 1) + [-1]
    for _ in range(index):
        want = _imul(want, base)
    return prod == want

def _imul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out

def linear_rows(table):
    return [nu for nu in range(table.r) if table.degrees[nu] == 1]

def row_of_values(table, vals):
    for nu in range(table.r):
        if all(table.values[nu][k] == vals[k] for k in range(table.r)):
            return nu
    return None

def kernel_subgroup(table, cl, G, nu):
    return PermGroup([g for g in G.elements() if table.values[nu][cl.class_of(g)] == 1], n=G.n)

def resolvent_for_subgroup(G, H, f, p, policy, seed=0):
    """Integer resolvent R_{G,H,F_H} (orbit sum on the base of G), squarefree by Tschirnhaus."""
    n = G.n
    T = coset_transversal(G, H)
    F0 = orbit_sum(n, H.elements(), G.base(), label=f"orbit sum H order {H.order()}")
    for Tr in TSCHIRNHAUS:
        F = Invariant(n, F0.terms, Tr, F0.label)
        trial = Policy(f, policy.Delta, mode=policy.mode)
        trial.register(F, len(T))
        ring, roots, _ = roots_at(f, p, trial.k(p), seed=seed)
        R, _ = resolvent(F, T, roots, ring)
        if squarefree_certificate(R) is not None:
            policy.register(F, len(T), label=f"subgroup resolvent {F.label}, m={len(T)}")
            return R
    raise RuntimeError("no squarefree resolvent for the subgroup")

def unramified_euler(table, nu, k_class, cl):
    """det(1 - rho(Frob) T) from the eigenvalue multiplicities at class k."""
    e = table.e
    o = cl.orders[k_class]
    P = [Cyc.one(e)]
    for j, m in enumerate(eigen_multiplicities(table, nu, k_class)):
        for _ in range(m):
            P = cpoly_mul(P, [Cyc.one(e), -Cyc.root_of_unity(e, o, j)])
    return P

def run_euler(G, cl, table, f, p, policy, objects, FILT, models, assigner, ramified_json, log=print, unram_primes=3, seed=0):
    out = {"ramified": {}, "determinant_checks": {}, "zeta_identities": {}}
    n = G.n
    # ---- subgroup family for the Euler-factor identity system (with their resolvents and residue decompositions per prime)
    lin = linear_rows(table)
    family = []
    seen = set()
    def add(H, label):
        key = frozenset(H.elements())
        if key in seen or H.order() == 1 and G.order() > 24:
            return
        seen.add(key)
        family.append((H, label))
    pt = PermGroup([g for g in G.elements() if g[0] == 0], n=n)
    add(pt, "point stabilizer")
    for k in range(1, cl.r):
        add(PermGroup([cl.reps[k]], n=n), f"cyclic <g_{k+1}>")
    for nu in lin:
        add(kernel_subgroup(table, cl, G, nu), f"ker chi_{nu+1}")
    add(G, "G")
    fam = []
    for H, label in family:
        if G.order() // H.order() > 30:
            continue
        R = f if (label == "point stabilizer" and len(ramified_json["factors"]) == 1) else resolvent_for_subgroup(G, H, f, p, policy, seed)
        pair = Pair("subgroup", H.elements(), {x: 0 for x in H.elements()}, 1, label)
        fam.append({"H": H, "label": label, "resolvent": R, "pair": pair, "index": G.order() // H.order()})
    log(f"subgroup family for the zeta identities: {[(x['label'], x['index']) for x in fam]}")

    # ---- ramified tame primes
    for ell, (rec, local, D, I, Fr) in objects.items():
        FL = FILT[str(ell)]
        fl = {c["chi"] - 1: (c["f_ell"], c["swan"]) for c in FL["conductor_exponents"]}
        N_by_class = {int(k) - 1: v for k, v in FL.get("N_by_class", {}).items()}
        delta = FL["filtration"]["delta"]
        e_in = I.order(); fdeg = rec["f"]
        Del, Iel = D.elements(), I.elements()
        primes_rec = {}
        polys = []
        prodP = None
        for nu in range(table.r):
            dimVI = multiplicity(table, nu, Pair("I", Iel, {x: 0 for x in Iel}, 1, "I"))
            P_C, mults = euler_factor_C(table, nu, Del, Iel, Fr, cl, fdeg)
            P_B = euler_factor_B(table, nu, Iel, Fr, cl, fdeg, dimVI)
            if not cpoly_eq(P_C, P_B):
                raise HardFailure(f"the trace and coset routes for the Euler factor disagree at {ell} for chi_{nu+1}")
            f_ell, swan = fl[nu]
            if len(P_C) - 1 != dimVI or dimVI != table.degrees[nu] - f_ell + swan:
                raise HardFailure(f"Euler-factor dimension check fails at {ell} for chi_{nu+1}: deg {len(P_C)-1}, dim V^I {dimVI}, chi(1)-f+sw {table.degrees[nu]-f_ell+swan}")
            # coefficients in Q(chi): invariance under Gal(Q(zeta_E)/Q(chi))
            E = P_C[0].e
            stab = set(table.stabilizers[nu])
            field_ok = all(all(c.galois(t) == c for c in P_C) for t in range(1, E + 1) if gcd(t, E) == 1 and (t % table.e) in stab)
            if not field_ok:
                raise HardFailure(f"coefficients of P_ell(chi_{nu+1}) not in Q(chi) at {ell}")
            entry = {"chi": nu + 1, "P": cpoly_json(P_C), "field_E": E, "dim_V^I": dimVI, "eigen_mults": mults,
                     "f_ell": f_ell, "swan": swan, "routes_B_C_agree": True}
            # route A for degree <= 2 characters with a model
            rho = models.get(nu)
            if rho is not None and table.degrees[nu] <= 2:
                import json, os
                # model info: multiplicity and field from the model builder are needed; recompute dimension from rho(identity)
                Mid = rho(identity(n))
                dim_model = len(Mid)
                Fm = Mid[0][0].e
                P_A, dA = euler_factor_A(rho, dim_model, Fm, Iel, Fr)
                a = dim_model // table.degrees[nu]
                EE = E * Fm // gcd(E, Fm)
                lhs = [c.embed(EE) for c in P_A]
                rhs = [c.embed(EE) for c in cpoly_pow(P_C, a)]
                agree = cpoly_eq(lhs, rhs) and dA == a * dimVI
                entry["route_A"] = {"model_multiplicity": a, "dim_e_I_image": dA, "agrees_with_C": agree}
                if not agree:
                    raise HardFailure(f"the matrix-model route disagrees with the coset route at {ell} for chi_{nu+1}")
            polys.append(P_C)
            primes_rec[nu + 1] = entry
        # product identity prod P^{chi(1)} = (1 - T^f)^{[G:D]}: Galois orbits of characters first
        # (their products have integer coefficients), then integer polynomial arithmetic
        if not product_identity(polys, table, fdeg, G.order() // D.order()):
            raise HardFailure(f"the Euler-factor product identity fails at {ell}")
        out["ramified"][str(ell)] = {"characters": primes_rec, "product_identity_ok": True,
                                     "route_A_checked": [nu + 1 for nu in range(table.r) if "route_A" in primes_rec[nu + 1]]}
        # ---- determinant characters and conductor-discriminant identities at ell
        a14 = {}
        for nu in range(table.r):
            dv = det_character(table, nu)
            row = row_of_values(table, dv)
            if row is None:
                raise HardFailure("det rho_chi is not a row of the table")
            f_det = fl[row][0]
            if f_det > fl[nu][0]:
                raise HardFailure(f"the determinant conductor bound fails at {ell}: f(det) > f(chi) for chi_{nu+1}")
            a14[nu + 1] = {"det_row": row + 1, "f_ell(det)": f_det, "f_ell(chi)": fl[nu][0]}
        # conductor-discriminant for the fields N_lambda against the discriminant of the corresponding field
        lam_checks = []
        for nu in lin:
            # rows of lambda^k, k = 1..o-1 (the trivial character contributes 0)
            vals = list(table.values[nu])
            powers, cur = [], [Cyc.one(table.e)] * cl.r
            while True:
                cur = [x * y for x, y in zip(cur, vals)]
                if all(c == 1 for c in cur):
                    break
                r2 = row_of_values(table, cur)
                if r2 is None:
                    raise HardFailure("power of a linear character is not a row")
                powers.append(r2)
            s = sum(fl[r][0] for r in powers)
            H = kernel_subgroup(table, cl, G, nu)
            entry = next((x for x in fam if x["label"] == f"ker chi_{nu+1}"), None)
            if entry is None:
                continue
            w = decide_prime_for_factor(entry["resolvent"], ell) if discriminant(entry["resolvent"]) % ell == 0 else {"v_disc_O": 0}
            if w["v_disc_O"] != s:
                raise HardFailure(f"the conductor-discriminant identity fails at {ell} for lambda = chi_{nu+1}: sum f = {s}, v(disc N_lambda) = {w['v_disc_O']}")
            lam_checks.append({"lambda": nu + 1, "powers": [r + 1 for r in powers], "sum_f": s, "v_disc_N_lambda": w["v_disc_O"]})
        out["determinant_checks"][str(ell)] = {"det": a14, "N_lambda_checks": lam_checks}
        # ---- the Euler-factor identity system at ell
        idents = []
        for x in fam:
            lhs, ns = subgroup_zeta_factor(x["pair"], table, cl, polys, table.degrees)
            R = x["resolvent"]
            if x["label"] == "point stabilizer" and len(ramified_json["factors"]) == 1:
                w = ramified_json["per_prime"].get(str(ell), {"factors": {}})["factors"].get("0")
                dec = [tuple(t) for t in w["residue_decomposition"]] if w else [(1, len(g) - 1) for g, e in factor_mod(f, ell)]
            else:
                dR = discriminant(R)
                dec = [tuple(t) for t in decide_prime_for_factor(R, ell)["residue_decomposition"]] if dR % ell == 0 else [(1, len(g) - 1) for g, e in factor_mod(R, ell)]
            rhs = direct_zeta_factor(dec, lhs[0].e)
            ok = cpoly_eq(lhs, rhs)
            idents.append({"H": x["label"], "index": x["index"], "n_chi": ns, "residue_decomposition": dec, "ok": ok})
            if not ok:
                raise HardFailure(f"zeta identity fails at {ell} for H = {x['label']}")
        out["zeta_identities"][str(ell)] = idents
        log(f"ell = {ell}: Euler factors for {table.r} characters (routes B=C; route A on {out['ramified'][str(ell)]['route_A_checked']}), product identity ok; determinant and conductor-discriminant checks ok; zeta identities ok for {len(idents)} subgroups")
    # ---- the Euler-factor identity system at a few unramified primes, with the class assignment's classes
    from .chartable import is_prime
    q, cnt = 2, 0
    unram = {}
    while cnt < unram_primes:
        q += 1
        if not is_prime(q) or policy.Delta % q == 0 or q == p:
            continue
        k_class = assigner.assign(q)["class"] - 1
        polys = [unramified_euler(table, nu, k_class, cl) for nu in range(table.r)]
        idents = []
        for x in fam:
            lhs, ns = subgroup_zeta_factor(x["pair"], table, cl, polys, table.degrees)
            dec = [(1, len(g) - 1) for g, e in factor_mod(x["resolvent"], q)] if discriminant(x["resolvent"]) % q else [tuple(t) for t in decide_prime_for_factor(x["resolvent"], q)["residue_decomposition"]]
            rhs = direct_zeta_factor(dec, lhs[0].e)
            ok = cpoly_eq(lhs, rhs)
            if not ok:
                raise HardFailure(f"zeta identity fails at unramified {q} for H = {x['label']}")
            idents.append({"H": x["label"], "ok": ok, "residue_degrees": sorted(fd for e_, fd in dec)})
        unram[str(q)] = {"class": k_class + 1, "identities": idents}
        cnt += 1
    out["zeta_identities_unramified"] = unram
    out["_family"] = [{"H": x["H"], "resolvent": x["resolvent"], "label": x["label"], "index": x["index"]} for x in fam]
    log(f"zeta identities at unramified {sorted(unram)}: ok")
    return out