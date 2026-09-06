"""The functional-equation test and the subfield zeta identities.
"""
from __future__ import annotations
import math, cmath
from math import gcd

from .analytic import (Kernel, closed_form, dirichlet_coefficients, functional_equation_defect,
                       tail_estimate, zeta_subfield_coefficients, euler_factor_coeffs, theta)
from .stages_euler import unramified_euler
from .cyclo import Cyc
from .chartable import is_prime
from .filtration import HardFailure
from .ramified import factor_mod, decide_prime_for_factor, discriminant
from .schur import multiplicity, Pair
from .archimedean import real_root_count

def choose_X(kernel, d, sqrt_f, t_min, eps, cap):
    X = 50
    while X < cap:
        if tail_estimate(d, kernel, sqrt_f, t_min, X, factor=3) < eps:
            return X
        X = int(X * 1.25) + 1
    return cap

def run_analytic(G, cl, table, f, assigner, EU, CJ, RN, AR, fam_resolvents, log=print,
                 eps=1e-12, cap=20000, ts=(1.1, 1.3), seed=0):
    out = {"characters": [], "subfield_identities": {}}
    wild = CJ["wild_primes"]
    conds = {c["chi"]: c["partial_conductor"] for c in CJ["conductors"]}
    ab = {d["chi"]: (d["a"], d["b"]) for d in AR["characters"]}
    Wc = {c["chi"]: complex(*c["W_complex"]) for c in RN["characters"] if c.get("W_complex") is not None}
    conj_row = table.galois_action[table.e - 1] if table.e > 1 else list(range(table.r))
    ramified_euler = {int(ell): {int(k) - 1: [Cyc.from_json(c) for c in v["P"]] for k, v in rec["characters"].items()}
                      for ell, rec in EU["ramified"].items()}
    # X for the test
    kernels = {}
    X_needed = 0
    for nu in range(1, table.r):
        if wild:
            break
        a, b = ab[nu + 1]
        key = (a, b)
        if key not in kernels:
            kernels[key] = Kernel(a, b)
        X_needed = max(X_needed, choose_X(kernels[key], a + b, math.sqrt(conds[nu + 1]), min(ts) / max(ts), eps, cap))
    X = max(X_needed, 200)
    log(f"functional equation: X = {X} (eps = {eps}, cap = {cap}); assigning classes at the extra primes")
    # Euler factors at all primes <= X
    euler = {nu: {} for nu in range(table.r)}
    classes = {}
    q = 1
    while True:
        q += 1
        if q > X:
            break
        if not is_prime(q):
            continue
        if assigner.policy.Delta % q == 0:
            if q in ramified_euler:
                for nu in range(table.r):
                    euler[nu][q] = euler_factor_coeffs(ramified_euler[q][nu])
            elif not wild:
                raise HardFailure(f"no Euler factor at {q}")
            continue
        k = assigner.assign(q)["class"] - 1
        classes[q] = k
        for nu in range(table.r):
            euler[nu][q] = euler_factor_coeffs(unramified_euler(table, nu, k, cl))
    coeffs = {nu: dirichlet_coefficients(euler[nu], X) for nu in range(table.r)}
    # ---- functional-equation defects
    for nu in range(1, table.r):
        if wild:
            out["characters"].append({"chi": nu + 1, "status": f"no test: wild primes {wild}"})
            continue
        a, b = ab[nu + 1]
        g = kernels[(a, b)]
        sf = math.sqrt(conds[nu + 1])
        W = Wc[nu + 1]
        nb = conj_row[nu]
        rec = {"chi": nu + 1, "conductor": conds[nu + 1], "gamma": {"a": a, "b": b}, "W": [W.real, W.imag], "tests": []}
        for t in ts:
            lhs, rhs, defect = functional_equation_defect(coeffs[nu], coeffs[nb], g, sf, W, t, X)
            bound = tail_estimate(a + b, g, sf, min(t, 1 / t), X) * (1 + 1 / t)
            scale = max(abs(lhs), abs(rhs), 1e-300)
            ok = defect <= bound + 1e-9 * scale
            rec["tests"].append({"t": t, "Theta": [lhs.real, lhs.imag], "W_Theta_bar": [rhs.real, rhs.imag],
                                 "defect": defect, "bound": bound, "relative": defect / scale, "ok": ok,
                                 "informative": bound < 1e-6 * scale})
            if not ok:
                raise HardFailure(f"functional-equation defect above bound for chi_{nu+1} at t = {t}: defect {defect}, bound {bound}")
        # wrong-sign / wrong-W diagnostics: the defect with -W must be large
        lhs, rhs, d2 = functional_equation_defect(coeffs[nu], coeffs[nb], g, sf, -W, ts[0], X)
        rec["defect_with_minus_W"] = d2
        rec["status"] = "ok"
        out["characters"].append(rec)
    if not wild:
        log(f"functional equation: defects {[(r['chi'], '%.2e' % r['tests'][0]['relative']) for r in out['characters'] if r.get('tests')]}; with -W: {[('%.2e' % r['defect_with_minus_W']) for r in out['characters'] if r.get('tests')]}")
    # ---- the subfield zeta identities on the family (only when every prime <= X has an Euler factor, i.e. no wild prime)
    XA = min(X, 2000)
    if wild:
        log(f"subfield identities: skipped (Euler factors unknown at the wild primes {wild})")
        fam_resolvents = []
    for Hrec in fam_resolvents:
        H, R, label = Hrec["H"], Hrec["resolvent"], Hrec["label"]
        pair = Pair("H", H.elements(), {x: 0 for x in H.elements()}, 1, label)
        ns = [multiplicity(table, nu, pair) for nu in range(table.r)]
        # character side: Dirichlet convolution of the L(chi)^{n_chi} coefficients
        c = [0j] * (XA + 1); c[1] = 1
        for nu in range(table.r):
            for _ in range(ns[nu]):
                new = [0j] * (XA + 1)
                for i in range(1, XA + 1):
                    if c[i] == 0:
                        continue
                    for j in range(1, XA // i + 1):
                        new[i * j] += c[i] * coeffs[nu][j]
                c = new
        # direct side: ideals of norm m in N^H
        rd = {}
        dR = discriminant(R)
        for q in range(2, XA + 1):
            if not is_prime(q):
                continue
            if dR % q == 0:
                rd[q] = [fd for e_, fd in decide_prime_for_factor(R, q)["residue_decomposition"]]
            else:
                rd[q] = [len(g) - 1 for g, e in factor_mod(R, q)]
        bcoef = zeta_subfield_coefficients(rd, XA)
        maxdiff = max(abs(c[m] - bcoef[m]) for m in range(1, XA + 1))
        # r_1 = sum n_chi chi(c) vs real roots of R_H
        r1 = sum(ns[nu] * (ab[nu + 1][0] - ab[nu + 1][1]) for nu in range(table.r))
        rr = real_root_count(R)
        # discriminant valuations at the tame ramified primes: sum n_chi f_ell(chi) = v_ell(disc O_{N^H})
        dv = {}
        for ell_s, prec in CJ["primes"].items():
            fl = {cc["chi"] - 1: cc["f_ell"] for cc in prec["conductor_exponents"]}
            s = sum(ns[nu] * fl[nu] for nu in range(table.r))
            ell = int(ell_s)
            w = decide_prime_for_factor(R, ell)["v_disc_O"] if dR % ell == 0 else 0
            dv[ell_s] = {"sum_n_f": s, "v_disc": w, "ok": s == w}
        ok = maxdiff < 1e-6 and r1 == rr and all(v["ok"] for v in dv.values())
        out["subfield_identities"][label] = {"index": Hrec["index"], "n_chi": ns, "X": XA, "max_coefficient_defect": maxdiff,
                             "r1_character_side": r1, "real_roots_of_R": rr, "discriminant_valuations": dv, "ok": ok}
        if not ok:
            raise HardFailure(f"subfield identity fails for H = {label}: {out['subfield_identities'][label]}")
    log(f"subfield identities: (coefficients to {XA}, r_1, discriminant valuations) hold for {len(out['subfield_identities'])} subgroups")
    out["X"] = X
    out["classes_extra"] = len(classes)
    out["_coefficients"] = {nu + 1: [[round(c.real, 12), round(c.imag, 12)] for c in coeffs[nu]] for nu in range(table.r)}
    out["_raw"] = {"coeffs": coeffs, "euler": euler, "kernels": kernels}
    return out