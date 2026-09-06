"""Re-checking a run directory in place.

Recomputes the checkable parts of a run from its own records: the character
table certificate, the ramified-prime witnesses, the precision consultations,
the Frobenius classes against the resolvent factorisations, the local data,
the conductor exponents, the Euler factors, the archimedean data, the root
numbers and the functional-equation defects.  For a check that shares nothing
with this package, use artinverify instead.
"""
import os
import sys

from .certificate import load_json, DescentCertificate
from .chartable import verify_certificate
from .ramified import verify_ramified
from .precision import mu, ilog

def verify_run(run_dir):
    report = {}
    D = DescentCertificate(load_json(os.path.join(run_dir, "descent.json")), source=run_dir)
    cl = D.G.classes()
    tab = load_json(os.path.join(run_dir, "chartable.json"))
    ok_t, msg_t = verify_certificate(cl, tab["certificate"])
    report["character_table"] = {"ok": ok_t, "message": msg_t, "S": tab["certificate"]["S"]}
    R = load_json(os.path.join(run_dir, "ramified.json"))
    ok_r, msgs = verify_ramified(R)
    report["ramified"] = {"ok": ok_r, "messages": msgs}
    # precision policy: recompute M* from the recorded pairs and check every consultation
    pj = load_json(os.path.join(run_dir, "precision.json"))
    # M* as a function of the number of pairs registered so far (registration only raises it)
    prefix = [1]
    for pr in pj["pairs"]:
        B = pr["norm1"] * pj["root_bound_R"] ** pr["degree"]
        prefix.append(max(prefix[-1], (2 * B + 2) ** (mu(pr["m"]) if pj["mode"] == "conservative" else pr["m"])))
    ok_p = (prefix[-1].bit_length() - 1) == pj["log2_Mstar"]
    bad = []
    Delta = abs(int(load_json(os.path.join(run_dir, "ramified.json"))["disc_f"]))
    for entry in pj["consultation_log"]:
        ell, k = entry["ell"], entry["k"]
        Mstar = prefix[entry["pairs"]]
        v = 0
        d = Delta
        while d % ell == 0:
            d //= ell
            v += 1
        need = max(ilog(ell, Mstar) + 1, v + 1)
        if k < need:
            ok_p = False
            bad.append(f"ell={ell}: recorded {k} < required {need} with {entry['pairs']} pairs")
    report["precision"] = {"ok": ok_p, "problems": bad, "mode": pj["mode"]}
    # class assignment (the verifier specification V6(a),(b) on the recorded data)
    cpath = os.path.join(run_dir, "classes.json")
    if os.path.exists(cpath):
        from .classes import candidate_subgroups, coset_action, signature
        from .ramified import factor_mod
        from .resolvent import divides_disc
        from .perm import cycle_type
        C = load_json(cpath)
        subs = {str(S.label): S for S in candidate_subgroups(D.G, cl)}
        fam = []
        for t in C["family"]:
            S = subs[str(t["subgroup"])]
            T, act = coset_action(D.G, S.H)
            fam.append((t, {k: signature(cl, act, k) for k in range(cl.r)}))
        ok_c, bad_c = True, []
        for rec in C["primes"]:
            ell, k = rec["ell"], rec["class"] - 1
            facs = factor_mod(D.f, ell)
            ct = sorted((len(g) - 1 for g, e in facs for _ in range(e)), reverse=True)
            if ct != rec["block"] or list(cl.cycle_types[k]) != ct:
                ok_c = False; bad_c.append(f"ell={ell}: block/cycle type mismatch")
            for t, sig in fam:
                if divides_disc(t["resolvent"], ell):
                    continue
                fR = factor_mod(t["resolvent"], ell)
                ctR = tuple(sorted((len(g) - 1 for g, e in fR for _ in range(e)), reverse=True))
                if sig[k] != ctR:
                    ok_c = False; bad_c.append(f"ell={ell}: class {k+1} inconsistent with resolvent {t['subgroup']}")
        report["classes"] = {"ok": ok_c, "problems": bad_c[:10], "primes": len(C["primes"])}
    else:
        ok_c = True
        report["classes"] = {"ok": True, "problems": [], "primes": 0}
    # local data (the verifier specification V2(b)-(c), V3(e) Check 2 on the recorded data)
    lpath = os.path.join(run_dir, "local.json")
    ok_l, bad_l = True, []
    if os.path.exists(lpath):
        from .perm import PermGroup, from_json, mul, inverse
        L = load_json(lpath)
        for ell, r in L["ramified"].items():
            if r["status"] != "ok":
                continue
            Dg = [from_json(g) for g in r["D_generators"]]
            Ig = [from_json(g) for g in r["I_generators"]]
            if not all(D.G.contains(g) for g in Dg):
                ok_l = False; bad_l.append(f"ell={ell}: D not in G")
            Dl, Il = PermGroup(Dg, n=D.n), PermGroup(Ig, n=D.n)
            if Dl.order() != r["|D|"] or Il.order() != r["e"] or Dl.order() != r["e"] * r["f"]:
                ok_l = False; bad_l.append(f"ell={ell}: orders")
            for g in Dg:
                for h in Ig:
                    if not Il.contains(mul(g, mul(h, inverse(g)))):
                        ok_l = False; bad_l.append(f"ell={ell}: I not normal in D")
            if not r["inertia_orbits_match_factors"] or not r["decomposition_orbits_match_factors"]:
                ok_l = False; bad_l.append(f"ell={ell}: orbit sizes vs the factor data")
            if not D.G.contains(mul(inverse(from_json(r["matching"])), from_json(r["matching2"]))):
                ok_l = False; bad_l.append(f"ell={ell}: matchings not in one G-coset")
    report["local"] = {"ok": ok_l, "problems": bad_l}
    # conductors (the verifier specification V3(c),(e), V4(b),(c),(d)) recomputed from the recorded i_L and the table
    cpath2 = os.path.join(run_dir, "conductors.json")
    ok_f, bad_f = True, []
    if os.path.exists(cpath2) and os.path.exists(lpath):
        from fractions import Fraction
        from .cyclo import Cyc
        from .chartable import CharacterTable
        T = CharacterTable(cl)
        CJ = load_json(cpath2)
        for ell_s, rec in CJ["primes"].items():
            lrec = L["ramified"][ell_s]
            rho = from_json(lrec["matching"])
            fil = rec["filtration"]
            i_L = {tuple(int(x) - 1 for x in k.strip("[]").split(",")): v for k, v in fil["i_L"].items()} if fil["i_L"] else {}
            e = lrec["e"]
            delta = sum(i_L.values())
            if delta != fil["delta"]:
                ok_f = False; bad_f.append(f"ell={ell_s}: delta")
            # P must be the ell-Sylow: elements with i_L >= 2 form a group of ell-power order, index prime to ell
            ell = int(ell_s)
            P = [s for s, v in i_L.items() if v >= 2]
            m = len(P) + 1
            while m % ell == 0:
                m //= ell
            if m != 1 or (e // (len(P) + 1)) % ell == 0:
                ok_f = False; bad_f.append(f"ell={ell_s}: D_1 not the ell-Sylow")
            # Hasse-Arf from the recorded orders (abelian case only checkable via the group; use orders on the whole I when abelian)
            N = {}
            for s, v in i_L.items():
                k = cl.class_of(mul(inverse(rho), mul(s, rho)))
                N[k] = N.get(k, 0) + v
            tot = 0
            for c in rec["conductor_exponents"]:
                nu = c["chi"] - 1
                val = T.values[nu][0] * delta
                for k, Nk in N.items():
                    val = val - T.values[nu][k] * Nk
                q = val.rational() / e
                if q.denominator != 1 or int(q) != c["f_ell"]:
                    ok_f = False; bad_f.append(f"ell={ell_s}: f_ell(chi_{nu+1}) recomputed {q} != {c['f_ell']}")
                tot += T.degrees[nu] * c["f_ell"]
            vdN = (D.G.order() // (e * lrec["f"])) * lrec["f"] * delta
            if tot != vdN or rec["identities"]["v_ell_d_N"] != vdN:
                ok_f = False; bad_f.append(f"ell={ell_s}: conductor identity")
            for fa in rec["identities"]["factors"]:
                if fa["v_ell(disc O)"] is not None and fa["f_ell(pi)"] != fa["v_ell(disc O)"]:
                    ok_f = False; bad_f.append(f"ell={ell_s}: factor conductor identity, factor {fa['factor']}")
    report["conductors"] = {"ok": ok_f, "problems": bad_f}
    # Euler factors (the verifier specification V5(a)-(e)) recomputed by route (C) from the recorded D, I, Frobenius
    epath = os.path.join(run_dir, "euler.json")
    ok_e, bad_e = True, []
    if os.path.exists(epath) and os.path.exists(lpath):
        from .euler import euler_factor_C, cpoly_eq, cpoly_pow, cpoly_mul
        from .schur import multiplicity, Pair
        from .cyclo import Cyc
        from .chartable import CharacterTable
        T = CharacterTable(cl)
        EJ = load_json(epath)
        for ell_s, rec in EJ["ramified"].items():
            lrec = L["ramified"][ell_s]
            Dl = PermGroup([from_json(g) for g in lrec["D_generators"]], n=D.n)
            Il = PermGroup([from_json(g) for g in lrec["I_generators"]], n=D.n)
            Fr = from_json(lrec["frobenius_rep"])
            Del, Iel = Dl.elements(), Il.elements()
            Ps_ = {}
            for chi_s, c in rec["characters"].items():
                nu = int(chi_s) - 1
                P, mults = euler_factor_C(T, nu, Del, Iel, Fr, cl, lrec["f"])
                Prec = [Cyc.from_json(x) for x in c["P"]]
                if not cpoly_eq(P, Prec):
                    ok_e = False; bad_e.append(f"ell={ell_s}: P(chi_{nu+1}) differs")
                dimVI = multiplicity(T, nu, Pair("I", Iel, {x: 0 for x in Iel}, 1, "I"))
                if len(P) - 1 != dimVI or dimVI != T.degrees[nu] - c["f_ell"] + c["swan"]:
                    ok_e = False; bad_e.append(f"ell={ell_s}: dimension check chi_{nu+1}")
                Ps_[nu] = P
            from .stages_euler import product_identity
            if not product_identity([Ps_[nu] for nu in range(T.r)], T, lrec["f"], D.G.order() // Dl.order()):
                ok_e = False; bad_e.append(f"ell={ell_s}: product identity")
        for ell_s, rec in EJ["determinant_checks"].items():
            for chi_s, c in rec["det"].items():
                if c["f_ell(det)"] > c["f_ell(chi)"]:
                    ok_e = False; bad_e.append(f"ell={ell_s}: f(det) > f(chi) for chi_{chi_s}")
        for ell_s, idents in EJ["zeta_identities"].items():
            if not all(x["ok"] for x in idents):
                ok_e = False; bad_e.append(f"ell={ell_s}: recorded zeta identity failure")
    report["euler"] = {"ok": ok_e, "problems": bad_e}
    apath = os.path.join(run_dir, "archimedean.json")
    ok_a, bad_a = True, []
    if os.path.exists(apath):
        from .archimedean import real_root_count
        from .euler import det_character
        from .chartable import CharacterTable
        T = CharacterTable(cl)
        AJ = load_json(apath)
        r = real_root_count(D.f)
        k_c = AJ["class_of_c"] - 1
        ct = tuple(sorted([2] * ((D.n - r) // 2) + [1] * r, reverse=True))
        if r != AJ["real_roots"] or cl.cycle_types[k_c] != ct:
            ok_a = False; bad_a.append("real root count / cycle type of c")
        for d in AJ["characters"]:
            nu = d["chi"] - 1
            chic = T.values[nu][k_c].rational()
            if d["a"] + d["b"] != T.degrees[nu] or d["a"] - d["b"] != chic:
                ok_a = False; bad_a.append(f"chi_{nu+1}: a, b")
            if det_character(T, nu)[k_c].rational() != (-1) ** d["b"]:
                ok_a = False; bad_a.append(f"chi_{nu+1}: parity")
    report["archimedean"] = {"ok": ok_a, "problems": bad_a}
    rpath = os.path.join(run_dir, "rootnumbers.json")
    ok_w, bad_w = True, []
    if os.path.exists(rpath):
        from .cyclo import Cyc
        from .rootnumber import sqrt_prime
        from .chartable import CharacterTable
        T = CharacterTable(cl)
        RJ = load_json(rpath)
        for c in RJ["characters"]:
            if c.get("W_complex") is None:
                continue
            nu = c["chi"] - 1
            k4 = (-c["b"]) % 4
            mu4 = True
            for ell_s, loc in c["local"].items():
                val = complex(*loc["complex"])
                if abs(abs(val) - 1) > 1e-7:
                    ok_w = False; bad_w.append(f"chi_{nu+1} at {ell_s}: modulus")
                if loc["i_exponent"] is None:
                    mu4 = False
                else:
                    k4 = (k4 + loc["i_exponent"]) % 4
                    if abs(val - 1j ** loc["i_exponent"]) > 1e-6:
                        ok_w = False; bad_w.append(f"chi_{nu+1} at {ell_s}: recorded i-exponent")
            real, fs = T.is_real[nu], T.indicator[nu]
            if real and fs == 1 and (not mu4 or k4 != 0):
                ok_w = False; bad_w.append(f"chi_{nu+1}: orthogonal W != 1")
            if real and fs == -1 and (not mu4 or k4 % 2):
                ok_w = False; bad_w.append(f"chi_{nu+1}: symplectic W not +-1")
            if abs(abs(complex(*c["W_complex"])) - 1) > 1e-6:
                ok_w = False; bad_w.append(f"chi_{nu+1}: |W| != 1")
    report["rootnumbers"] = {"ok": ok_w, "problems": bad_w}
    # functional-equation defects recomputed from the recorded coefficients, own kernel
    npath = os.path.join(run_dir, "analytic.json")
    ok_n, bad_n = True, []
    if os.path.exists(npath):
        import math as _m
        from .analytic import Kernel, functional_equation_defect, tail_estimate
        from .chartable import CharacterTable
        T = CharacterTable(cl)
        AN = load_json(npath)
        CO = load_json(os.path.join(run_dir, "coefficients.json"))
        X = CO["X"]
        conj_row = T.galois_action[T.e - 1] if T.e > 1 else list(range(T.r))
        kern = {}
        for c in AN["characters"]:
            if c.get("status") != "ok":
                continue
            nu = c["chi"] - 1
            a, b = c["gamma"]["a"], c["gamma"]["b"]
            if (a, b) not in kern:
                kern[(a, b)] = Kernel(a, b)
            g = kern[(a, b)]
            co = [complex(x, y) for x, y in CO["coefficients"][str(nu + 1)]]
            cob = [complex(x, y) for x, y in CO["coefficients"][str(conj_row[nu] + 1)]]
            W = complex(*c["W"])
            for tst in c["tests"]:
                lhs, rhs, dft = functional_equation_defect(co, cob, g, _m.sqrt(c["conductor"]), W, tst["t"], X)
                bnd = tail_estimate(a + b, g, _m.sqrt(c["conductor"]), min(tst["t"], 1 / tst["t"]), X) * (1 + 1 / tst["t"])
                if dft > bnd + 1e-9 * max(abs(lhs), abs(rhs), 1e-300):
                    ok_n = False; bad_n.append(f"chi_{nu+1} t={tst['t']}: defect {dft:.2e} > bound {bnd:.2e}")
        for label, v in AN["subfield_identities"].items():
            if not v["ok"]:
                ok_n = False; bad_n.append(f"subfield identity {label}")
    report["analytic"] = {"ok": ok_n, "problems": bad_n}
    report["descent_verified"] = D.verified
    report["ok"] = ok_t and ok_r and ok_p and ok_c and ok_l and ok_f and ok_e and ok_a and ok_w and ok_n
    return report

def main(argv=None):
    run_dir = (argv or sys.argv[1:])[0]
    rep = verify_run(run_dir)
    print(f"character table certificate: {'OK' if rep['character_table']['ok'] else 'FAIL'} ({rep['character_table']['message']})")
    for m in rep["ramified"]["messages"]:
        print("  " + m)
    print(f"ramified-prime witnesses: {'OK' if rep['ramified']['ok'] else 'FAIL'}")
    print(f"class records: {'OK' if rep['classes']['ok'] else 'FAIL'} ({rep['classes']['primes']} primes) {rep['classes']['problems']}")
    print(f"local data: {'OK' if rep['local']['ok'] else 'FAIL'} {rep['local']['problems']}")
    print(f"conductors: {'OK' if rep['conductors']['ok'] else 'FAIL'} {rep['conductors']['problems']}")
    print(f"euler factors: {'OK' if rep['euler']['ok'] else 'FAIL'} {rep['euler']['problems']}")
    print(f"archimedean: {'OK' if rep['archimedean']['ok'] else 'FAIL'} {rep['archimedean']['problems']}")
    print(f"root numbers: {'OK' if rep['rootnumbers']['ok'] else 'FAIL'} {rep['rootnumbers']['problems']}")
    print(f"functional equations: {'OK' if rep['analytic']['ok'] else 'FAIL'} {rep['analytic']['problems']}")
    print(f"precision policy ({rep['precision']['mode']}): {'OK' if rep['precision']['ok'] else 'FAIL'} {rep['precision']['problems']}")
    print(f"descent certificate verified by the descent artifact: {rep['descent_verified']}")
    print("VERDICT:", "ACCEPT" if rep["ok"] else "REJECT")
    return 0 if rep["ok"] else 1

if __name__ == "__main__":
    sys.exit(main())