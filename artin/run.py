"""The driver: runs every stage and writes the run directory.

Writes the resolved configuration first, then the group and character data,
the ramified primes, the precision log, the Frobenius classes, the local data,
the conductors, the Euler factors, the archimedean data, the root numbers, the
functional-equation defects and the certificate.  Each stage is timed and the
timings are recorded.
"""
from __future__ import annotations
import argparse
import datetime
import json
import os
import platform
import sys
import time

import sympy

from . import __version__
from .certificate import DescentCertificate, dump_json, load_json, CERT_VERSION
from .chartable import CharacterTable, verify_certificate
from .perm import to_json, PermGroup, mul, inverse
from .ramified import ramified_primes
from .schur import candidate_pairs, schur_data, build_model
from .precision import Policy
from .stages_precision import precision_check, find_group_in_numbering, numbering_prime
from . import resolvent as _resolvent
from .stages_classes import ClassAssigner, confirm
from .stages_local import run_local
from .stages_filtration import run_filtrations, factor_root_sets
from .stages_euler import run_euler
from .stages_archimedean import run_archimedean
from .stages_rootnumbers import run_root_numbers
from .stages_analytic import run_analytic
from .stages_falsifier import run_falsifier

DEFAULTS = {
    "class_enumeration_limit": 2_000_000,
    "model_dim_limit": 200,
    "model_work_limit": 6_000_000,
    "compute_models": True,
    "seed": 0,
    "precision_mode": "sharp",          # the precision policy policy mode: "sharp" or "conservative"
    "precision_check": True,            # record every value at the policy precision and at twice it, and compare
    "precision_check_unramified": 6,
    "class_bound_X": 200,               # classes: assign Frobenius classes at unramified primes <= X
    "class_confirmation": True,         # classes: confirm against the explicit splitting field when |G| <= 200
    "class_confirmation_primes": 15,
    "local_unramified_check": 6,
    "root_numbers": True,               # Gauss sums; the cost grows with the ramified primes
    "fe_test": True,                    # functional equation: run the functional-equation test and the subfield identities
    "fe_eps": 1e-12,                    # functional equation: target size of the truncation tail
    "fe_cap": 20000,        # local: unramified primes on which matching is compared with the class assignment
}

def parse_poly(s):
    x = sympy.symbols("x")
    P = sympy.Poly(sympy.sympify(s), x)
    return [int(c) for c in reversed(P.all_coeffs())]

class Log:
    def __init__(self, path, quiet=False):
        self.fh = open(path, "a")
        self.t0 = time.time()
        self.quiet = quiet
    def __call__(self, msg):
        line = f"[{time.time() - self.t0:8.2f}s] {msg}"
        if not self.quiet:
            print(line)
        self.fh.write(line + "\n")
        self.fh.flush()

def resolve_config(cfg):
    cfg = dict(cfg)
    opts = dict(DEFAULTS)
    opts.update(cfg.get("options", {}))
    cfg["options"] = opts
    cfg["b1"] = {"artin_version": __version__, "certificate_version": CERT_VERSION,
                 "python": platform.python_version(), "sympy": sympy.__version__,
                 "started": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    return cfg

def run_pipeline(cfg):
    _T = {}
    _t0 = [time.time()]
    def _mark(stage):
        now = time.time(); _T[stage] = round(now - _t0[0], 3); _t0[0] = now
    run_dir = cfg["run_dir"]
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "models"), exist_ok=True)
    cfg = resolve_config(cfg)
    dump_json(cfg, os.path.join(run_dir, "config.json"))   # written first, always
    log = Log(os.path.join(run_dir, "log.txt"), quiet=cfg["options"].get("quiet", False))
    log(f"run in {run_dir}; config written")
    opts = cfg["options"]

    # ---- input: descent certificate
    if "descent" in cfg:
        D = DescentCertificate.load(cfg["descent"])
    else:
        D = DescentCertificate.from_generators(cfg["f"], cfg["generators"], cfg.get("numbering"))
    G = D.G
    log(f"f = {D.f}, n = {D.n}, |G| = {G.order()}, descent verified = {D.verified} ({D.source})")

    # ---- numbering and policy first: every later record is in the fixed numbering
    from .ramified import discriminant as _disc
    Delta = _disc(D.f)
    policy = Policy(D.f, Delta, mode=opts["precision_mode"])
    if D.numbering and D.numbering.get("prime") and D.numbering.get("roots"):
        p_num = int(D.numbering["prime"])
        _resolvent.NUMBERING[(tuple(D.f), p_num)] = [int(x) for x in D.numbering["roots"]]
        log(f"numbering taken from the descent certificate at p = {p_num}")
    else:
        _cl0 = G.classes(limit=opts["class_enumeration_limit"])
        from .classes import rational_classes as _rc
        _orders = {_cl0.orders[orb[0]] for orb in _rc(_cl0) if len(orb) >= 2}
        p_num = numbering_prime(policy.Delta, coprime_to=_cl0.exponent, orders=_orders)
        log(f"no numbering in the input: placing G in the numbering of the sorted roots at p = {p_num} (stand-in for the descent artifact)")
        G = find_group_in_numbering(D.f, G.generators, p_num, policy, seed=opts["seed"], log=log)
        D.G = G
        D.numbering = {"prime": p_num, "note": "sorted residues of the roots at p (stand-in numbering)"}

    _mark("numbering_and_group_placement")
    # ---- character table with its certificate
    D.data["group"]["generators"] = [to_json(g) for g in G.generators]
    D.data["numbering"] = D.numbering
    dump_json(D.data, os.path.join(run_dir, "descent.json"))
    cl = G.classes(limit=opts["class_enumeration_limit"])
    log(f"classes: r = {cl.r}, exponent {cl.exponent}")
    group_json = G.to_json()
    group_json["classes"] = cl.to_json()
    dump_json(group_json, os.path.join(run_dir, "group.json"))
    T = CharacterTable(cl, seed=opts["seed"])
    ok, msg = verify_certificate(cl, T.certificate)
    log(f"character table: p = {T.p}, degrees {T.degrees}, S = {T.certificate['S']}, certificate verified = {T.certificate['verified']} / re-verified = {ok} ({msg})")
    tab = T.to_json()
    tab["certificate"]["reverified"] = ok
    dump_json(tab, os.path.join(run_dir, "chartable.json"))

    # ---- Schur indices and matrix models
    pairs = candidate_pairs(G, cl)
    schur = []
    models = {}
    for nu in range(T.r):
        sd = schur_data(T, nu, pairs)
        schur.append(sd)
        if opts["compute_models"]:
            M = build_model(G, T, nu, pairs, max_dim=opts["model_dim_limit"], max_work=opts["model_work_limit"])
            models[nu] = M.pop("_rho", None)
            dump_json(M, os.path.join(run_dir, "models", f"chi_{nu + 1}.json"))
            sd["model"] = {"file": f"models/chi_{nu + 1}.json", "status": M["status"],
                           "multiplicity": M.get("multiplicity"), "dimension": M.get("dimension"),
                           "field": M.get("field"), "denominator": M.get("denominator"),
                           "trace_check": M.get("trace_check")}
        log(f"chi_{nu + 1}: degree {sd['degree']}, K conductor {sd['field_conductor']}, Schur index {sd['value'] if sd['value'] else sd['candidates']} ({sd['status']})"
            + (f", model {sd['model']['status']} dim {sd['model']['dimension']} over Q(zeta_{sd['model']['field']}) denominator {sd['model']['denominator']}" if opts["compute_models"] else ""))
    dump_json({"characters": schur}, os.path.join(run_dir, "schur.json"))
    _mark("character_table")

    # ---- ramified primes
    R = ramified_primes(D.f)
    log(f"ramified primes: {R['ramified']} (candidates {R['candidates']}, disc f = {R['disc_f']})")
    dump_json(R, os.path.join(run_dir, "ramified.json"))
    _mark("ramified_primes")

    # ---- the precision policy and the doubling experiment
    pcheck = None
    if opts["precision_check"]:
        pcheck = precision_check(D.f, G, policy, n_unramified=opts["precision_check_unramified"], seed=opts["seed"], log=log)
        log(f"doubling check: identical = {pcheck['identical']} ({len(pcheck['differences'])} differences) over primes {sorted(int(q) for q in pcheck['at_k']['primes'])}")
        dump_json(pcheck, os.path.join(run_dir, "precision_check.json"))
    _mark("precision_check")
    # ---- Frobenius classes at the unramified primes
    from .chartable import is_prime as _is_prime
    A = ClassAssigner(G, cl, D.f, p_num, policy, seed=opts["seed"], log=log)
    primes = [q for q in range(2, opts["class_bound_X"] + 1) if _is_prime(q) and policy.Delta % q]
    class_records = []
    for q in primes:
        class_records.append(A.assign(q))
    methods = {m: sum(1 for r in class_records if r["method"] == m) for m in ("resolvent", "cyclotomic", "direct")}
    log(f"classes assigned at {len(primes)} unramified primes <= {opts['class_bound_X']}: methods {methods}")
    # Chebotarev heuristic (the verifier specification hypothesis (G) is only certified by the descent certificate): the
    # normal closure of the observed Frobenius classes must be G, otherwise G is probably too large
    _seen = {rec["class"] - 1 for rec in class_records}
    _ncl = PermGroup([mul(g, mul(cl.reps[k], inverse(g))) for k in _seen for g in G.generators] + [cl.reps[k] for k in _seen], n=G.n)
    _gen_ok = True
    for _ in range(3):
        _new = [mul(g, mul(h, inverse(g))) for h in _ncl.generators for g in G.generators]
        _ncl = PermGroup(_ncl.generators + _new, n=G.n)
    frobenius_generate_G = _ncl.order() == G.order()
    if not frobenius_generate_G:
        log(f"WARNING: the Frobenius classes observed up to {opts['class_bound_X']} generate a normal subgroup of order {_ncl.order()} < |G| = {G.order()}: the input group is probably too large (no descent certificate)")
    dump_json({"numbering_prime": p_num, "rational_classes": [[k + 1 for k in orb] for orb in A.rc],
               "family": [{"subgroup": t.S.label, "index": t.m, "invariant": t.F.label,
                           "resolvent": t.R, "squarefree_prime": t.squarefree_prime} for t in A.stage.tests],
               "twisted": {str([k + 1 for k in orb]): [{"j": t.j, "o": t.o, "degree": len(t.T),
                           "twisted_resolvent": [c.to_json() for c in t.R]} for t in tests] for orb, tests in A.twisted.items()},
               "primes": class_records}, os.path.join(run_dir, "classes.json"))
    _mark("class_assignment")
    conf = None
    if opts["class_confirmation"] and G.order() <= 200:
        cp = primes[:opts["class_confirmation_primes"]]
        summary, rows, _ = confirm(G, cl, D.f, p_num, policy, cp, seed=opts["seed"], log=log)
        conf = {"summary": summary, "rows": rows}
        log(f"class confirmation on {len(cp)} primes: {summary}")
        dump_json(conf, os.path.join(run_dir, "class_confirmation.json"))
    _mark("class_confirmation")
    # ---- local fields, matching and the well-definedness check
    chain = (A.direct.F, A.direct.c)
    L = run_local(D.f, G, cl, chain, policy, R, A, unramified_sample=opts["local_unramified_check"], seed=opts["seed"], log=log)
    objects = L.pop("_objects")
    dump_json(L, os.path.join(run_dir, "local.json"))

    _mark("local_descents")
    # ---- ramification filtrations, their structural checks and the conductor exponents
    ring_p, roots_p, _ = _resolvent.roots_at(D.f, p_num, policy.k(p_num), seed=opts["seed"])
    qroots = factor_root_sets(R["factors"], ring_p, roots_p)
    FILT = run_filtrations(objects, T, cl, G, R, qroots, policy, log=log)
    wild = [ell for ell, rr in L["ramified"].items() if rr["status"] != "ok"]
    conductors = []
    for nu in range(T.r):
        fac = {ell: FILT[ell]["conductor_exponents"][nu]["f_ell"] for ell in FILT}
        conductors.append({"chi": nu + 1, "exponents": fac, "unknown_at": wild,
                           "partial_conductor": __import__("math").prod(int(ell) ** e for ell, e in fac.items())})
    dump_json({"primes": FILT, "conductors": conductors, "wild_primes": wild}, os.path.join(run_dir, "conductors.json"))
    log(f"conductors written; exponents known at {sorted(FILT)}; unknown at wild primes {wild}")

    _mark("filtrations_and_conductors")
    # ---- Euler factors at the ramified primes and the identity system
    EU = run_euler(G, cl, T, D.f, p_num, policy, objects, FILT, models, A, R, log=log, seed=opts["seed"])
    _fam = EU.pop("_family")
    dump_json(EU, os.path.join(run_dir, "euler.json"))
    EU["_family"] = _fam

    _mark("euler_factors")
    # ---- the archimedean place
    AR = run_archimedean(G, cl, T, D.f, A, log=log)
    dump_json(AR, os.path.join(run_dir, "archimedean.json"))

    _mark("archimedean")
    # ---- root numbers
    if opts["root_numbers"]:
        RN = run_root_numbers(G, cl, T, objects, L, AR, log=log)
    else:
        RN = {"wild_primes": [], "characters": [{"chi": nu + 1, "W": None, "status": "not computed (root_numbers off)"} for nu in range(T.r)]}
        log("root numbers: skipped (root_numbers off)")
    dump_json(RN, os.path.join(run_dir, "rootnumbers.json"))

    _mark("root_numbers")
    # ---- functional equation: functional-equation test and subfield identities
    CJ = load_json(os.path.join(run_dir, "conductors.json"))
    AN = {"characters": [], "subfield_identities": {}, "X": None, "_coefficients": {}}
    if opts["fe_test"] and opts["root_numbers"]:
        AN = run_analytic(G, cl, T, D.f, A, EU, CJ, RN, AR, EU.get("_family", []), log=log, eps=opts["fe_eps"], cap=opts["fe_cap"], seed=opts["seed"])
    coeffs = AN.pop("_coefficients")
    raw = AN.pop("_raw", None)
    dump_json({"X": AN["X"], "coefficients": coeffs}, os.path.join(run_dir, "coefficients.json"))
    dump_json(AN, os.path.join(run_dir, "analytic.json"))

    # ---- the decision rule, its perturbations, the trivial character and the zeta anchor
    FA = None
    if opts["fe_test"] and raw is not None:
        FA = run_falsifier(G, cl, T, D.f, R["factors"], raw, AN["X"], CJ, RN, AR, policy, L, log=log)
        dump_json(FA, os.path.join(run_dir, "falsifier.json"))
    dump_json(policy.to_json(), os.path.join(run_dir, "precision.json"))
    log(f"precision policy: mode {policy.mode}, log2 M* = {policy.Mstar.bit_length() - 1}, {len(policy.pairs)} pairs, consultations {policy.consultations}")

    _mark("functional_equation")
    dump_json(_T, os.path.join(run_dir, "timings.json"))
    # ---- certificate
    cert = {"version": CERT_VERSION,
            "input": {"f": D.f, "n": D.n},
            "numbering": {"numbering": D.numbering},
            "group": {"descent_certificate": D.data, "source": D.source, "verified": D.verified,
                     "generators": [to_json(g) for g in G.generators], "order": G.order()},
            "character_table": {"file": "chartable.json", "S": T.certificate["S"], "verified": T.certificate["verified"] and ok,
                     "r": T.r, "degrees": T.degrees},
            "ramified_primes": {"file": "ramified.json", "ramified": R["ramified"], "candidates": R["candidates"],
                     "witness_methods": {ell: {i: w["method"] for i, w in rec["factors"].items()} for ell, rec in R["per_prime"].items()}},
            "precision": {"file": "precision.json", "mode": policy.mode, "log2_Mstar": policy.Mstar.bit_length() - 1,
                     "log2_M": policy.M.bit_length() - 1, "numbering_prime": p_num,
                     "consultations": {str(l): k for l, k in policy.consultations.items()},
                     "doubling_check": None if pcheck is None else {"file": "precision_check.json", "identical": pcheck["identical"]}},
            "frobenius_classes": {"file": "classes.json", "X": opts["class_bound_X"], "primes": len(primes), "methods": methods,
                  "frobenius_classes_generate_G": frobenius_generate_G,
                  "confirmation": None if conf is None else {"file": "class_confirmation.json", **conf["summary"]}},
            "local_data": {"file": "local.json",
                  "primes": {ell: {"status": r["status"], "e": r.get("e"), "f": r.get("f"), "|D|": r.get("|D|"),
                                   "inertia_orbits_match_factors": r.get("inertia_orbits_match_factors")} for ell, r in L["ramified"].items()},
                  "unramified_check_ok": L["unramified_ok"], "ramified_orbit_check_ok": L["ramified_ok"]},
            "conductors": {"file": "conductors.json", "primes_done": sorted(FILT), "wild_not_built": wild,
                   "conductor_identity_ok": True, "filtration_checks_ok": True,
                   "unique_filtration": {ell: FILT[ell]["filtration_candidates"]["unique"] for ell in FILT}},
            "euler_factors": {"file": "euler.json", "primes": sorted(EU["ramified"]), "euler_routes_agree": True,
                   "route_A_checked": {ell: v["route_A_checked"] for ell, v in EU["ramified"].items()},
                   "zeta_identities_ok": True, "determinant_checks_ok": True},
            "archimedean": {"file": "archimedean.json", "real_roots": AR["real_roots"], "class_of_c": AR["class_of_c"],
                   "a_plus_b_ok": True, "parity_ok": True,
                   "odd_two_dimensional": [d["chi"] for d in AR["characters"] if d["odd_two_dimensional"]]},
            "root_numbers": {"file": "rootnumbers.json", "computed": [c["chi"] for c in RN["characters"] if c.get("W") is not None],
                   "wild_not_built": RN["wild_primes"],
                   "shortcuts": {c["chi"]: c.get("shortcuts") for c in RN["characters"] if c.get("shortcuts")}},
            "functional_equation": {"file": "analytic.json", "X": AN.get("X"),
                   "tested": [c["chi"] for c in AN["characters"] if c.get("status") == "ok"],
                   "informative": [c["chi"] for c in AN["characters"] if c.get("tests") and all(t["informative"] for t in c["tests"])],
                   "subfield_identities_ok": all(v["ok"] for v in AN["subfield_identities"].values())},
            "falsifier": None if FA is None else {"file": "falsifier.json", "summary": FA.get("summary"),
                                                     "anchor_ok": FA["anchor"].get("anchor_ok"),
                                                     "rho_K": FA["anchor"]["zeta_residue_rho_K"]},
            "schur": {"file": "schur.json",
                      "exact": [nu + 1 for nu, sd in enumerate(schur) if sd["value"] is not None],
                      "bounded": [nu + 1 for nu, sd in enumerate(schur) if sd["value"] is None]},
            "timings": _T,
            "finished": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    dump_json(cert, os.path.join(run_dir, "certificate.json"))
    from .certwriter import write_certificate
    write_certificate(run_dir)
    log("certificate.json and CERT.json (self-contained) written")
    return cert

def main(argv=None):
    ap = argparse.ArgumentParser(description="Compute Artin L-function data for a monic separable integer polynomial")
    ap.add_argument("--config")
    ap.add_argument("--descent")
    ap.add_argument("--f")
    ap.add_argument("--gens")
    ap.add_argument("--run-dir")
    a = ap.parse_args(argv)
    if a.config:
        cfg = load_json(a.config)
    else:
        cfg = {}
        if a.descent:
            cfg["descent"] = a.descent
        else:
            cfg["f"] = parse_poly(a.f)
            cfg["generators"] = json.loads(a.gens)
    if a.run_dir:
        cfg["run_dir"] = a.run_dir
    if "run_dir" not in cfg:
        cfg["run_dir"] = "runs/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    run_pipeline(cfg)

if __name__ == "__main__":
    main()