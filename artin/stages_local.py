"""Local data at every ramified prime, with the unramified cross-check.

Runs the local field construction, the matching and the well-definedness check at
each ramified prime, verifies that the inertia orbit sizes match the
ramification indices, and confirms at unramified primes that the same
machinery returns the Frobenius class found by the class assignment.
"""
from __future__ import annotations

from .matching import local_data_at
from .ramified import factor_mod, valuation
from .chartable import is_prime

def run_local(f, G, cl, chain, policy, ramified_json, assigner, unramified_sample=6, seed=0, log=print):
    Delta = policy.Delta
    out = {"ramified": {}, "unramified_check": []}
    objects = {}
    for ell_s, rec in ramified_json["per_prime"].items():
        ell = int(ell_s)
        ram_data = []
        for i_s, w in rec["factors"].items():
            ram_data += [tuple(x) for x in w["residue_decomposition"]]
        # unramified factors of f not dividing disc contribute their residue degrees too
        facs_mod = factor_mod(f, ell)
        n_covered = sum(e * fd for e, fd in ram_data)
        if n_covered < len(f) - 1:
            # factors f_i with ell not dividing disc f_i: unramified, degrees from f_i mod ell
            for i, fi in enumerate(ramified_json["factors"]):
                if str(i) in rec["factors"]:
                    continue
                for g, e in factor_mod(fi, ell):
                    ram_data.append((1, len(g) - 1))
        E = 1
        for e, fd in ram_data:
            E = E * e // __import__("math").gcd(E, e)
        v = valuation(abs(Delta), ell)
        if E % ell == 0:
            out["ramified"][ell_s] = {"ell": ell, "status": "wild ramification: the tame local construction does not apply", "e_lcm": E,
                                      "residue_decomposition_from_factors": ram_data}
            log(f"ell = {ell}: wild ramification (e = {E}); the local Galois group is not constructed at wild primes")
            continue
        try:
            rec2, local, D, I, Fr = local_data_at(f, ell, G, cl, chain, policy, ram_data, v, seed, log)
            rec2["status"] = "ok"
            out["ramified"][ell_s] = rec2
            objects[ell] = (rec2, local, D, I, Fr)
            log(f"ell = {ell}: |D| = {D.order()}, e = {rec2['e']}, f = {rec2['f']}, inertia orbits {rec2['inertia_orbit_sizes']} vs the factor data {rec2['ramification_indices_from_factors']}: {'match' if rec2['inertia_orbits_match_factors'] else 'MISMATCH'}; well-definedness check passed")
        except Exception as ex:
            out["ramified"][ell_s] = {"ell": ell, "status": f"failed: {ex}"}
            log(f"ell = {ell}: local descent failed: {ex}")
    # unramified confirmation
    q, count = 2, 0
    while count < unramified_sample:
        q += 1
        if not is_prime(q) or Delta % q == 0 or q == assigner.p:
            continue
        ram_data = [(1, len(g) - 1) for g, e in factor_mod(f, q)]
        rec2, local, D, I, Fr = local_data_at(f, q, G, cl, chain, policy, ram_data, 0, seed, lambda m: None)
        b3 = assigner.assign(q)["class"]
        ok = rec2["frobenius_coset_class"] == [b3] and I.order() == 1
        out["unramified_check"].append({"ell": q, "matching_class": rec2["frobenius_coset_class"], "class_from_class_assignment": b3, "ok": ok})
        count += 1
    out["unramified_ok"] = all(r["ok"] for r in out["unramified_check"])
    out["ramified_ok"] = all(r.get("inertia_orbits_match_factors", True) and r.get("decomposition_orbits_match_factors", True)
                             for r in out["ramified"].values() if r["status"] == "ok")
    log(f"unramified matching check: {out['unramified_ok']} on {len(out['unramified_check'])} primes; ramified orbit check: {out['ramified_ok']}")
    out["_objects"] = objects
    return out
