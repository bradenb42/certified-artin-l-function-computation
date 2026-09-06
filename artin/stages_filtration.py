"""Filtrations, conductor exponents and their identities at every prime.
"""
from __future__ import annotations
from .filtration import (Filtration, HardFailure, enumerate_candidates, discriminating_subgroups,
                         discriminating_resolvent_check)
from .local import _orbits
from .perm import from_json, mul, inverse
from .ramified import valuation

def factor_root_sets(f_factors, roots_ring, roots):
    """Indices of the global roots belonging to each irreducible factor of f over Q."""
    ring = roots_ring
    out = []
    for fi in f_factors:
        out.append([j for j, r in enumerate(roots) if ring.valuation(ring.eval_poly(fi, r)) >= ring.k])
    assert sorted(x for o in out for x in o) == list(range(len(roots)))
    return out

def run_filtrations(local_results, table, cl, G, ramified_json, qfactor_roots, policy, log=print, max_i=6):
    """local_results: dict ell -> (record, LocalGalois, D, I, Fr) from the local stage.  Returns per-prime records
    and the partial global conductors."""
    out = {}
    for ell, (rec, local, D, I, Fr) in local_results.items():
        rho = from_json(rec["matching"])
        Fl = Filtration(local, log=log)
        Fl.hasse_arf()
        n = len(local.roots)
        # D-orbits in the local numbering = Q_ell-factors; polygons for each
        Dorbs = _orbits(local.D, n)
        polys = Fl.polygons(Dorbs)
        conds = Fl.conductor_exponents(table, cl, rho)
        # discriminant valuations per factor of f.  Convention of the matching test F(beta^rho) = F(alpha):
        # alpha_i = beta_{rho(i)}, so the global root set orb corresponds to the local set rho(orb)
        prec = ramified_json["per_prime"].get(str(ell), {"factors": {}})
        factor_disc_valuations = {}
        for j in range(len(qfactor_roots)):
            w = prec["factors"].get(str(j))
            factor_disc_valuations[j] = w["v_disc_O"] if w else 0
        local_orbits = [[rho[a] for a in orb] for orb in qfactor_roots]
        ids = Fl.identities(table, cl, conds, G, rho, factor_disc_valuations, local_orbits)
        cands = enumerate_candidates(Fl.I_elems, Fl.D_elems, ell, [(fa["root"] - 1, fa["polygon"]) for fa in polys["factors"]], max_i, n)
        true_key = {str([x + 1 for x in s]): v for s, v in Fl.i_L.items()}
        contains_true = any(all(true_key.get(k) == v for k, v in c.items()) for c in cands) or not Fl.i_L
        subs = discriminating_subgroups(Fl)
        S = subs[0]
        dres = discriminating_resolvent_check(Fl, S, G.base())
        if not dres["relation_ok"]:
            raise HardFailure(f"discriminating resolvent relation fails at {ell}")
        out[str(ell)] = {"filtration": Fl.to_json(), "polygons": polys, "N_by_class": Fl.N_by_class,
                         "filtration_candidates": {"count": len(cands), "contains_true": contains_true, "unique": len(cands) == 1},
                         "discriminating_subgroup": {"order": S.order(), **dres},
                         "conductor_exponents": conds, "identities": ids}
        log(f"ell = {ell}: conductor exponents {[(c['chi'], c['f_ell']) for c in conds]}; conductor identity {ids['sum_chi1_f']} = {ids['v_ell_d_N']}; factor identities {[ (x['f_ell(pi)'], x['v_ell(disc O)']) for x in ids['factors']]}; filtration candidates {len(cands)} (true included: {contains_true})")
    return out