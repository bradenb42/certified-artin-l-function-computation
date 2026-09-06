"""Assignment of the Frobenius class at unramified primes.

Combines the rational-class tests, the twisted-resolvent refinement and the
direct route, choosing between them by cost, and confirms the result against
an explicitly constructed splitting field when the group is small enough.
"""
from __future__ import annotations
from math import gcd

from .classes import RationalClassStage, rational_classes, TSCHIRNHAUS
from .twisted import TwistedResolvent
from .direct import DirectRoute
from .stages_precision import default_chain, make_squarefree
from .invariants import Invariant
from .perm import symmetric, cycle_type
from .precision import Policy
from .resolvent import coset_transversal, roots_at, resolvent, invariant_value, evaluate_permuted, squarefree_certificate, divides_disc
from .ramified import factor_mod, factorint
from .chartable import is_prime

class ClassAssigner:
    def __init__(self, G, cl, f, p, policy: Policy, seed=0, log=print, twisted=True, twisted_work_budget=5e8):
        self.G, self.cl, self.f, self.p, self.policy, self.seed = G, cl, f, p, policy, seed
        Sn, F = default_chain(G)
        T_G = coset_transversal(Sn, G)
        F, RG = make_squarefree(F, T_G, f, p, policy, len(T_G), seed)
        ring, roots, _ = roots_at(f, p, policy.k(p), seed=seed)
        cG = invariant_value(F, roots, ring)
        self.direct = DirectRoute(G, cl, f, F, cG, RG, policy, seed)
        self.stage = RationalClassStage(G, cl, f, p, policy, seed, log)
        self.rc = self.stage.rc
        self.twisted = {}
        if twisted:
            for orb in self.rc:
                if len(orb) < 2:
                    continue
                k0 = orb[0]
                o = cl.orders[k0]
                # the cost analysis routing: the twisted resolvents cost ~ [G:C] evaluations at a precision of
                # m log2(2B+2) bits with B ~ phi(o)^{phi(o)/2} o^o R^{o deg m_t}; when that exceeds the
                # budget, the direct route (residue-field arithmetic) is used for this rational class
                import math as _m
                from .cyclo import phi as _phi
                deg_mt = sum(range(1, len(G.base()) + 1))
                m = G.order() // o
                bits = m * (_m.log2(2 * (_phi(o) ** (_phi(o) / 2) * o ** o) * policy.R ** (o * deg_mt) + 2))
                from .twisted import multiplicative_order as _mo
                from .ramified import factor_mod as _fm
                r_f = 1
                for gg, ee in _fm(f, p):
                    dd = len(gg) - 1; r_f = r_f * dd // gcd(r_f, dd)
                ro = _mo(p, o); R_ring = r_f * ro // gcd(r_f, ro)
                work = bits * R_ring ** 2 * m * o        # ~ bigint operations per twisted resolvent
                if work > twisted_work_budget:
                    log(f"twisted resolvents for rational class {[k+1 for k in orb]} skipped: estimated work {work:.1e} > budget {twisted_work_budget:.0e} (precision {int(bits)} bits, ring degree {R_ring}); direct route will be used")
                    continue
                try:
                    tests = [TwistedResolvent(G, cl, k0, j, f, p, policy, seed) for j in range(1, o) if gcd(j, o) == 1]
                except ValueError as ex:
                    log(f"twisted resolvents for rational class {[k+1 for k in orb]} skipped: {ex} (direct route will be used)")
                    continue
                self.twisted[tuple(orb)] = tests
                log(f"twisted resolvents for rational class {[k+1 for k in orb]}: {len(tests)} of degree {len(tests[0].T)} over Z[zeta_{o}]")

    def assign(self, ell, with_direct=False):
        rec = {"ell": ell}
        cands, ct, recs = self.stage.candidates(ell)
        rec["block"] = list(ct)
        rec["resolvent_candidates"] = [k + 1 for k in cands]
        rec["resolvent_tests"] = recs
        method = "resolvent"
        if len(cands) > 1:
            orb = next(o for o in self.rc if cands[0] in o)
            tests = self.twisted.get(tuple(orb), [])
            for t in tests:
                keep, info = t.consistent_classes(ell, cands)
                if keep is None:
                    rec.setdefault("cyclotomic_refinement", []).append({"j": t.j, "status": "excluded at ell"})
                    continue
                rec.setdefault("cyclotomic_refinement", []).append({"j": t.j, **info, "candidates": [k + 1 for k in keep]})
                cands = keep
                if len(cands) == 1:
                    break
            rec["refined_candidates"] = [k + 1 for k in cands]
            method = "cyclotomic" if len(cands) == 1 else "direct"
        if len(cands) != 1 or with_direct:
            k, info = self.direct.frobenius_class(ell)
            rec["direct_route"] = {"class": k + 1, **info}
            if len(cands) != 1:
                if k not in cands:
                    raise RuntimeError(f"direct route class {k+1} outside candidates {cands} at {ell}")
                cands = [k]
            elif k != cands[0]:
                raise RuntimeError(f"inconsistent: the resolvent and cyclotomic tests gave {cands[0]+1}, direct {k+1} at {ell}")
        rec["class"] = cands[0] + 1
        rec["method"] = method
        return rec

# ------------------------------------------------------------------ explicit splitting field

class SplittingField:
    """N = Q(theta), theta = sum c_i alpha_i, minimal polynomial of degree |G| recovered
    exactly at the numbering prime."""
    def __init__(self, G, f, p, policy, seed=0, log=print):
        n = G.n
        self.G = G
        for attempt in range(6):
            c = [(i + 1) + attempt * (i * i % 5) for i in range(n)]
            F = Invariant(n, {tuple(1 if i == j else 0 for i in range(n)): c[j] for j in range(n)}, label="theta")
            trial = Policy(f, policy.Delta, mode=policy.mode)
            trial.register(F, G.order())
            ring, roots, _ = roots_at(f, p, trial.k(p), seed=seed)
            T = G.elements()
            R, vals = resolvent(F, T, roots, ring)
            if squarefree_certificate(R) is not None:
                self.c, self.minpoly = c, R
                policy.register(F, G.order(), label="theta (splitting field)")
                log(f"splitting field: theta = sum c_i alpha_i, c = {c}, minimal polynomial of degree {len(R)-1}")
                return
        raise RuntimeError("no primitive element among the candidates")

    def residue_degrees(self, ell):
        """Degrees of the irreducible factors of the minimal polynomial mod ell (ell not
        dividing its discriminant): all equal to the order of Frob_ell."""
        if divides_disc(self.minpoly, ell):
            return None
        facs = factor_mod(self.minpoly, ell)
        return sorted(len(g) - 1 for g, e in facs for _ in range(e))

# ------------------------------------------------------------------ the experiment

def confirm(G, cl, f, p, policy, primes, seed=0, log=print, splitting_field=True):
    A = ClassAssigner(G, cl, f, p, policy, seed, log)
    SF = SplittingField(G, f, p, policy, seed, log) if splitting_field else None
    rc_of = {k: i for i, orb in enumerate(A.rc) for k in orb}
    rows, ok_all = [], True
    for ell in primes:
        rec = A.assign(ell, with_direct=True)
        true = rec["direct_route"]["class"] - 1
        # (i) assigned class equals the direct class (checked inside assign); (ii) rational-class stage
        # leaves exactly the rational class of the true class
        a5 = sorted(k - 1 for k in rec["resolvent_candidates"])
        rec["rational_class_prediction_ok"] = (a5 == sorted(A.rc[rc_of[true]]))
        if SF is not None:
            degs = SF.residue_degrees(ell)
            rec["splitting_field_degrees"] = degs
            rec["splitting_field_ok"] = degs is None or all(d == cl.orders[true] for d in degs)
        else:
            rec["splitting_field_ok"] = True
        ok = rec["rational_class_prediction_ok"] and rec["splitting_field_ok"]
        ok_all = ok_all and ok
        rows.append(rec)
    summary = {"primes": len(rows),
               "methods": {m: sum(1 for r in rows if r["method"] == m) for m in ("resolvent", "cyclotomic", "direct")},
               "rational_class_prediction_ok": all(r["rational_class_prediction_ok"] for r in rows),
               "splitting_field_ok": all(r["splitting_field_ok"] for r in rows),
               "all_ok": ok_all}
    return summary, rows, A