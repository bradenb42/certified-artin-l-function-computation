"""The precision experiment: recording every value at two precisions.

Records the values that the precision policy governs at the policy precision
and at twice that precision, at the numbering prime, at a sample of
unramified primes and at every ramified prime, and reports whether the two
records agree.
"""
from __future__ import annotations
import random
from math import gcd

from .invariants import orbit_sum, set_product, stabilizer_in, Invariant
from .padic import count_Zl_roots, frobenius_perm, PrecisionExhausted
from .perm import symmetric, PermGroup, mul, inverse, cycle_type, to_json
from .precision import Policy
from .ramified import discriminant, valuation, factorint
from .resolvent import (coset_transversal, roots_at, resolvent, invariant_value,
                        alignment_single_level, evaluate_permuted, squarefree_certificate, divides_disc)

def is_prime(n):
    from .chartable import is_prime as ip
    return ip(n)

def default_chain(G: PermGroup):
    """The one-level chain S_n > G with the orbit-sum invariant on the
    base (1,...,n-1) of S_n, so that Stab_{S_n}(F) = G."""
    n = G.n
    Sn = symmetric(n)
    F = orbit_sum(n, G.elements(), list(range(n - 1)), label="F_G (orbit sum)")
    return Sn, F

TSCHIRNHAUS = [None, [0, 1, 1], [0, 1, 2], [0, 2, 1], [1, 1, 1, 1], [0, 1, -1, 1], [0, 3, 1, 1]]

def make_squarefree(F, transversal, f, p, policy, m, seed=0):
    """Replace F by F o T until the resolvent over the transversal is
    squarefree (disc != 0, decided exactly at the policy precision).  Returns the
    (possibly transformed) invariant and its resolvent at p."""
    for T in TSCHIRNHAUS:
        Ft = Invariant(F.n, F.terms, T, F.label + (f" o T{T}" if T else ""))
        trial = Policy(f, policy.Delta, mode=policy.mode)
        trial.register(Ft, m)
        k = trial.k(p)
        ring, roots, _ = roots_at(f, p, k, seed=seed)
        R, _ = resolvent(Ft, transversal, roots, ring)
        if squarefree_certificate(R) is not None:
            policy.register(Ft, m, label=f"{Ft.label}, m={m}")
            return Ft, R
    raise RuntimeError("no squarefree resolvent found among the Tschirnhaus candidates")

def find_group_in_numbering(f, abstract_gens, p, policy, seed=0, log=print):
    """Stand-in for the descent artifact when no descent certificate is supplied: among the
    conjugates H of <abstract_gens> in S_n, the Galois group in the numbering of the sorted
    p-adic roots is the H with F_H(alpha) integral and Res_{S_n,H,F_H} squarefree
    (Stauduhar step: a squarefree resolvent with an integral root forces Gal <= H, and
    |H| = |Gal| is assumed from the abstract input).  Exact at the policy precision."""
    n = len(f) - 1
    G0 = PermGroup(abstract_gens, n=n)
    Sn = symmetric(n)
    tried = set()
    # Frob_p in the numbering is known exactly from the residues: a necessary filter
    ring0, roots0, _ = roots_at(f, p, 1, seed=seed)
    from .padic import frobenius_perm
    phi_p = frobenius_perm(ring0, roots0)
    for tau in Sn.elements():
        H = PermGroup([mul(mul(tau, g), inverse(tau)) for g in G0.generators], n=n)
        if not H.contains(phi_p):
            continue
        key = frozenset(H.elements())
        if key in tried:
            continue
        tried.add(key)
        F = orbit_sum(n, H.elements(), list(range(n - 1)), label="F_G (orbit sum)")
        T_H = coset_transversal(Sn, H)
        m = len(T_H)
        for T in TSCHIRNHAUS:
            Ft = Invariant(n, F.terms, T, F.label)
            pol2 = Policy(f, policy.Delta, mode=policy.mode)
            pol2.register(Ft, m)
            k = pol2.k(p)
            ring, roots, _ = roots_at(f, p, k, seed=seed)
            vals = [evaluate_permuted(Ft, roots, tau2, ring) for tau2 in T_H]
            if ring.centered_int(vals[0]) is None:
                break
            R, _ = resolvent(Ft, T_H, roots, ring)
            if squarefree_certificate(R) is None:
                continue
            if ring.centered_int(vals[0]) is not None and abs(ring.centered_int(vals[0])) <= Ft.height_bound(policy.R):
                log(f"group located in the numbering at p={p}: conjugate with generators {[to_json(g) for g in H.generators]}")
                return H
            break
    raise RuntimeError("no conjugate of the abstract group has an integral invariant value: wrong group or numbering")

def auxiliary_pairs(G: PermGroup):
    """(F, H, label) with H = Stab_G(F): a point stabilizer and a 2-set stabilizer."""
    n = G.n
    out = []
    F1 = set_product(n, [0], label="x_1")
    H1 = PermGroup(stabilizer_in(F1, G.elements()) or [tuple(range(n))], n=n)
    out.append((F1, H1))
    if n >= 3:
        F2 = set_product(n, [0, 1], label="x_1 x_2")
        H2 = PermGroup(stabilizer_in(F2, G.elements()) or [tuple(range(n))], n=n)
        out.append((F2, H2))
    return out

def numbering_prime(Delta, avoid=(), coprime_to=1, orders=()):
    """Smallest odd prime not dividing Delta (nor coprime_to, e.g. the exponent of G), not in avoid,
    and, when possible within the first 300 candidates, with Phi_o irreducible mod p for every o in
    orders."""
    from .twisted import multiplicative_order
    from .cyclo import phi as _phi
    def ok_basic(p):
        return Delta % p and coprime_to % p and p not in avoid
    p, tried, first = 3, 0, None
    while True:
        if is_prime(p) and ok_basic(p):
            if first is None:
                first = p
            if all(multiplicative_order(p, o) == _phi(o) for o in orders if o > 2):
                return p
            tried += 1
            if tried > 300:
                return first
        p += 2

def record_values(f, G, Sn, F, c_pairs, policy, factor, p, unram, ram, seed=0):
    """All the precision policy-governed values at precision factor * k_ell.  c_pairs: list of (F_aux, H_aux)."""
    rec = {"factor": factor, "primes": {}}
    n = G.n
    T_G = coset_transversal(Sn, G)
    aux = [(Fa, coset_transversal(G, Ha), Fa.label) for Fa, Ha in c_pairs]
    Tcache = {}
    # numbering prime
    kp = policy.k(p) * factor
    ring, roots, facs = roots_at(f, p, kp, seed=seed)
    assert len(roots) == n
    cG = invariant_value(F, roots, ring)
    ResG, _ = resolvent(F, T_G, roots, ring)
    entry = {"k": kp, "c_G": cG, "Res_G": ResG,
             "disc_Res_G_nonzero": squarefree_certificate(ResG) is not None,
             "roots_mod_k": [ring.truncate(x, policy.k(p)) for x in roots], "frobenius": to_json(frobenius_perm(ring, roots))}
    for Fa, Ta, label in aux:
        R, _ = resolvent(Fa, Ta, roots, ring)
        entry[f"Res[{label}]"] = R
        entry[f"disc[{label}]_nonzero"] = squarefree_certificate(R) is not None
    rec["primes"][str(p)] = entry
    resolvents = {"Res_G": ResG}
    resolvents.update({f"Res[{label}]": entry[f"Res[{label}]"] for _, _, label in aux})
    resolvents["f"] = list(f)
    # unramified sample
    for ell in unram:
        k = policy.k(ell) * factor
        ring, roots, facs = roots_at(f, ell, k, seed=seed)
        assert len(roots) == n
        phi = frobenius_perm(ring, roots)
        surv = alignment_single_level(F, cG, T_G, roots, ring)
        e = {"k": k, "roots_mod_k": [ring.truncate(x, policy.k(ell)) for x in roots],
             "frobenius_beta": to_json(phi), "alignment_survivors": surv,
             "factorization_type": sorted(len(g) - 1 for g, _ in facs)}
        if len(surv) == 1:
            rho = T_G[surv[0]]
            sigma = mul(inverse(rho), mul(phi, rho))
            assert G.contains(sigma)
            e["frobenius_class"] = G.classes().class_of(sigma) + 1
            e["frobenius_cycle_type"] = list(cycle_type(sigma))
        for name, R in resolvents.items():
            e[f"Zl_roots[{name}]"] = _count(R, ell, k, policy, factor)
        rec["primes"][str(ell)] = e
    # ramified sample
    for ell in ram:
        v = valuation(abs(discriminant(f)), ell)
        k = policy.k(ell, v_disc=v) * factor
        ring, roots, facs = roots_at(f, ell, k, seed=seed)
        e = {"k": k, "v_disc": v, "simple_roots_mod_k": [ring.truncate(x, policy.k(ell, v_disc=v)) for x in roots],
             "factorization_type_mod_ell": sorted((len(g) - 1, m) for g, m in facs)}
        for name, R in resolvents.items():
            e[f"Zl_roots[{name}]"] = _count(R, ell, k, policy, factor)
        rec["primes"][str(ell)] = e
    return rec

def _count(R, ell, k, policy, factor):
    """T5 count at the given precision; if Panayi needs more (disc R not registered), the
    precision is raised by doubling, keeping the factor between the two records."""
    kk = k
    while True:
        try:
            return count_Zl_roots(R, ell, kk)
        except PrecisionExhausted:
            kk *= 2
            policy.extra.setdefault(str(ell), []).append({"k": kk // factor, "factor": factor})
            if kk > 10 ** 6:
                raise

def compare(rec1, rec2):
    diffs = []
    for ell, e1 in rec1["primes"].items():
        e2 = rec2["primes"][ell]
        for key, v1 in e1.items():
            if key == "k":
                continue
            v2 = e2.get(key)
            if v1 != v2:
                diffs.append({"prime": ell, "key": key, "at_k": v1, "at_2k": v2})
    return diffs

def precision_check(f, G, policy: Policy, n_unramified=6, seed=0, log=print):
    Delta = abs(discriminant(f))
    ram = sorted(factorint(Delta).keys())
    p = numbering_prime(Delta)
    Sn, F = default_chain(G)
    mG = Sn.order() // G.order()
    T_G = coset_transversal(Sn, G)
    F, RG = make_squarefree(F, T_G, f, p, policy, mG, seed)
    policy.register_disc("f", Delta)
    pairs = []
    for Fa, Ha in auxiliary_pairs(G):
        Ta = coset_transversal(G, Ha)
        Fa, Ra = make_squarefree(Fa, Ta, f, p, policy, len(Ta), seed)
        pairs.append((Fa, Ha))
    unram, q = [], p
    while len(unram) < n_unramified:
        q += 2
        while not is_prime(q):
            q += 2
        if Delta % q:
            unram.append(q)
    log(f"policy: log2 M* = {policy.Mstar.bit_length() - 1}, log2 M = {policy.M.bit_length() - 1}; numbering prime {p} (k={policy.k(p)}), unramified sample {unram}, ramified {ram}")
    rec1 = record_values(f, G, Sn, F, pairs, policy, 1, p, unram, ram, seed)
    # excluded primes of the chain resolvent: add the small ones to the sample
    RG = rec1["primes"][str(p)]["Res_G"]
    excluded = []
    q = 2
    while q < 200 and len(excluded) < 5:
        if is_prime(q) and Delta % q and q != p and divides_disc(RG, q):
            excluded.append(q)
        q += 1
    if excluded:
        log(f"excluded primes of Res_G in sample: {excluded[:5]}")
        rec1 = record_values(f, G, Sn, F, pairs, policy, 1, p, unram + excluded[:5], ram, seed)
    rec2 = record_values(f, G, Sn, F, pairs, policy, 2, p, unram + excluded[:5], ram, seed)
    diffs = compare(rec1, rec2)
    return {"numbering_prime": p, "unramified_sample": unram + excluded[:5], "excluded_primes_sampled": excluded[:5],
            "ramified": ram, "at_k": rec1, "at_2k": rec2, "differences": diffs, "identical": not diffs}