import json, math, os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pytest
from artin.perm import PermGroup, symmetric, alternating, from_cycles
from artin.cyclo import Cyc
from artin.chartable import CharacterTable, verify_certificate
from artin.schur import candidate_pairs, schur_data, build_model
from artin.ramified import ramified_primes, verify_ramified
from artin.run import run_pipeline
from artin.verify import verify_run

def q8():
    q = ['1','-1','i','-i','j','-j','k','-k']
    tab = {('1','1'):(1,'1'),('1','i'):(1,'i'),('1','j'):(1,'j'),('1','k'):(1,'k'),
           ('i','1'):(1,'i'),('i','i'):(-1,'1'),('i','j'):(1,'k'),('i','k'):(-1,'j'),
           ('j','1'):(1,'j'),('j','i'):(-1,'k'),('j','j'):(-1,'1'),('j','k'):(1,'i'),
           ('k','1'):(1,'k'),('k','i'):(1,'j'),('k','j'):(-1,'i'),('k','k'):(-1,'1')}
    def qmul(a, b):
        sa = -1 if a.startswith('-') else 1; sb = -1 if b.startswith('-') else 1
        s, r = tab[(a.strip('-'), b.strip('-'))]
        s *= sa * sb
        return ('-' if s < 0 else '') + r
    perm = lambda x: tuple(q.index(qmul(x, y)) for y in q)
    return PermGroup([perm('i'), perm('j')])

def test_perm_orders():
    for n in range(3, 9):
        assert symmetric(n).order() == math.factorial(n)
        assert alternating(n).order() == math.factorial(n) // 2
    assert symmetric(5).classes().r == 7
    assert alternating(5).classes().r == 5

def test_cyclo():
    z = Cyc.zeta(5); s = z + z ** 4
    assert (s * s + s - 1).is_zero()
    assert s.conductor() == 5 and s.conj() == s
    w = Cyc.zeta(8) + Cyc.zeta(8, 7)
    assert (w * w) == 2 and w.inverse() * w == 1

@pytest.mark.parametrize("G,degrees", [
    (symmetric(3), [1, 1, 2]), (symmetric(4), [1, 1, 2, 3, 3]), (alternating(5), [1, 3, 3, 4, 5]),
    (q8(), [1, 1, 1, 1, 2]), (symmetric(5), [1, 1, 4, 4, 5, 5, 6])])
def test_chartable(G, degrees):
    cl = G.classes()
    T = CharacterTable(cl)
    assert T.degrees == degrees
    assert T.certificate["verified"]
    assert verify_certificate(cl, T.certificate)[0]
    assert sum(d * d for d in T.degrees) == G.order()

def test_a5_values():
    T = CharacterTable(alternating(5).classes())
    # a degree-3 character takes the values (1 +- sqrt5)/2 on the two classes of 5-cycles
    row = T.values[1]
    vals = [v for k, v in enumerate(row) if T.cl.orders[k] == 5]
    assert all(not v.is_rational() for v in vals)
    assert (vals[0] + vals[1]) == 1 and (vals[0] * vals[1]) == -1

def test_schur_q8():
    G = q8(); cl = G.classes(); T = CharacterTable(cl); pairs = candidate_pairs(G, cl)
    sd = [schur_data(T, nu, pairs) for nu in range(T.r)]
    assert [s["value"] for s in sd] == [1, 1, 1, 1, 2]
    M = build_model(G, T, 4, pairs)
    assert M["status"] == "ok" and M["trace_check"] and M["dimension"] == 2

def test_models_s5():
    G = symmetric(5); cl = G.classes(); T = CharacterTable(cl); pairs = candidate_pairs(G, cl)
    for nu in range(T.r):
        assert schur_data(T, nu, pairs)["value"] == 1
        M = build_model(G, T, nu, pairs)
        assert M["status"] == "ok" and M["trace_check"] and M["dimension"] == T.degrees[nu]

@pytest.mark.parametrize("f,ram", [
    ([1, 0, 0, 0, 1], [2]), ([-2, 0, 0, 1], [2, 3]), ([-5, 0, 1], [5]),
    ([-8, -2, -1, 1], [503]), ([-1, -1, 0, 0, 0, 1], [19, 151]), ([1, 1, 1, 1, 1, 1, 1], [7]),
    ([-6, 0, -1, 0, 1], [2, 3]), ([9, 0, -2, 0, 1], [2])])
def test_ramified(f, ram):
    R = ramified_primes(f)
    assert R["ramified"] == ram
    assert verify_ramified(R)[0]

def test_dedekind_example_decomposition():
    R = ramified_primes([-8, -2, -1, 1])
    w = R["per_prime"]["2"]["factors"]["0"]
    assert w["method"] == "Round 2" and w["index_valuation"] == 1
    assert sorted(w["residue_decomposition"]) == [(1, 1), (1, 1), (1, 1)]

def test_run_and_verify():
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [-1, -1, 0, 0, 0, 1], "generators": [[2, 3, 4, 5, 1], [2, 1, 3, 4, 5]], "run_dir": d,
               "options": {"fe_test": False, "class_confirmation_primes": 4, "precision_check_unramified": 2}}
        cert = run_pipeline(cfg)
        assert os.path.exists(os.path.join(d, "config.json"))
        c = json.load(open(os.path.join(d, "config.json")))
        assert "b1" in c and "started" in c["b1"] and c["options"]["model_dim_limit"] > 0
        assert cert["character_table"]["verified"] and cert["ramified_primes"]["ramified"] == [19, 151]
        rep = verify_run(d)
        assert rep["ok"]

# ---------------------------------------------------------------- the precision stage
from artin.precision import Policy, PolicyViolation
from artin.padic import count_Zl_roots
from artin.stages_precision import precision_check, find_group_in_numbering, numbering_prime
from artin.ramified import discriminant

def test_panayi():
    assert count_Zl_roots([-17, 0, 1], 2, 10) == 2
    assert count_Zl_roots([-5, 0, 1], 2, 10) == 0
    assert count_Zl_roots([-8, -2, -1, 1], 2, 10) == 3
    assert count_Zl_roots([-2, 0, 0, 1], 3, 10) == 0

def test_policy_monotone_and_violation():
    from artin.invariants import set_product
    pol = Policy([-1, -1, 0, 0, 0, 1], 2869, mode="sharp")
    k0 = pol.k(7)
    pol.register(set_product(5, [0, 1]), 10)
    assert pol.k(7) >= k0
    with pytest.raises(PolicyViolation):
        pol.check(7, 1)

@pytest.mark.parametrize("f,gens,mode", [
    ([1, 0, 0, 0, 1], [(1, 0, 3, 2), (2, 3, 0, 1)], "sharp"),
    ([1, 0, 0, 0, 1], [(1, 0, 3, 2), (2, 3, 0, 1)], "conservative"),
    ([-8, -2, -1, 1], [(1, 2, 0), (1, 0, 2)], "sharp"),
    ([-2, 0, 0, 0, 1], [(1, 2, 3, 0), (0, 3, 2, 1)], "conservative")])
def test_doubling(f, gens, mode):
    D = discriminant(f)
    pol = Policy(f, D, mode=mode)
    p = numbering_prime(abs(D))
    G = find_group_in_numbering(f, gens, p, pol, log=lambda m: None)
    assert G.order() == PermGroup(gens).order()
    R = precision_check(f, G, pol, n_unramified=3, log=lambda m: None)
    assert R["identical"], R["differences"][:2]
    for ell, e in R["at_k"]["primes"].items():
        if "alignment_survivors" in e:
            assert len(e["alignment_survivors"]) == 1
            assert e["frobenius_cycle_type"] == sorted(e["factorization_type"], reverse=True)

def test_group_placed_before_table():
    """Input generators not in the numbering: the run relocates G first, and every
    recorded object (group.json, descent.json, certificate G0.2) is in the numbering."""
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [-2, 0, 0, 0, 1], "generators": [[2, 3, 4, 1], [1, 4, 3, 2]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "fe_test": False}}
        cert = run_pipeline(cfg)
        g = json.load(open(os.path.join(d, "group.json")))
        desc = json.load(open(os.path.join(d, "descent.json")))
        assert g["generators"] == cert["group"]["generators"] == desc["group"]["generators"]
        assert cert["precision"]["doubling_check"]["identical"]
        assert verify_run(d)["ok"]

def test_char2_extension_roots():
    import random
    from artin.padic import GF, irreducible_poly, roots_in_GF, _eval_over
    rng = random.Random(3)
    K = GF(2, irreducible_poly(2, 10, rng))
    h = irreducible_poly(2, 5, rng)
    rt = roots_in_GF(h, K, rng)
    assert len(rt) == 5 and all(K.is_zero(_eval_over(h, a, K)) for a in rt)

# ---------------------------------------------------------------- the class-assignment stage
from artin.stages_classes import confirm

@pytest.mark.parametrize("f,gens,nprimes", [
    ([1, 3, -3, -4, 1, 1], [(1, 2, 3, 4, 0)], 8),                       # C5, conductor 11
    ([12, -5, 0, 0, 0, 1], [(1, 2, 3, 4, 0), (0, 4, 3, 2, 1)], 10),      # D5
    ([-2, 0, 0, 0, 1], [(1, 2, 3, 0), (0, 3, 2, 1)], 10),                # D4: a block with two rational classes
    ([16, 20, 0, 0, 0, 1], [(1, 2, 0, 3, 4), (1, 2, 3, 4, 0)], 8)])      # the separating-subgroup construction, order 60
def test_class_assignment(f, gens, nprimes):
    from artin.chartable import is_prime
    D = discriminant(f)
    G0 = PermGroup(gens)
    pol = Policy(f, D, mode="sharp")
    p = numbering_prime(abs(D), coprime_to=G0.classes().exponent)
    G = find_group_in_numbering(f, gens, p, pol, log=lambda m: None)
    cl = G.classes()
    primes, q = [], p
    while len(primes) < nprimes:
        q += 2
        while not is_prime(q):
            q += 2
        if D % q:
            primes.append(q)
    S, rows, A = confirm(G, cl, f, p, pol, primes, log=lambda m: None)
    assert S["all_ok"], S
    # the rational-class stage alone leaves exactly the rational class
    assert all(r["rational_class_prediction_ok"] for r in rows)
    # a block with several rational classes was actually separated by the the separating-subgroup construction family for D4
    if G.order() == 8:
        assert len(A.stage.tests) >= 1 and all(len(r["resolvent_candidates"]) == 1 for r in rows)

def test_run_with_classes():
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [12, -5, 0, 0, 0, 1], "generators": [[2, 3, 4, 5, 1], [1, 5, 4, 3, 2]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 60, "class_confirmation_primes": 5, "fe_test": False}}
        cert = run_pipeline(cfg)
        assert cert["frobenius_classes"]["primes"] > 0 and cert["frobenius_classes"]["confirmation"]["all_ok"]
        assert verify_run(d)["ok"]

# ---------------------------------------------------------------- the local stage
from artin.local import LocalGalois
from artin.matching import local_data_at
from artin.stages_precision import default_chain, make_squarefree
from artin.resolvent import coset_transversal, roots_at, invariant_value
from artin.ramified import ramified_primes, valuation, factor_mod

def _setup(f, gens):
    D = discriminant(f)
    pol = Policy(f, D, mode="sharp")
    p = numbering_prime(abs(D), coprime_to=PermGroup(gens).classes().exponent)
    G = find_group_in_numbering(f, gens, p, pol, log=lambda m: None)
    cl = G.classes()
    Sn, F = default_chain(G)
    T = coset_transversal(Sn, G)
    F, RG = make_squarefree(F, T, f, p, pol, len(T))
    ring, roots, _ = roots_at(f, p, pol.k(p))
    c = invariant_value(F, roots, ring)
    return D, pol, p, G, cl, (F, c)

def test_local_tame_x5_x_1():
    f, gens = [-1, -1, 0, 0, 0, 1], [(1, 2, 3, 4, 0), (1, 0, 2, 3, 4)]
    D, pol, p, G, cl, chain = _setup(f, gens)
    R = ramified_primes(f)
    for ell, want_e, want_f in [(19, 2, 3), (151, 2, 2)]:
        rd = [tuple(x) for x in R["per_prime"][str(ell)]["factors"]["0"]["residue_decomposition"]]
        rec, local, Dl, Il, Fr = local_data_at(f, ell, G, cl, chain, pol, rd, valuation(abs(D), ell), log=lambda m: None)
        assert (rec["e"], rec["f"]) == (want_e, want_f)
        assert rec["inertia_orbits_match_factors"] and rec["decomposition_orbits_match_factors"]
        assert G.contains(Fr) and all(G.contains(g) for g in Il.generators)

def test_local_c5_at_11_and_unramified_matches_b3():
    f, gens = [1, 3, -3, -4, 1, 1], [(1, 2, 3, 4, 0)]
    D, pol, p, G, cl, chain = _setup(f, gens)
    rec, local, Dl, Il, Fr = local_data_at(f, 11, G, cl, chain, pol, [(5, 1)], 4, log=lambda m: None)
    assert rec["e"] == 5 and rec["f"] == 1 and rec["inertia_orbit_sizes"] == [5]
    from artin.stages_classes import ClassAssigner
    A = ClassAssigner(G, cl, f, p, pol, log=lambda m: None)
    for ell in [7, 13, 17]:
        rd = [(1, len(g) - 1) for g, e in factor_mod(f, ell)]
        rec, local, Dl, Il, Fr = local_data_at(f, ell, G, cl, chain, pol, rd, 0, log=lambda m: None)
        assert Il.order() == 1 and rec["frobenius_coset_class"] == [A.assign(ell)["class"]]

def test_run_with_local():
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [-8, -2, -1, 1], "generators": [[2, 3, 1], [2, 1, 3]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 40, "class_confirmation_primes": 4, "local_unramified_check": 3, "fe_test": False}}
        cert = run_pipeline(cfg)
        assert cert["local_data"]["primes"]["503"]["e"] == 2 and cert["local_data"]["primes"]["2"]["e"] == 1
        assert cert["local_data"]["unramified_check_ok"] and cert["local_data"]["ramified_orbit_check_ok"]
        assert verify_run(d)["ok"]

# ---------------------------------------------------------------- the filtration stage
from artin.filtration import Filtration, enumerate_candidates, HardFailure
from artin.perm import from_cycles, from_json

def test_a12_enumeration_synthetic():
    s = from_cycles(4, [[0, 1, 2, 3]]); C4 = PermGroup([s]).elements()
    assert len(enumerate_candidates(C4, C4, 2, [(0, [2, 2, 6])], 8, 4)) == 1      # Hasse-Arf consistent
    assert enumerate_candidates(C4, C4, 2, [(0, [2, 2, 3])], 8, 4) == []           # rejected by Hasse-Arf
    V4 = PermGroup([from_cycles(4, [[0, 1], [2, 3]]), from_cycles(4, [[0, 2], [1, 3]])]).elements()
    assert len(enumerate_candidates(V4, V4, 2, [(0, [2, 2, 4])], 6, 4)) == 3       # the filtration-uniqueness theory: polygons do not label

def test_conductors_x5_x_1():
    f, gens = [-1, -1, 0, 0, 0, 1], [(1, 2, 3, 4, 0), (1, 0, 2, 3, 4)]
    D, pol, p, G, cl, chain = _setup(f, gens)
    T = CharacterTable(cl)
    R = ramified_primes(f)
    for ell in (19, 151):
        rd = [tuple(x) for x in R["per_prime"][str(ell)]["factors"]["0"]["residue_decomposition"]]
        rec, local, Dl, Il, Fr = local_data_at(f, ell, G, cl, chain, pol, rd, valuation(abs(D), ell), log=lambda m: None)
        Fl = Filtration(local, log=lambda m: None)
        Fl.hasse_arf()
        conds = Fl.conductor_exponents(T, cl, from_json(rec["matching"]))
        ids = Fl.identities(T, cl, conds, G, from_json(rec["matching"]), {0: 1}, [list(range(5))])
        assert ids["sum_chi1_f"] == ids["v_ell_d_N"] == 60
        # sign character of S_5 has exponent 1 at each of 19, 151 (quadratic field of disc 19*151)
        assert conds[1]["f_ell"] == 1 and conds[0]["f_ell"] == 0

def test_run_reducible_with_conductors():
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [65, 0, -18, 0, 1], "generators": [[2, 1, 3, 4], [1, 2, 4, 3]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 40, "class_confirmation_primes": 3, "local_unramified_check": 2, "fe_test": False}}
        cert = run_pipeline(cfg)
        assert cert["conductors"]["primes_done"] == ["13", "2", "5"]
        assert verify_run(d)["ok"]
        C = json.load(open(os.path.join(d, "conductors.json")))
        # (x^2-5)(x^2-13): the four characters of V4 have conductors 1, 5, 13, 65
        assert sorted(c["partial_conductor"] for c in C["conductors"]) == [1, 5, 13, 65]

# ---------------------------------------------------------------- the Euler-factor stage
def test_run_s3_tame_e3_euler_and_identities():
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [14, -7, 0, 1], "generators": [[2, 3, 1], [2, 1, 3]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 40, "class_confirmation_primes": 3, "local_unramified_check": 2, "fe_test": False}}
        cert = run_pipeline(cfg)
        E = json.load(open(os.path.join(d, "euler.json")))
        # at 7 (e = 3, I = C_3) the degree-2 character has V^I = 0; at 5 (e = 2) it has dim 1
        assert E["ramified"]["7"]["characters"]["3"]["dim_V^I"] == 0
        assert E["ramified"]["5"]["characters"]["3"]["dim_V^I"] == 1
        assert E["ramified"]["7"]["characters"]["3"]["route_A"]["agrees_with_C"]
        assert cert["euler_factors"]["route_A_checked"]["7"] == [1, 2, 3]
        assert all(x["ok"] for ell in E["zeta_identities"] for x in E["zeta_identities"][ell])
        assert verify_run(d)["ok"]

def test_euler_routes_agree_s4():
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [-1, -1, 0, 0, 1], "generators": [[2, 3, 4, 1], [2, 1, 3, 4]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 40, "class_confirmation_primes": 3, "local_unramified_check": 2, "fe_test": False}}
        cert = run_pipeline(cfg)
        E = json.load(open(os.path.join(d, "euler.json")))
        assert cert["euler_factors"]["route_A_checked"]["283"] == [1, 2, 3]
        assert E["ramified"]["283"]["product_identity_ok"]
        assert verify_run(d)["ok"]

# ---------------------------------------------------------------- the archimedean stage
from artin.archimedean import real_root_count, isolate_real_roots

def test_sturm():
    assert real_root_count([-1, -1, 0, 0, 0, 1]) == 1          # x^5-x-1
    assert real_root_count([1, 0, 0, 0, 1]) == 0               # x^4+1
    assert real_root_count([1, 3, -3, -4, 1, 1]) == 5          # totally real C5 quintic
    assert real_root_count([-2, 0, 0, 0, 1]) == 2
    iso = isolate_real_roots([-2, 0, 0, 0, 1])
    assert len(iso) == 2 and all(a < b for a, b in iso)

def test_archimedean_run_s4_and_reducible():
    for f, gens, want_odd in [([-1, -1, 0, 0, 1], [[2, 3, 4, 1], [2, 1, 3, 4]], [3]),
                              ([65, 0, -18, 0, 1], [[2, 1, 3, 4], [1, 2, 4, 3]], [])]:
        with tempfile.TemporaryDirectory() as d:
            cfg = {"f": f, "generators": gens, "run_dir": d,
                   "options": {"precision_check_unramified": 2, "class_bound_X": 40, "class_confirmation_primes": 3, "local_unramified_check": 2, "fe_test": False}}
            cert = run_pipeline(cfg)
            AR = json.load(open(os.path.join(d, "archimedean.json")))
            assert all(c["a+b=chi(1)"] for c in AR["characters"]) and AR["permutation_check"]["ok"]
            assert cert["archimedean"]["odd_two_dimensional"] == want_odd
            assert verify_run(d)["ok"]

# ---------------------------------------------------------------- the root-number stage
def test_root_numbers_quadratic_and_quartic():
    import cmath
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [65, 0, -18, 0, 1], "generators": [[2, 1, 3, 4], [1, 2, 4, 3]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 30, "class_confirmation_primes": 3, "local_unramified_check": 2, "fe_test": False}}
        cert = run_pipeline(cfg)
        RJ = json.load(open(os.path.join(d, "rootnumbers.json")))
        assert all(c["W_exact"] == "i^0" for c in RJ["characters"])      # quadratic characters: W = 1 at 5, 13 (finite candidate set)
        assert verify_run(d)["ok"]
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [1, 1, 1, 1, 1], "generators": [[2, 3, 4, 1]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 30, "class_confirmation_primes": 3, "local_unramified_check": 2, "fe_test": False}}
        cert = run_pipeline(cfg)
        RJ = json.load(open(os.path.join(d, "rootnumbers.json")))
        # independent Dirichlet computation for the two quartic characters mod 5
        z = cmath.exp(2j * cmath.pi / 5)
        want = set()
        for c2 in (1j, -1j):
            chi = {1: 1, 2: c2, 4: c2 ** 2, 3: c2 ** 3}
            tau = sum(chi[a] * z ** a for a in chi)
            want.add(complex(round((tau / (1j * 5 ** 0.5)).real, 6), round((tau / (1j * 5 ** 0.5)).imag, 6)))
        got = {complex(round(c["W_complex"][0], 6), round(c["W_complex"][1], 6)) for c in RJ["characters"] if not c["real"]}
        assert got == want
        assert verify_run(d)["ok"]

def test_root_numbers_s4_orthogonal():
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [-1, -1, 0, 0, 1], "generators": [[2, 3, 4, 1], [2, 1, 3, 4]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 30, "class_confirmation_primes": 3, "local_unramified_check": 2, "fe_test": False}}
        cert = run_pipeline(cfg)
        assert all(v.get("orthogonal_W=1") for v in cert["root_numbers"]["shortcuts"].values())
        assert verify_run(d)["ok"]

# ---------------------------------------------------------------- the functional-equation stage
from artin.analytic import Kernel, closed_form

def test_kernels_against_closed_forms():
    for a, b in [(1, 0), (0, 1), (1, 1), (2, 0)]:
        K = Kernel(a, b); cf = closed_form(a, b)
        assert max(abs(K(x) - cf(x)) for x in (0.2, 1.0, 2.5)) < 1e-12
        assert all(abs(K(x)) <= K.bound(x) * 1.001 + 1e-15 for x in (0.5, 1.0, 3.0))

def test_functional_equation_c4_and_v4():
    for f, gens in [([1, 1, 1, 1, 1], [[2, 3, 4, 1]]), ([65, 0, -18, 0, 1], [[2, 1, 3, 4], [1, 2, 4, 3]])]:
        with tempfile.TemporaryDirectory() as d:
            cfg = {"f": f, "generators": gens, "run_dir": d,
                   "options": {"precision_check_unramified": 2, "class_bound_X": 30, "class_confirmation_primes": 3, "local_unramified_check": 2}}
            cert = run_pipeline(cfg)
            AN = json.load(open(os.path.join(d, "analytic.json")))
            for c in AN["characters"]:
                assert c["status"] == "ok" and all(t["ok"] and t["informative"] and t["relative"] < 1e-9 for t in c["tests"])
                assert c["defect_with_minus_W"] > 0.1
            assert all(v["ok"] for v in AN["subfield_identities"].values())
            FA = json.load(open(os.path.join(d, "falsifier.json")))
            for c in FA["characters"]:
                assert not any(e["rejected"] for e in c["true"])
                assert all(e["verdict"] == "rejected" for e in c["wrong_W"])
                assert all(e["verdict"] == "rejected" for e in c["wrong_frobenius"])
            assert all(x["ok"] for x in FA["trivial"])
            assert FA["anchor"]["anchor_ok"] and FA["anchor"]["relative_difference"] < 1e-10
            if f == [1, 1, 1, 1, 1]:
                assert abs(FA["anchor"]["zeta_residue_rho_K"] - 0.339837278241) < 1e-9   # class number formula for Q(zeta_5)
            assert verify_run(d)["ok"]

# ---------------------------------------------------------------- the certificate stage
from artinverify.verify import Verifier
from artinverify.mutate import mutations

def test_standalone_verifier_accepts_and_rejects_mutations():
    with tempfile.TemporaryDirectory() as d:
        cfg = {"f": [65, 0, -18, 0, 1], "generators": [[2, 1, 3, 4], [1, 2, 4, 3]], "run_dir": d,
               "options": {"precision_check_unramified": 2, "class_bound_X": 30, "class_confirmation_primes": 3,
                           "local_unramified_check": 2, "fe_test": False}}
        run_pipeline(cfg)
        cert = json.load(open(os.path.join(d, "CERT.json")))
        ok, rep = Verifier(cert).run()
        assert ok, rep
        muts = mutations(cert)
        assert set(muts) == {"filtration", "conductor", "frobenius", "euler", "rootnumber"}
        for name, c in muts.items():
            ok, rep = Verifier(c).run()
            assert not ok, name

# ---------------------------------------------------------------- the newform check
def test_check3_eta_product_level_23():
    from artin.check_newforms import check_dihedral, class_group, characters_of, theta_coefficients
    forms, idx, tab, e = class_group(-23)
    psi = [c for c in characters_of(tab, e) if any(abs(c[i] - 1) > 1e-9 for i in c)][0]
    a = [round(x.real) for x in theta_coefficients(-23, psi, idx, 30)]
    assert a[1:14] == [1, -1, -1, 0, 0, 1, 0, 1, 0, 0, 0, 0, -1]      # eta(z) eta(23z)
    r = check_dihedral([-1, -1, 0, 1], [[2, 3, 1], [2, 1, 3]], -23, X=120, log=lambda m: None)
    assert r["all_ok"] and r["characters"][0]["n_mismatches"] == 0