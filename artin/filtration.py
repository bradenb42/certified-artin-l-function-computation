"""Ramification filtrations and conductor exponents.

Computes the lower ramification filtration from a uniformiser of the local
field, the Herbrand function and the upper breaks, and the ramification
polygons of the compositum and of the individual factors.  Includes the
structural checks (normality, the wild inertia subgroup, elementary abelian
quotients, integrality of the upper breaks on abelian subquotients), the
enumeration of the filtrations consistent with the factor polygons, and the
Artin conductor exponents with the identities that tie their sum to the
discriminant valuation of the Galois closure.
"""
from __future__ import annotations
from fractions import Fraction
from itertools import product
from math import gcd

from .cyclo import Cyc
from .perm import PermGroup, mul, inverse, identity
from .local import LocalGalois, _orbits

class HardFailure(Exception):
    pass

# ------------------------------------------------------------------ automorphisms of K -> D

def automorphism_words(local: LocalGalois):
    """All elements of Gal(K/Q_ell) as words phi^a tau^b with their permutations of the roots."""
    K = local.K
    out = {}
    for a in range(local.F):
        for b in range(local.E):
            perm = local.phi
            perm = identity(len(local.roots))
            for _ in range(a):
                perm = mul(local.phi, perm)
            for _ in range(b):
                perm = mul(local.tau, perm)
            out.setdefault(perm, (a, b))
    return out

def apply_word(local: LocalGalois, word, y):
    a, b = word
    K = local.K
    for _ in range(a):
        y = K.frob(y)
    for _ in range(b):
        y = K.tau(y)
    return y

# ------------------------------------------------------------------ the filtration

class Filtration:
    """Lower filtration of D = Gal(L/Q_ell) computed in K (L = splitting field inside K)."""
    def __init__(self, local: LocalGalois, log=print):
        self.local = local
        K = local.K
        n = len(local.roots)
        self.words = automorphism_words(local)
        self.D_elems = list(self.words.keys())
        self.I_elems = local.I.elements()
        self.e, self.f = local.e, local.f
        self.eKL = local.E // self.e      # e(K/L): v_K = e(K/L) v_L on L
        # uniformizer of L: a polynomial in the roots with v_K = e(K/L)
        self.pi_L, self.pi_L_desc = self._uniformizer()
        # i_L(sigma) = v_L(sigma pi_L - pi_L), sigma in I \ 1
        self.i_L = {}
        for s in self.I_elems:
            if s == identity(n):
                continue
            w = self.words[s]
            d = K.sub(apply_word(local, w, self.pi_L), self.pi_L)
            vK = K.val(d)
            if vK >= K.prec - K.E:
                raise HardFailure("precision too low to read i_L")
            if vK % self.eKL:
                raise HardFailure("v_K(sigma pi_L - pi_L) not a multiple of e(K/L): pi_L not in L?")
            self.i_L[s] = vK // self.eKL
        self.delta = sum(self.i_L.values())   # v_L of the different of L/L^ur (Hilbert)
        self.lower = self._lower()
        self.P = self.lower[1] if len(self.lower) > 1 else [identity(n)]
        self.checks = self._structural_checks()
        self.herbrand, self.lower_breaks, self.upper_breaks = self._herbrand()
        log(f"filtration at {local.ell}: e = {self.e}, f = {self.f}, |P| = {len(self.P)}, i_L values {sorted(set(self.i_L.values()))}, delta = {self.delta}, breaks lower {self.lower_breaks} upper {self.upper_breaks}")

    def _uniformizer(self):
        K, local = self.local.K, self.local
        n = len(local.roots)
        target = self.eKL
        cands = []
        for i in range(n):
            cands.append((local.roots[i], f"beta_{i+1}"))
        for i in range(n):
            for j in range(i + 1, n):
                cands.append((K.sub(local.roots[i], local.roots[j]), f"beta_{i+1}-beta_{j+1}"))
        cands.append((K.from_int(local.ell), "ell"))
        # products and quotients by ell of the above, valuations reduced modulo e(K/L)*e... search by valuation
        base = [(y, d) for (y, d) in cands if 0 < K.val(y) < K.prec]
        base.sort(key=lambda t: K.val(t[0]))
        for y, d in base:
            if K.val(y) == target:
                return y, d
        # combine: y1^a * y2^b / ell^c with small exponents to reach valuation target
        vals = [(K.val(y), y, d) for y, d in base][:6]
        for (v1, y1, d1) in vals:
            for (v2, y2, d2) in vals:
                for a in range(0, 4):
                    for b in range(0, 4):
                        for c in range(0, 3):
                            if a * v1 + b * v2 - c * local.E == target and a + b > 0:
                                y = K.mul(K.power(y1, a), K.power(y2, b))
                                if c:
                                    y = K.div_pi(y, c * local.E)
                                if K.val(y) == target:
                                    return y, f"({d1})^{a} ({d2})^{b} / ell^{c}"
        raise HardFailure("no uniformizer of L found among the candidates")

    def _lower(self):
        """D_i = {sigma in I: i_L(sigma) >= i+1}, i = 0, 1, ..., until trivial."""
        n = len(self.local.roots)
        e1 = identity(n)
        levels = []
        i = 0
        while True:
            Di = [s for s in self.I_elems if s == e1 or self.i_L[s] >= i + 1]
            levels.append(Di)
            if len(Di) == 1:
                break
            i += 1
        return levels

    def _structural_checks(self):
        ell = self.local.ell
        n = len(self.local.roots)
        e1 = identity(n)
        D = self.D_elems
        checks = {}
        # normality of each D_i in D
        for i, Di in enumerate(self.lower):
            S = set(Di)
            if any(mul(g, mul(s, inverse(g))) not in S for g in D for s in Di):
                raise HardFailure(f"D_{i} not normal in D")
        # P = D_1 is the ell-Sylow subgroup of I: an ell-group of index prime to ell
        P = self.lower[1] if len(self.lower) > 1 else [e1]
        if len(P) & (len(P) - 1) == 0 and len(P) == 1:
            pass
        m = len(P)
        while m % ell == 0:
            m //= ell
        if m != 1 or (len(self.I_elems) // len(P)) % ell == 0:
            raise HardFailure("D_1 is not the ell-Sylow subgroup of the inertia group")
        checks["P_is_ell_Sylow"] = True
        # D_i / D_{i+1} elementary abelian ell-groups for i >= 1
        for i in range(1, len(self.lower) - 1):
            Di, Di1 = self.lower[i], set(self.lower[i + 1])
            q = len(Di) // len(Di1)
            mm = q
            while mm % ell == 0:
                mm //= ell
            if mm != 1:
                raise HardFailure(f"D_{i}/D_{i+1} not an ell-group")
            for s in Di:
                for t in Di:
                    comm = mul(mul(s, t), mul(inverse(s), inverse(t)))
                    if comm not in Di1 or mul(*[s] * 1) is None:
                        raise HardFailure(f"D_{i}/D_{i+1} not abelian")
                sp = identity(n)
                for _ in range(ell):
                    sp = mul(sp, s)
                if sp not in Di1:
                    raise HardFailure(f"D_{i}/D_{i+1} not elementary abelian")
        checks["elementary_abelian_quotients"] = True
        # I / P cyclic of order prime to ell (tame quotient): |I/P| prime to ell already; cyclicity
        checks["tame_quotient_order"] = len(self.I_elems) // len(P)
        return checks

    def _herbrand(self):
        """phi_L(u) = int_0^u dt/[D_0:D_t]; lower breaks (i with D_i != D_{i+1}) and upper breaks."""
        orders = [len(Di) for Di in self.lower]
        e = orders[0]
        lower_breaks = [i for i in range(len(orders) - 1) if orders[i] != orders[i + 1]]
        # phi at integers: phi(u) = (1/e) sum_{i=1}^{u} |D_i|
        def phi(u):
            return Fraction(sum(orders[i] if i < len(orders) else 1 for i in range(1, u + 1)), e)
        upper_breaks = [phi(i) for i in lower_breaks]
        return {"orders": orders, "phi_at_breaks": [(i, str(phi(i))) for i in lower_breaks]}, lower_breaks, [str(u) for u in upper_breaks]

    # -------------------------------------------------------------- Hasse-Arf on abelian subquotients
    def hasse_arf(self):
        """For every pair D_a >= D_b (b > a) with D_a/D_b abelian, the upper breaks of the
        quotient filtration must be integers.  Raises HardFailure."""
        n = len(self.local.roots)
        orders = self.herbrand["orders"]
        levels = self.lower
        # the quotient D_a/D_b has filtration (D_a/D_b)_i = D_{a+i} D_b / D_b for the lower numbering
        # of the field L^{D_b}/L^{D_a}; its Herbrand function uses orders |D_{a+i} D_b / D_b| = |D_{max(a+i,b)}| ... /|D_b|
        for b in range(1, len(levels)):
            for a in range(0, b):
                Da, Db = levels[a], set(levels[b])
                # abelian?
                if any(mul(mul(s, t), mul(inverse(s), inverse(t))) not in Db for s in Da for t in Da):
                    continue
                # ramification groups of Da/Db in the lower numbering of L^{Db}: (Da/Db)_u = D_{a+u}Db/Db  (upper numbering
                # of the quotient is that of D restricted; we use Herbrand of L^{Db}/L^{Da} via the quotient orders)
                qorders = [len(levels[min(a + u, len(levels) - 1)]) // len(Db) if a + u < len(levels) else 1 for u in range(len(levels) - a + 1)]
                e_q = qorders[0]
                breaks = [u for u in range(len(qorders) - 1) if qorders[u] != qorders[u + 1]]
                for u in breaks:
                    val = Fraction(sum(qorders[i] for i in range(1, u + 1)), e_q)
                    if val.denominator != 1:
                        raise HardFailure(f"Hasse-Arf fails on D_{a}/D_{b}: upper break {val} not integral")
        return True

    # -------------------------------------------------------------- polygons
    def polygons(self, factor_root_sets):
        """Compositum polygon (multiset of i_L on I\\1) and, for each factor j with root beta_j,
        the factor polygon {v_L(sigma beta_j - beta_j): sigma in I, sigma beta_j != beta_j},
        checked against  (>=, equality when beta_j generates the ring of integers of
        the factor field over O_{L^ur})."""
        K, local = self.local.K, self.local
        comp = sorted(self.i_L.values())
        facs = []
        for roots_idx in factor_root_sets:
            j = roots_idx[0]
            beta = local.roots[j]
            H = [s for s in self.I_elems if s[j] == j]
            poly = []
            lemma_ok = True
            for s in self.I_elems:
                if s[j] == j:
                    continue
                vK = K.val(K.sub(apply_word(local, self.words[s], beta), beta))
                if vK % self.eKL:
                    raise HardFailure("factor polygon valuation not in L")
                v = vK // self.eKL
                poly.append(v)
                coset_sum = sum(self.i_L[mul(s, h)] for h in H if mul(s, h) != identity(len(local.roots)))
                if v < coset_sum:
                    lemma_ok = False
            facs.append({"root": j + 1, "stabilizer_order": len(H), "polygon": sorted(poly), "lemma_3_1_ok": lemma_ok})
        return {"compositum": comp, "factors": facs}

    # -------------------------------------------------------------- conductor exponents
    def conductor_exponents(self, table, cl, rho_to_global):
        """f_ell(chi) = (1/|I|)(chi(1) delta - sum_k N_k chi(g_k)), N_k = sum_{sigma in (I\\1) cap C_k} i_L(sigma),
        with class labels in G through the matching rho (sigma_global = rho^{-1} sigma rho)."""
        n = len(self.local.roots)
        N = {}
        for s, i in self.i_L.items():
            sg = mul(inverse(rho_to_global), mul(s, rho_to_global))
            k = cl.class_of(sg)
            N[k] = N.get(k, 0) + i
        e = len(self.I_elems)
        out = []
        for nu in range(table.r):
            val = table.values[nu][0] * self.delta
            for k, Nk in N.items():
                val = val - table.values[nu][k] * Nk
            q = val.rational() / e if val.is_rational() else None
            if q is None or q.denominator != 1 or q < 0:
                raise HardFailure(f"conductor exponent of chi_{nu+1} at {self.local.ell} is not a nonnegative integer: {q}")
            fl = int(q)
            # Swan part: sum over i >= 1 of |D_i|/|D_0| codim V^{D_i}
            swan = 0
            for i, Di in enumerate(self.lower):
                if i == 0:
                    continue
                inv = Cyc.zero(table.e)
                for s in Di:
                    sg = mul(inverse(rho_to_global), mul(s, rho_to_global))
                    inv = inv + table.values[nu][cl.class_of(sg)]
                dimfix = inv.rational() / len(Di)
                assert dimfix.denominator == 1
                swan += Fraction(len(Di), e) * (table.degrees[nu] - int(dimfix))
            assert swan.denominator == 1
            out.append({"chi": nu + 1, "f_ell": fl, "swan": int(swan), "tame": fl - int(swan)})
        self.N_by_class = {k + 1: v for k, v in N.items()}
        return out

    def identities(self, table, cl, conds, G, rho_to_global, A2_records, factor_orbits):
        """The two identities, both hard failures: sum_chi chi(1) f_ell(chi) = v_ell(d_N) = [G:D] f delta,
        and, for each irreducible factor f_j of f, the conductor exponent of the permutation
        character of its root set equals the discriminant valuation of its ring of integers."""
        e = len(self.I_elems)
        lhs = sum(table.degrees[c["chi"] - 1] * c["f_ell"] for c in conds)
        vdN = (G.order() // (e * self.f)) * self.f * self.delta
        if lhs != vdN:
            raise HardFailure(f"the conductor identity fails at {self.local.ell}: sum chi(1) f_ell = {lhs}, [G:D] f delta = {vdN}")
        # permutation characters of the factor orbits: pi_j(sigma) = fixed points on the orbit
        res = {"sum_chi1_f": lhs, "v_ell_d_N": vdN, "factors": []}
        for j, orb in enumerate(factor_orbits):
            # f_ell(pi) = (1/|I|)(|orb| delta - sum_{sigma in I\1} i_L(sigma) fix_orb(sigma))
            s = len(orb) * self.delta
            for sg, i in self.i_L.items():
                s -= i * sum(1 for x in orb if sg[x] == x)
            assert s % e == 0
            f_pi = s // e
            want = A2_records.get(j)
            res["factors"].append({"factor": j, "f_ell(pi)": f_pi, "v_ell(disc O)": want})
            if want is not None and f_pi != want:
                raise HardFailure(f"the factor-conductor identity fails at {self.local.ell} for factor {j}: f_ell(pi_H) = {f_pi} but v_ell(disc O) = {want}")
        return res

    def to_json(self):
        n = len(self.local.roots)
        return {"e": self.e, "f": self.f, "uniformizer_of_L": self.pi_L_desc, "e(K/L)": self.eKL,
                "i_L": {str([i + 1 for i in s]): v for s, v in self.i_L.items()},
                "delta": self.delta, "lower_orders": self.herbrand["orders"],
                "lower_breaks": self.lower_breaks, "upper_breaks": self.upper_breaks,
                "checks": self.checks}

# ------------------------------------------------- filtrations consistent with the polygons, and resolvents that discriminate them

def enumerate_candidates(I_elems, D_elems, ell, factor_data, max_i, n):
    """All class functions i: I\\1 -> [1, max_i], constant on D-classes,
    defining a chain of normal subgroups with D_1 the ell-Sylow of I and elementary abelian
    quotients, whose coset sums reproduce the given factor polygons (equality form) and pass
    Hasse-Arf.  factor_data: list of (root index j, polygon sorted list).  Returns the list of
    candidate dicts."""
    e1 = identity(n)
    elems = [s for s in I_elems if s != e1]
    # D-classes on I \ 1
    classes, seen = [], set()
    for s in elems:
        if s in seen:
            continue
        cls = {mul(g, mul(s, inverse(g))) for g in D_elems}
        seen |= cls
        classes.append(sorted(cls))
    # elements of ell-power order must have i >= 2 (they lie in P), others exactly 1
    def order(s):
        k, x = 1, s
        while x != e1:
            x = mul(x, s); k += 1
        return k
    ranges = []
    for cls in classes:
        o = order(cls[0])
        m = o
        while m % ell == 0:
            m //= ell
        ranges.append(range(2, max_i + 1) if m == 1 and o > 1 else range(1, 2))
    out = []
    for choice in product(*ranges):
        i = {}
        for cls, v in zip(classes, choice):
            for s in cls:
                i[s] = v
        # factor polygons via coset sums
        ok = True
        for j, poly in factor_data:
            H = [s for s in I_elems if s[j] == j]
            got = sorted(sum(i[mul(s, h)] for h in H if mul(s, h) != e1) for s in I_elems if s[j] != j)
            if got != sorted(poly):
                ok = False
                break
        if not ok:
            continue
        # chain and Hasse-Arf via a lightweight filtration object
        try:
            _check_chain(i, I_elems, D_elems, ell, n)
        except HardFailure:
            continue
        out.append({str([x + 1 for x in cls[0]]): v for cls, v in zip(classes, choice)})
    return out

def _check_chain(i, I_elems, D_elems, ell, n):
    e1 = identity(n)
    levels = []
    u = 0
    while True:
        Du = [s for s in I_elems if s == e1 or i[s] >= u + 1]
        levels.append(Du)
        if len(Du) == 1:
            break
        u += 1
    for Du in levels:
        S = set(Du)
        if any(mul(g, mul(s, inverse(g))) not in S for g in D_elems for s in Du):
            raise HardFailure("normality")
    for u in range(1, len(levels) - 1):
        Du, Du1 = levels[u], set(levels[u + 1])
        for s in Du:
            sp = e1
            for _ in range(ell):
                sp = mul(sp, s)
            if sp not in Du1:
                raise HardFailure("elementary abelian")
            for t in Du:
                if mul(mul(s, t), mul(inverse(s), inverse(t))) not in Du1:
                    raise HardFailure("abelian")
    # Hasse-Arf on abelian subquotients
    orders = [len(L) for L in levels]
    for b in range(1, len(levels)):
        for a in range(b):
            Da, Db = levels[a], set(levels[b])
            if any(mul(mul(s, t), mul(inverse(s), inverse(t))) not in Db for s in Da for t in Da):
                continue
            qorders = [len(levels[min(a + w, len(levels) - 1)]) // len(Db) if a + w < len(levels) else 1 for w in range(len(levels) - a + 1)]
            for w in range(len(qorders) - 1):
                if qorders[w] != qorders[w + 1]:
                    val = Fraction(sum(qorders[t] for t in range(1, w + 1)), qorders[0])
                    if val.denominator != 1:
                        raise HardFailure("Hasse-Arf")

def discriminating_subgroups(filt: Filtration):
    """Subgroups S of D with S meet P = 1, as element lists, sorted by index."""
    D = filt.D_elems
    P = set(filt.P)
    n = len(filt.local.roots)
    e1 = identity(n)
    subs = {}
    # cyclic subgroups and the subgroups generated by pairs (enough for small D)
    gens_sets = [[e1]] + [[g] for g in D] + [[g, h] for g in D for h in D if g < h]
    for gs in gens_sets:
        H = PermGroup(gs, n=n) if any(g != e1 for g in gs) else PermGroup([], n=n)
        els = frozenset(H.elements())
        if els in subs:
            continue
        if len(els & P) == 1:
            subs[els] = H
    return sorted(subs.values(), key=lambda H: len(D) // H.order())

def discriminating_resolvent_check(filt: Filtration, S: PermGroup, monomial_base):
    """y_S = sum_{s in S} s.(monomial), J_S(sigma) = v_L(sigma y_S - y_S) for sigma in D \\ S; check
    's relation J_S(sigma) >= |S meet I| + i_L(sigma) - 1 on P \\ 1 and the coset-sum bound
    J_S(sigma) >= sum_{tau in sigma S, tau in I} i_L(tau) (equalities when y_S generates)."""
    K, local = filt.local.K, filt.local
    n = len(local.roots)
    e1 = identity(n)
    def mono(perm):
        t = K.one()
        for j, b in enumerate(monomial_base):
            t = K.mul(t, K.power(local.roots[perm[b]], j + 1))
        return t
    y = K.zero()
    for s in S.elements():
        y = K.add(y, mono(s))
    Sset = set(S.elements())
    SI = len(Sset & set(filt.I_elems))
    rec = {"index": len(filt.D_elems) // S.order(), "|S meet I|": SI, "J": {}, "relation_ok": True}
    for sg in filt.I_elems:
        if sg in Sset or sg == e1:
            continue
        d = K.sub(apply_word(local, filt.words[sg], y), y)
        vK = K.val(d)
        if vK >= K.prec - K.E:
            continue
        J = Fraction(vK, filt.eKL)
        rec["J"][str([x + 1 for x in sg])] = str(J)
        cs = sum(filt.i_L[mul(sg, s)] for s in Sset if mul(sg, s) in filt.i_L)
        if J < cs:
            rec["relation_ok"] = False
        if sg in set(filt.P) and J < SI + filt.i_L[sg] - 1:
            rec["relation_ok"] = False
    return rec