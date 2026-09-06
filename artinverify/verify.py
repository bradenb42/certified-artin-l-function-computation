"""The checker.

Reads a certificate and walks its claims in order, recomputing each from the
recorded data: the character table, the ramified-prime witnesses, the
precision consultations, the Frobenius classes, the local groups, the
filtrations and conductor exponents, the Euler factors, the archimedean data,
the root numbers and the functional-equation defects.  The first claim that
fails is reported and the certificate is rejected.
"""
from __future__ import annotations
import json, math, sys, cmath
from fractions import Fraction
from math import gcd

from artin.perm import PermGroup, from_json, mul, inverse, identity, cycle_type, cycles
from.arith import Z, root_of_unity, factorization_type, pgcd, ptrim, pmul, pdivmod, pmod, ppowmod, is_prime, sturm_count, disc_mod_p_zero, phi
from.fields import Order, val

class Reject(Exception):
    pass

def need(cond, msg):
    if not cond:
        raise Reject(msg)

# ------------------------------------------------------------------ helpers

def resultant_disc(f):
    """disc of monic integer f by the Sylvester determinant over Q."""
    n = len(f) - 1
    df = [i * c for i, c in enumerate(f)][1:]
    N = 2 * n - 1
    M = []
    for i in range(n - 1):
        row = [Fraction(0)] * N
        for j, c in enumerate(reversed(f)):
            row[i + j] = Fraction(c)
        M.append(row)
    for i in range(n):
        row = [Fraction(0)] * N
        for j, c in enumerate(reversed(df)):
            row[i + j] = Fraction(c)
        M.append(row)
    # determinant
    A = M; d = Fraction(1)
    for c in range(N):
        piv = next((r for r in range(c, N) if A[r][c] != 0), None)
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]; d = -d
        d *= A[c][c]
        for r in range(c + 1, N):
            if A[r][c] != 0:
                fct = A[r][c] / A[c][c]
                A[r] = [x - fct * y for x, y in zip(A[r], A[c])]
    sign = -1 if (n * (n - 1) // 2) % 2 else 1
    return int(sign * d)

def is_irreducible_mod(g, p):
    g = ptrim([c % p for c in g]); d = len(g) - 1
    if d <= 1:
        return d == 1
    x = [0, 1]; xp = x
    for i in range(1, d // 2 + 1):
        xp = ppowmod(xp, p, g, p)
        diff = ptrim([(a - b) % p for a, b in zip(xp + [0] * 2, x + [0] * len(xp))])
        if len(pgcd(g, diff, p)) > 1:
            return False
    return True

def _imul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out

def coset_action(G, H):
    Hel = H.elements()
    reps, index = [], {}
    for x in G.elements():
        if x in index:
            continue
        c = len(reps); reps.append(x)
        for h in Hel:
            index[mul(x, h)] = c
    return lambda g: tuple(index[mul(g, x)] for x in reps)

# ------------------------------------------------------------------ the verifier

class Verifier:
    def __init__(self, cert):
        self.C = cert
        S = cert["sections"]
        self.f = [int(c) for c in cert["input"]["f"]]
        self.n = len(self.f) - 1
        g = S["group"]
        self.G = PermGroup([from_json(x) for x in g["generators"]], n=self.n)
        self.cl = self.G.classes()
        self.report = {}

    def run(self):
        steps = [("V0_group_table", self.v0), ("V0f_ramified", self.v0f), ("V0e_precision", self.v0e),
                 ("V6_classes", self.v6), ("V2_local", self.v2), ("V3_V4_filtration_conductors", self.v34),
                 ("V5_euler", self.v5), ("V7a_archimedean", self.v7a), ("V7_rootnumbers", self.v7),
                 ("V9_functional_equation", self.v9)]
        for name, fn in steps:
            try:
                fn()
                self.report[name] = "ok"
            except Reject as e:
                self.report[name] = f"REJECT: {e}"
                return False, self.report
        return True, self.report

    # ---- V0: group, classes, table certificate
    def v0(self):
        S = self.C["sections"]; g = S["group"]; T = S["chartable"]
        need(g["order"] == self.G.order(), "group order")
        C = g["classes"]
        need(C["r"] == self.cl.r, "number of classes")
        reps = [from_json(x) for x in C["representatives"]]
        # recorded representatives must be in the recorded classes with recorded sizes: identify by lookup
        self.cls_map = []   # recorded class index -> our class index
        for k, rep in enumerate(reps):
            need(self.G.contains(rep), "class representative not in G")
            j = self.cl.class_of(rep)
            need(self.cl.sizes[j] == C["sizes"][k] and self.cl.orders[j] == C["orders"][k], "class size/order")
            self.cls_map.append(j)
        need(sorted(self.cls_map) == list(range(self.cl.r)), "class representatives not a system")
        self.e = C["exponent"]
        need(self.e == self.cl.exponent, "exponent")
        # the table, in the recorded class ordering
        self.r = T["r"]
        self.deg = T["degrees"]
        self.vals = [[Z.from_json(v) for v in row] for row in T["values"]]
        # class matrices in the recorded ordering
        cert = T["certificate"]
        Sset = [i - 1 for i in cert["S"]]
        omegas = [[Z.from_json(w) for w in om] for om in cert["omegas"]]
        need(len(omegas) == self.r == self.cl.r, "table size")
        for nu in range(self.r):
            for k in range(self.r):
                need(omegas[nu][k] == self.vals[nu][k] * C["sizes"][k] / self.deg[nu], "central character vs table")
        for i in Sset:
            M = self.class_matrix(i, reps)
            for w in omegas:
                for j in range(self.r):
                    lhs = Z.const(self.e, 0)
                    for k in range(self.r):
                        if M[j][k]:
                            lhs = lhs + w[k] * M[j][k]
                    need(lhs == w[j] * w[i], f"eigenvector identity M_{i+1}")
        need(len({tuple(w[i].c for i in Sset) for w in omegas}) == self.r, "tuples not distinct")
        for nu in range(self.r):
            s = Z.const(self.e, 0)
            for k in range(self.r):
                s = s + self.vals[nu][k] * self.vals[nu][k].conj() * C["sizes"][k]
            need(s == self.G.order(), "orthogonality")
        need(all(v == 1 for v in self.vals[0]), "row 1 not trivial")
        self.reps = reps
        # power maps in recorded ordering
        self.power_map = C["power_maps"]
        for k, rep in enumerate(reps):
            x = identity(self.n)
            for t in range(C["orders"][k]):
                need(self.rec_class(x) == self.power_map[k][t], "power map")
                x = mul(x, rep)
        self.sizes = C["sizes"]; self.orders = C["orders"]

    def rec_class(self, g):
        return self.cls_map.index(self.cl.class_of(g))

    def class_matrix(self, i, reps):
        r = self.r
        M = [[0] * r for _ in range(r)]
        members = [g for g in self.G.elements() if self.rec_class(g) == i]
        for x in members:
            xi = inverse(x)
            for k, gk in enumerate(reps):
                M[self.rec_class(mul(xi, gk))][k] += 1
        return M

    # ---- V0(f): ramified witnesses
    def v0f(self):
        R = self.C["sections"]["ramified"]
        factors = R["factors"]
        prod = [1]
        # factors multiply to f (over Z)
        for fi in factors:
            new = [0] * (len(prod) + len(fi) - 1)
            for i, a in enumerate(prod):
                for j, b in enumerate(fi):
                    new[i + j] += a * b
            prod = new
        need(prod == self.f, "factors do not multiply to f")
        discs = [resultant_disc(fi) for fi in factors]
        for i, d in enumerate(discs):
            need(d == R["disc_factors"][i], f"disc of factor {i}")
            need(d != 0, "factor not separable")
            fac = {int(k): v for k, v in R["disc_factorizations"][i].items()}
            need(all(is_prime(p) for p in fac), "non-prime in factorization")
            pr = 1
            for p, v in fac.items():
                pr *= p ** v
            need(pr == abs(d), f"factorization of disc {i}")
        cand = sorted({int(p) for fa in R["disc_factorizations"] for p in fa})
        need(cand == R["candidates"], "candidate primes")
        ram = []
        for ell_s, rec in R["per_prime"].items():
            ell = int(ell_s)
            anyram = False
            for i_s, w in rec["factors"].items():
                i = int(i_s); fi = factors[i]
                v = val(abs(discs[i]), ell)
                need(v == w["v_disc_f"], "v_ell(disc f_i)")
                facs = [(g["factor"], g["exponent"]) for g in w["factorization_mod_ell"]]
                pr = [1]
                for g, e in facs:
                    need(is_irreducible_mod(g, ell), "recorded factor not irreducible")
                    for _ in range(e):
                        pr = pmul(pr, [c % ell for c in g], ell)
                need(pr == ptrim([c % ell for c in fi]), "factorization mod ell")
                if w["method"] == "odd valuation":
                    need(v % 2 == 1 and w["ramified"], "odd valuation verdict")
                    anyram = True
                    continue
                if w["method"] == "Dedekind criterion":
                    # h = (f - prod g~^e)/ell with symmetric lifts; gcd(h, g_j) = 1 for e_j >= 2
                    P = [1]
                    for g, e in facs:
                        gl = [c if c <= ell // 2 else c - ell for c in g]
                        for _ in range(e):
                            new = [0] * (len(P) + len(gl) - 1)
                            for a, x in enumerate(P):
                                for b, y in enumerate(gl):
                                    new[a + b] += x * y
                            P = new
                    diff = [a - b for a, b in zip(fi, P + [0] * (len(fi) - len(P)))]
                    need(all(c % ell == 0 for c in diff), "Dedekind h not integral")
                    h = ptrim([(c // ell) % ell for c in diff])
                    for g, e in facs:
                        if e >= 2:
                            gg = pgcd(h, g, ell) if h else g
                            need(len(gg) == 1, "Z[alpha] not ell-maximal (Dedekind)")
                    ramified = any(e >= 2 for g, e in facs)
                    need(ramified == w["ramified"] and w["v_disc_O"] == v, "Dedekind verdict")
                    need(sorted([(e, len(g) - 1) for g, e in facs]) == sorted(tuple(x) for x in w["residue_decomposition"]), "Dedekind decomposition")
                    anyram = anyram or ramified
                    continue
                need(w["method"] == "Round 2", "unknown method")
                B = [[Fraction(x) for x in row] for row in w["ell_maximal_order_basis"]]
                try:
                    O = Order(fi, ell, B)
                except ValueError as ex:
                    raise Reject(f"recorded order: {ex}")
                need(O.index_val == w["index_valuation"] and v - 2 * O.index_val == w["v_disc_O"], "index valuation")
                need(O.is_maximal(), "recorded order not ell-maximal")
                C, nd = O.radical_lattice()
                need((nd > 0) == w["ramified"], "radical vs verdict")
                dec = sorted(O.residue_decomposition())
                need(dec == sorted(tuple(x) for x in w["residue_decomposition"]), "residue decomposition")
                anyram = anyram or w["ramified"]
            need(anyram == rec["ramified"], "aggregated verdict")
            if anyram:
                ram.append(ell)
        need(ram == R["ramified"], "ramified list")
        self.ramified = R

    # ---- V0(e): precision policy
    def v0e(self):
        P = self.C["sections"]["precision"]
        prefix = [1]
        for pr in P["pairs"]:
            B = pr["norm1"] * P["root_bound_R"] ** pr["degree"]
            m = pr["m"]
            mu = 1 if m <= 1 else m * (m - 1)
            prefix.append(max(prefix[-1], (2 * B + 2) ** (mu if P["mode"] == "conservative" else m)))
        need(prefix[-1].bit_length() - 1 == P["log2_Mstar"], "M*")
        Delta = abs(resultant_disc(self.f))
        for entry in P["consultation_log"]:
            ell, k = entry["ell"], entry["k"]
            Ms = prefix[entry["pairs"]]
            nd, y = 0, ell
            while y <= Ms:
                y *= ell; nd += 1
            need(k >= max(nd + 1, val(Delta, ell) + 1), f"consultation at {ell}")

    # ---- V6: class records
    def v6(self):
        C = self.C["sections"]["classes"]
        subs = {}
        for t in C["family"]:
            lab = t["subgroup"]
            H = self.subgroup_from_label(lab)
            act = coset_action(self.G, H)
            sig = {k: cycle_type(act(rep)) for k, rep in enumerate(self.reps)}
            subs[str(lab)] = (t, sig)
        Delta = abs(resultant_disc(self.f))
        self.class_at = {}
        for rec in C["primes"]:
            ell, k = rec["ell"], rec["class"] - 1
            need(Delta % ell, "class recorded at a ramified prime")
            ct = factorization_type(self.f, ell)
            need(ct == rec["block"] == list(cycle_type(self.reps[k])), f"block at {ell}")
            for t, sig in subs.values():
                R = t["resolvent"]
                if disc_mod_p_zero(R, ell):
                    continue
                need(tuple(factorization_type(R, ell)) == sig[k], f"resolvent type at {ell} vs class")
            self.class_at[ell] = k

    def subgroup_from_label(self, lab):
        n = self.n
        kind = lab["kind"]
        if kind == "cyclic":
            return PermGroup([self.reps[lab["class"] - 1]], n=n)
        if kind == "point stabilizer":
            i = lab["point"] - 1
            return PermGroup([g for g in self.G.elements() if g[i] == i] or [identity(n)], n=n)
        if kind == "2-set stabilizer":
            S = {lab["set"][0] - 1, lab["set"][1] - 1}
            return PermGroup([g for g in self.G.elements() if {g[i] for i in S} == S] or [identity(n)], n=n)
        if kind == "base stabilizer":
            pts = [b - 1 for b in lab["prefix"]]
            return PermGroup([g for g in self.G.elements() if all(g[i] == i for i in pts)] or [identity(n)], n=n)
        raise Reject("unknown subgroup label")

    # ---- V2: local data
    def v2(self):
        L = self.C["sections"]["local"]
        self.local = {}
        for ell_s, r in L["ramified"].items():
            if r["status"] != "ok":
                continue
            Dg = [from_json(x) for x in r["D_generators"]]; Ig = [from_json(x) for x in r["I_generators"]]
            need(all(self.G.contains(x) for x in Dg + Ig), "D not in G")
            D = PermGroup(Dg, n=self.n); I = PermGroup(Ig, n=self.n)
            need(D.order() == r["|D|"] == r["e"] * r["f"] and I.order() == r["e"], "local orders")
            for g in Dg:
                for h in Ig:
                    need(I.contains(mul(g, mul(h, inverse(g)))), "I not normal")
            Fr = from_json(r["frobenius_rep"])
            need(D.contains(Fr), "Frobenius not in D")
            need(self.G.contains(mul(inverse(from_json(r["matching"])), from_json(r["matching2"]))), "matchings")
            need(r["inertia_orbits_match_factors"] and r["decomposition_orbits_match_factors"], "orbit checks")
            self.local[ell_s] = (r, D, I, Fr)

    # ---- V3/V4: filtration and conductors
    def v34(self):
        CJ = self.C["sections"]["conductors"]
        self.f_ell = {}
        for ell_s, rec in CJ["primes"].items():
            r, D, I, Fr = self.local[ell_s]
            ell = int(ell_s)
            rho = from_json(r["matching"])
            fil = rec["filtration"]
            # i_L is recorded on the local permutations; transport to the global numbering by rho
            i_L = {}
            for k, v in fil["i_L"].items():
                s_loc = tuple(int(x) - 1 for x in k.strip("[]").split(","))
                i_L[mul(inverse(rho), mul(s_loc, rho))] = v
            e = r["e"]
            need(len(i_L) == e - 1, "i_L not defined on I \\ 1")
            for s in i_L:
                need(I.contains(s) and s != identity(self.n), "i_L key not in I")
            delta = sum(i_L.values())
            need(delta == fil["delta"], "delta")
            P = [s for s, v in i_L.items() if v >= 2]
            m = len(P) + 1
            while m % ell == 0: m //= ell
            need(m == 1 and (e // (len(P) + 1)) % ell, "D_1 not the ell-Sylow")
            # D_i normal in D and elementary abelian quotients, Hasse-Arf on the whole chain when I abelian
            levels = []
            u = 0
            while True:
                Du = [s for s in I.elements() if s == identity(self.n) or i_L[s] >= u + 1]
                levels.append(Du)
                if len(Du) == 1: break
                u += 1
            for Du in levels:
                Sset = set(Du)
                need(all(mul(g, mul(s, inverse(g))) in Sset for g in D.elements() for s in Du), "D_i not normal")
            N = {}
            for s, v in i_L.items():
                k = self.rec_class(s)
                N[k] = N.get(k, 0) + v
            tot = 0
            fl = {}
            for c in rec["conductor_exponents"]:
                nu = c["chi"] - 1
                v = self.vals[nu][0] * delta
                for k, Nk in N.items():
                    v = v - self.vals[nu][k] * Nk
                q = v.rational() / e
                need(q.denominator == 1 and int(q) == c["f_ell"], f"f_ell(chi_{nu+1}) at {ell}")
                fl[nu] = int(q)
                tot += self.deg[nu] * int(q)
            vdN = (self.G.order() // (e * r["f"])) * r["f"] * delta
            need(tot == vdN == rec["identities"]["v_ell_d_N"], "conductor identity")
            # the per-factor identity (roots of factor j: from the ramified witness via v_disc_O)
            for fa in rec["identities"]["factors"]:
                j = fa["factor"]
                w = self.ramified["per_prime"].get(ell_s, {"factors": {}})["factors"].get(str(j))
                want = w["v_disc_O"] if w else 0
                need(fa["v_ell(disc O)"] == want and fa["f_ell(pi)"] == want, "conductor identity")
            self.f_ell[ell] = fl
        for c in CJ["conductors"]:
            nu = c["chi"] - 1
            pc = 1
            for ell_s, ex in c["exponents"].items():
                need(self.f_ell[int(ell_s)][nu] == ex, "conductor exponent aggregation")
                pc *= int(ell_s) ** ex
            need(pc == c["partial_conductor"], "partial conductor")

    # ---- V5: Euler factors by route (C)
    def v5(self):
        EJ = self.C["sections"]["euler"]
        self.euler = {}
        for ell_s, rec in EJ["ramified"].items():
            r, D, I, Fr = self.local[ell_s]
            ell = int(ell_s); fdeg = r["f"]
            Del, Iel = D.elements(), I.elements()
            k_of = {}
            x = identity(self.n)
            for k in range(fdeg):
                for i in Iel:
                    k_of[mul(x, i)] = k
                x = mul(x, Fr)
            need(set(k_of) == set(Del), "coset decomposition")
            self.euler[ell] = {}
            Ps = {}
            for chi_s, c in rec["characters"].items():
                nu = int(chi_s) - 1
                E = self.e * fdeg // gcd(self.e, fdeg)
                P = [Z.const(E, 1)]
                dimVI = 0
                for j in range(fdeg):
                    s = Z.const(E, 0)
                    for sg in Del:
                        s = s + self.vals[nu][self.rec_class(sg)].lift(E) * root_of_unity(E, fdeg, (-j * k_of[sg]) % fdeg)
                    q = s.rational() / len(Del)
                    need(q.denominator == 1 and q >= 0, "multiplicity")
                    for _ in range(int(q)):
                        P = self.cmul(P, [Z.const(E, 1), -root_of_unity(E, fdeg, j)])
                    dimVI += int(q)
                Prec = [Z.from_json(x) for x in c["P"]]
                need(len(Prec) == len(P) and all(a == b for a, b in zip(P, Prec)), f"Euler factor chi_{nu+1} at {ell}")
                need(len(P) - 1 == dimVI == self.deg[nu] - c["f_ell"] + c["swan"] and c["f_ell"] == self.f_ell[ell][nu], "dimension check")
                self.euler[ell][nu] = P
                Ps[nu] = P
            # product identity over Z via Galois orbits of rows
            orbits = self.C["sections"]["chartable"]["galois_orbits"]
            prod = [1]
            for orb in orbits:
                Q = [Z.const(Ps[orb[0]][0].e, 1)]
                for nu in orb:
                    Q = self.cmul(Q, Ps[nu])
                Qi = []
                for c in Q:
                    q = c.rational()
                    need(q.denominator == 1, "orbit product not integral")
                    Qi.append(int(q))
                for _ in range(self.deg[orb[0]]):
                    prod = _imul(prod, Qi)
            want = [1]
            for _ in range(self.G.order() // D.order()):
                want = _imul(want, [1] + [0] * (fdeg - 1) + [-1])
            need(prod == want, "product identity")
        for ell_s, rec in EJ["determinant_checks"].items():
            for chi_s, c in rec["det"].items():
                need(c["f_ell(det)"] <= c["f_ell(chi)"], "conductor of the determinant character")
        for ell_s, idents in EJ["zeta_identities"].items():
            need(all(x["ok"] for x in idents), "the Euler-factor identity system verdicts")

    @staticmethod
    def cmul(a, b):
        E = max(a[0].e, b[0].e)
        out = [Z.const(E, 0) for _ in range(len(a) + len(b) - 1)]
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] = out[i + j] + x * y
        return out

    # ---- V7(a): archimedean
    def v7a(self):
        A = self.C["sections"]["archimedean"]
        r = sturm_count(self.f)
        need(r == A["real_roots"], "real root count")
        k = A["class_of_c"] - 1
        ct = tuple(sorted([2] * ((self.n - r) // 2) + [1] * r, reverse=True))
        need(cycle_type(self.reps[k]) == ct, "cycle type of c")
        self.ab = {}
        for d in A["characters"]:
            nu = d["chi"] - 1
            chic = self.vals[nu][k].rational()
            need(d["a"] + d["b"] == self.deg[nu] and d["a"] - d["b"] == chic, "a, b")
            # parity via eigenvalue multiplicities of c
            o = self.orders[k]
            s = 0
            for j in range(o):
                m = Z.const(self.e, 0)
                for t in range(o):
                    m = m + self.vals[nu][self.power_map[k][t]] * root_of_unity(self.e, o, (-j * t) % o)
                mj = m.rational() / o
                need(mj.denominator == 1, "multiplicity")
                s += j * int(mj)
            detc = root_of_unity(self.e, o, s % o).rational() if (s % o) in (0, o // 2 if o % 2 == 0 else -1) else None
            need(detc == (-1) ** d["b"], "parity")
            self.ab[nu] = (d["a"], d["b"])

    # ---- V7: root numbers (recorded local values: modulus, mu_4 identification, W, shortcuts;
    # the local factors themselves are not recomputed here; only the recorded values,
    # their identification as roots of unity, the product and the self-dual predictions
    def v7(self):
        RJ = self.C["sections"]["rootnumbers"]
        T = self.C["sections"]["chartable"]
        self.W = {}
        for c in RJ["characters"]:
            if c.get("W_complex") is None:
                continue
            nu = c["chi"] - 1
            k4 = (-c["b"]) % 4
            mu4 = True
            Wc = cmath.exp(-1j * math.pi / 2 * c["b"])
            for ell_s, loc in c["local"].items():
                val = complex(*loc["complex"])
                need(abs(abs(val) - 1) < 1e-7, f"modulus at {ell_s}")
                if loc["i_exponent"] is None:
                    mu4 = False
                    need(not loc["in_mu4"], "mu_4 flag")
                else:
                    need(abs(val - 1j ** loc["i_exponent"]) < 1e-6, f"i-exponent at {ell_s}")
                    k4 = (k4 + loc["i_exponent"]) % 4
                Wc *= val
            need(abs(Wc - complex(*c["W_complex"])) < 1e-7, "W value")
            real, fs = T["is_real"][nu], T["frobenius_schur"][nu]
            if real and fs == 1:
                need(mu4 and k4 == 0, f"orthogonal chi_{nu+1}: W != 1")
            if real and fs == -1:
                need(mu4 and k4 % 2 == 0, "symplectic")
            self.W[nu] = Wc

    def sqrt_prime(self, ell, M):
        if ell == 2:
            need(M % 8 == 0, "field for sqrt 2")
            return Z.zeta(M, M // 8) + Z.zeta(M, 7 * M // 8)
        need(M % (4 * ell) == 0, "field for sqrt ell")
        g = Z.const(M, 0)
        for a in range(1, ell):
            leg = 1 if pow(a, (ell - 1) // 2, ell) == 1 else -1
            g = g + Z.zeta(M, (M // ell) * a) * leg
        return g if ell % 4 == 1 else g * (-Z.zeta(M, M // 4))

    # ---- V9: functional equation from the recorded coefficients, own kernel and bound
    def v9(self):
        S = self.C["sections"]
        if "analytic" not in S or "coefficients" not in S:
            return
        AN, CO = S["analytic"], S["coefficients"]
        X = CO.get("X")
        if X is None or not CO.get("coefficients"):
            return
        # coefficients at primes must match the recorded classes / Euler factors
        for nu in range(self.r):
            co = CO["coefficients"][str(nu + 1)]
            for ell, k in self.class_at.items():
                if ell <= X:
                    need(abs(complex(*co[ell]) - self.vals[nu][k].to_complex()) < 1e-9, f"a_{ell}(chi_{nu+1}) vs class")
            for ell, Ps in self.euler.items():
                if ell <= X and nu in Ps:
                    a = -Ps[nu][1].to_complex() if len(Ps[nu]) > 1 else 0
                    need(abs(complex(*co[ell]) - a) < 1e-9, f"a_{ell}(chi_{nu+1}) vs Euler factor")
        import mpmath
        mpmath.mp.dps = 25
        kern = {}
        def kernel(a, b):
            if (a, b) in kern:
                return kern[(a, b)]
            h = 0.05; T = max(30.0, 20.0 + 60.0 / (a + b)); nn = int(T / h)
            us = [k * h for k in range(-nn, nn + 1)]
            def gam(s):
                lg = a * (-(s / 2) * mpmath.log(mpmath.pi) + mpmath.loggamma(s / 2)) + b * (-((s + 1) / 2) * mpmath.log(mpmath.pi) + mpmath.loggamma((s + 1) / 2))
                return mpmath.exp(lg)
            gs = [complex(gam(mpmath.mpc(1, u))) for u in us]
            def g(x):
                lx = math.log(x)
                return (sum(gv * cmath.exp(-(1 + 1j * u) * lx) for u, gv in zip(us, gs)) * h / (2 * math.pi)).real
            kern[(a, b)] = g
            return g
        def Mstar(a, b):
            d = a + b
            return 2 ** (d / 2) / math.pi ** (1 + b) * (4 / (math.pi * d) + math.pi ** (b / 2) * math.gamma(b / 2 + 1) * (4 / (math.pi * d)) ** (b / 2 + 1))
        def Echi(d, a, b, f, t, X):
            alpha = math.pi * d / 2 * (t / math.sqrt(f)) ** (2 / d); kappa = d * (d - 1) / 4; y = alpha * X ** (2 / d)
            if y < max(d * (d - 3) / 4, 2 * kappa, 1):
                return None
            K = 2 ** ((b + d) / 2) * Mstar(a, b)
            G = float(mpmath.e1(y)) if kappa == 0 else float(mpmath.gammainc(kappa, y, mpmath.inf))
            return 2 ** (d - 1) * K * math.sqrt(f) / t * d / 2 * alpha ** (-kappa) * G
        conj_row = self.conj_rows()
        for c in AN["characters"]:
            if c.get("status") != "ok":
                continue
            nu = c["chi"] - 1
            a, b = c["gamma"]["a"], c["gamma"]["b"]; d = a + b
            need((a, b) == self.ab[nu], "gamma factor vs archimedean")
            f = c["conductor"]
            need(f == self.C["sections"]["conductors"]["conductors"][nu]["partial_conductor"], "conductor in FE test")
            g = kernel(a, b)
            co = [complex(*x) for x in CO["coefficients"][str(nu + 1)]]
            cob = [complex(*x) for x in CO["coefficients"][str(conj_row[nu] + 1)]]
            W = self.W[nu]
            need(abs(W - complex(*c["W"])) < 1e-9, "W in FE test")
            sf = math.sqrt(f)
            for tst in c["tests"]:
                t = tst["t"]
                lhs = sum(co[m] * g(m * t / sf) for m in range(1, X + 1))
                rhs = W * sum(cob[m] * g(m / t / sf) for m in range(1, X + 1)) / t
                D = abs(lhs - rhs)
                E1, E2 = Echi(d, a, b, f, t, X), Echi(d, a, b, f, 1 / t, X)
                if E1 is None or E2 is None:
                    continue
                ev = 1e-13 * (1 + abs(lhs) + abs(rhs)) * math.log(X + 2)
                need(D <= E1 + E2 + ev, f"FE defect chi_{nu+1} at t={t}: {D:.2e} > {E1+E2:.2e}")

    def conj_rows(self):
        out = []
        for nu in range(self.r):
            row = [v.conj() for v in self.vals[nu]]
            j = next((m for m in range(self.r) if all(a == b for a, b in zip(self.vals[m], row))), None)
            need(j is not None, "conjugate row missing")
            out.append(j)
        return out

def main(argv=None):
    path = (argv or sys.argv[1:])[0]
    with open(path) as fh:
        cert = json.load(fh)
    ok, rep = Verifier(cert).run()
    for k, v in rep.items():
        print(f"{k}: {v}")
    print("VERDICT:", "ACCEPT" if ok else "REJECT")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())