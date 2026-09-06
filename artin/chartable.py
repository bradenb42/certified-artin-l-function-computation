"""Character tables by the Dixon-Schneider method.

The table is computed modulo a prime p = 1 mod exp(G) with p > 2 sqrt|G| and
lifted to exact values in Q(zeta_e) through eigenvalue multiplicities.  The
result carries a certificate: the central characters omega_chi in Z[zeta_e]
together with the classes whose class matrices were used, so that a
verifier can confirm the table by checking eigenvector identities rather than
by recomputing it.
"""
from __future__ import annotations
import math
import random
from fractions import Fraction
from math import gcd

from .cyclo import Cyc, phi, divisors
from .fpoly import charpoly, roots_mod_p, nullspace, rref, mat_vec
from .perm import ConjugacyClasses

def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def choose_prime(e, order):
    """p = 1 mod e, p > 2 sqrt(|G|), and an element z of order e mod p."""
    bound = int(2 * math.isqrt(order)) + 2
    p = e * ((bound // e) + 1) + 1
    while not is_prime(p):
        p += e
    # element of order e
    fac = [q for q in range(2, e + 1) if e % q == 0 and is_prime(q)]
    while True:
        g = random.randrange(2, p)
        z = pow(g, (p - 1) // e, p)
        if all(pow(z, e // q, p) != 1 for q in fac):
            return p, z


class CharacterTable:
    def __init__(self, classes: ConjugacyClasses, seed=0):
        self.cl = classes
        self.r = classes.r
        self.e = classes.exponent
        self.order = classes.G.order()
        random.seed(seed)
        self.p, self.z = choose_prime(self.e, self.order)
        self._class_matrices = {}
        self._compute()
        self.certificate = self._certify()
        self._postprocess()

    # ------------------------------------------------------------ Dixon-Schneider
    def class_matrix(self, i):
        if i not in self._class_matrices:
            self._class_matrices[i] = self.cl.class_matrix(i)
        return self._class_matrices[i]

    def _compute(self):
        cl, r, p = self.cl, self.r, self.p
        # parts: list of (basis rows (v x r) in rref, set of used classes)
        parts = [([[1 if i == j else 0 for i in range(r)] for j in range(r)], set())]
        order_classes = sorted(range(1, r), key=lambda i: (cl.sizes[i], i))
        done = []
        used_total = set()
        while parts:
            B, used = parts.pop()
            v = len(B)
            if v == 1:
                done.append(B[0])
                continue
            i = next((i for i in order_classes if i not in used), None)
            if i is None:
                raise RuntimeError("could not split a joint eigenspace: class matrices exhausted")
            used = used | {i}
            used_total.add(i)
            M = [[x % p for x in row] for row in self.class_matrix(i)]
            # restrict: M B^T = B^T A; with B in rref, A = rows of (M B^T) at pivot columns
            Bt_cols = B  # each row of B is a basis vector
            MB = [mat_vec(M, b, p) for b in B]  # images of basis vectors, as vectors of length r
            R, piv = rref(B, p)
            # express each image in the basis: coordinates = image at pivot positions (B is rref)
            A = [[MB[k][c] for c in piv] for k in range(v)]  # A[k][j]: coefficient of basis j in image of basis k
            # eigenvalues of A^T (action on coordinates); charpoly same
            cp = charpoly(A, p)
            lams = roots_mod_p(cp, p)
            if sum(1 for _ in lams) == 0:
                raise RuntimeError("no eigenvalue found")
            total = 0
            for lam in lams:
                # eigenvectors: coordinate vectors y with y^T A = lam y^T  (since image of sum y_k b_k = sum y_k A[k][:] b)
                At = [[(A[k][j] - (lam if j == k else 0)) % p for k in range(v)] for j in range(v)]
                ys = nullspace(At, p, ncols=v)
                if not ys:
                    continue
                vecs = [[sum(y[k] * B[k][c] for k in range(v)) % p for c in range(r)] for y in ys]
                Rv, _ = rref(vecs, p)
                total += len(Rv)
                parts.append((Rv, set(used)))
            if total != v:
                raise RuntimeError("eigenspace dimensions do not sum up; bad prime?")
        self.S = sorted(used_total)
        # normalize and lift
        rows = []
        for w in done:
            inv = pow(w[0], -1, p)
            w = [(x * inv) % p for x in w]
            rows.append(w)
        self._reduced_omegas = rows
        self.degrees = []
        self.values = []  # list of lists of Cyc
        for w in rows:
            d = self._degree(w)
            self.degrees.append(d)
            self.values.append(self._lift(w, d))
        # sort: trivial first, then by degree
        idx = sorted(range(self.r), key=lambda k: (self.degrees[k], [str(v.c) for v in self.values[k]]))
        # trivial character: all values 1
        triv = next(k for k in idx if all(v == 1 for v in self.values[k]))
        idx.remove(triv)
        idx = [triv] + idx
        self.degrees = [self.degrees[k] for k in idx]
        self.values = [self.values[k] for k in idx]
        self._reduced_omegas = [self._reduced_omegas[k] for k in idx]

    def _degree(self, w):
        cl, p = self.cl, self.p
        s = 0
        for k in range(self.r):
            s = (s + w[k] * w[cl.inverse_class[k]] * pow(cl.sizes[k], -1, p)) % p
        target = (self.order * pow(s, -1, p)) % p
        root = None
        for d in range(1, math.isqrt(self.order) + 1):
            if d * d % p == target:
                root = d
                break
        if root is None:
            raise RuntimeError("degree not recovered")
        return root

    def _lift(self, w, d):
        cl, p, e, z = self.cl, self.p, self.e, self.z
        # theta(chi(g_k)) = d * w_k / |C_k|
        th = [(d * w[k] * pow(cl.sizes[k], -1, p)) % p for k in range(self.r)]
        vals = []
        for k in range(self.r):
            o = cl.orders[k]
            zo = pow(z, e // o, p)
            zo_inv = pow(zo, -1, p)
            oinv = pow(o, -1, p)
            coeffs = [0] * e
            total = 0
            for j in range(o):
                s = 0
                for t in range(o):
                    s = (s + th[cl.power_map(k, t)] * pow(zo_inv, (j * t) % o, p)) % p
                m = (s * oinv) % p
                if m > d:
                    raise RuntimeError("multiplicity out of range; lifting failed")
                total += m
                coeffs[(e // o) * j] += m
            if total != d:
                raise RuntimeError("multiplicities do not sum to the degree")
            vals.append(Cyc(e, coeffs))
        return vals

    # ------------------------------------------------------------ certificate
    def omega(self, k_chi):
        """Central character omega_chi(K_k) = |C_k| chi(g_k)/chi(1), exact in Z[zeta_e]."""
        d = self.degrees[k_chi]
        return [(self.values[k_chi][k] * self.cl.sizes[k]).divexact(d) for k in range(self.r)]

    def _certify(self):
        r, e = self.r, self.e
        omegas = [self.omega(nu) for nu in range(r)]
        ok = True
        problems = []
        tuples = []
        for nu, w in enumerate(omegas):
            if not (w[0] == 1):
                ok = False
                problems.append(f"omega_{nu}(K_1) != 1")
            for i in self.S:
                M = self.class_matrix(i)
                lhs = [sum((w[k] * M[j][k] for k in range(r) if M[j][k]), Cyc.zero(e)) for j in range(r)]
                rhs = [w[j] * w[i] for j in range(r)]
                if any(not (a == b) for a, b in zip(lhs, rhs)):
                    ok = False
                    problems.append(f"M_{i} omega_{nu} != omega_{nu}(K_{i}) omega_{nu}")
            tuples.append(tuple(w[i] for i in self.S))
        if len(set(tuples)) != r:
            ok = False
            problems.append("joint eigenvalue tuples not pairwise distinct")
        # orthogonality (independent check)
        for nu in range(r):
            s = Cyc.zero(e)
            for k in range(r):
                s = s + self.values[nu][k] * self.values[nu][k].conj() * self.cl.sizes[k]
            if not (s == self.order):
                ok = False
                problems.append(f"row {nu}: sum |C_k||chi(g_k)|^2 != |G|")
        return {"S": [i + 1 for i in self.S], "verified": ok, "problems": problems,
                "omegas": [[w.to_json() for w in om] for om in omegas]}

    # ------------------------------------------------------------ derived data
    def _postprocess(self):
        e, r = self.e, self.r
        units = [t for t in range(1, e + 1) if gcd(t, e) == 1]
        # Galois action on rows: chi^{sigma_t}(g_k) = chi(g_k^t)
        row_index = {tuple(v.c for v in row): i for i, row in enumerate(self.values)}
        self.galois_action = {}
        for t in units:
            perm = []
            for nu in range(r):
                img = tuple(self.values[nu][self.cl.power_map(k, t)].c for k in range(r))
                perm.append(row_index[img])
            self.galois_action[t] = perm
        seen = set()
        self.orbits = []
        for nu in range(r):
            if nu in seen:
                continue
            orb = sorted({self.galois_action[t][nu] for t in units})
            seen.update(orb)
            self.orbits.append(orb)
        self.stabilizers = [[t for t in units if self.galois_action[t][nu] == nu] for nu in range(r)]
        self.conductors = []
        for nu in range(r):
            c = min(c for c in divisors(e)
                    if all(self.galois_action[t][nu] == nu for t in units if t % c == 1 % c))
            self.conductors.append(c)
        self.field_degrees = [len(units) // len(st) for st in self.stabilizers]
        self.is_real = [self.galois_action[e - 1][nu] == nu for nu in range(r)]
        # Frobenius-Schur indicator nu_2 = |G|^{-1} sum chi(g^2)
        self.indicator = []
        for nu in range(r):
            s = Cyc.zero(e)
            for k in range(r):
                s = s + self.values[nu][self.cl.power_map(k, 2)] * self.cl.sizes[k]
            q = s.rational() / self.order
            assert q.denominator == 1 and int(q) in (-1, 0, 1)
            self.indicator.append(int(q))

    def inner(self, a, b):
        """<a, b> for class functions given as lists of Cyc (same length r)."""
        e = self.e
        s = Cyc.zero(e)
        for k in range(self.r):
            s = s + a[k] * b[k].conj() * self.cl.sizes[k]
        q = s.rational() / self.order
        assert q.denominator == 1
        return int(q)

    def to_json(self):
        return {"r": self.r, "exponent": self.e, "group_order": self.order,
                "dixon_prime": self.p, "z": self.z,
                "degrees": self.degrees,
                "values": [[v.to_json() for v in row] for row in self.values],
                "galois_orbits": self.orbits,
                "field_conductors": self.conductors,
                "field_degrees": self.field_degrees,
                "is_real": self.is_real,
                "frobenius_schur": self.indicator,
                "certificate": self.certificate}


def verify_certificate(classes: ConjugacyClasses, cert):
    """Independent verifier for the table certificate (the verifier specification V0(d))."""
    e = classes.exponent
    r = classes.r
    S = [i - 1 for i in cert["S"]]
    omegas = [[Cyc.from_json(w) for w in om] for om in cert["omegas"]]
    if len(omegas) != r:
        return False, "wrong number of central characters"
    for i in S:
        M = classes.class_matrix(i)
        for w in omegas:
            for j in range(r):
                lhs = sum((w[k] * M[j][k] for k in range(r) if M[j][k]), Cyc.zero(e))
                if not (lhs == w[j] * w[i]):
                    return False, f"eigenvector identity fails at M_{i+1}"
    tuples = {tuple(w[i] for i in S) for w in omegas}
    if len(tuples) != r:
        return False, "tuples not distinct"
    if any(not (w[0] == 1) for w in omegas):
        return False, "omega(K_1) != 1"
    return True, "ok"