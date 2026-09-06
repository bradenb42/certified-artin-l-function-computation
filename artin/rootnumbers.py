"""Local root numbers at tamely ramified primes.

For each character of the inertia group appearing in chi, the local root
number is a determinant times a Gauss sum divided by the square root of the
residue field size.  The module fixes the normalisation through Lubin-Tate
theory, evaluates the determinant exactly and the Gauss sum through Gaussian
periods, and checks the modulus of every Gauss sum.
"""
from __future__ import annotations
from fractions import Fraction
from math import gcd
import cmath, math

from .cyclo import Cyc, cyclotomic_poly
from .padic import hensel_lift, roots_in_GF
from .perm import mul, inverse, identity
from .schur import Pair, multiplicity
from .filtration import HardFailure
from .euler import det_character

def lcm(a, b):
    return a * b // gcd(a, b)

def mult_order(a, m):
    if m == 1:
        return 1
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k

# ------------------------------------------------------------------ exact sqrt(ell) and Gauss sums in cyclotomic fields

def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def sqrt_prime(ell, M):
    """sqrt(ell) as an element of Q(zeta_M) (M must be divisible by 8 if ell = 2, by 4 ell otherwise)."""
    if ell == 2:
        assert M % 8 == 0
        return Cyc.zeta(M, M // 8) + Cyc.zeta(M, 7 * M // 8)
    assert M % (4 * ell) == 0
    # build in Q(zeta_ell), embed once
    g = Cyc(ell, [0] + [legendre(a, ell) for a in range(1, ell)]).embed(M)
    if ell % 4 == 1:
        return g
    # g = i sqrt(ell): sqrt(ell) = -i g
    return g * (-Cyc.zeta(M, M // 4))

def gauss_sum(theta_values, trace_values, ell, M):
    """g(theta) = sum_{u in k^x} theta(u) zeta_ell^{Tr(u)} for theta given by its values (Cyc in Q(zeta_M))
    indexed like trace_values (integers mod ell)."""
    s = Cyc.zero(M)
    for th, tr in zip(theta_values, trace_values):
        s = s + th * Cyc.zeta(M, (M // ell) * (tr % ell))
    return s

# ------------------------------------------------------------------ the tame root-number construction at a tame prime

class TameRootNumbers:
    def __init__(self, local, rec, table, cl, rho_match, log=print):
        """local: LocalGalois (K with phi, tau); rec: the local stage record; rho_match: matching (global i <-> local rho(i))."""
        self.local, self.table, self.cl = local, table, cl
        K = local.K
        self.ell, self.E, self.e = local.ell, local.E, local.e
        self.rho = rho_match
        # eta with eta^E = -1 in the Galois ring
        R = K.ring
        if self.E % 2:
            self.eta = R.from_int(-1)
        else:
            poly = [1] + [0] * (self.E - 1) + [1]
            w = roots_in_GF([c % self.ell for c in poly], R.res)
            if not w:
                raise HardFailure("mu_{2E} not in the residue field: enlarge F")
            self.eta = hensel_lift(R, poly, w[0])
        # residue field data: generator of F_{ell^F}^x
        self.res = R.res
        self.gen = self._primitive_element()
        self.M = None

    def _primitive_element(self):
        """A generator of F_{ell^F}^x by random sampling (order test on the prime factors of q-1)."""
        import random
        K = self.res
        q = K.q
        m, fac = q - 1, []
        d = 2
        while d * d <= m:
            if m % d == 0:
                fac.append(d)
                while m % d == 0:
                    m //= d
            d += 1
        if m > 1:
            fac.append(m)
        rng = random.Random(0)
        for _ in range(10000):
            cand = [rng.randrange(K.p) for _ in range(K.r)]
            if K.is_zero(cand):
                continue
            if all(not K.eq(K.power(cand, (q - 1) // p), K.one()) for p in fac):
                return cand
        raise RuntimeError("no primitive element")

    def to_global(self, perm):
        return mul(inverse(self.rho), mul(perm, self.rho))

    def epsilon(self, nu, log=print):
        """epsilon_ell(chi_nu) as a complex number (double precision), with the orbit data.
        The determinant on V_theta is exact (Q(zeta_e) arithmetic, then evaluated); the Gauss sums
        are evaluated numerically (error < q * 1e-15, far below the separation of the finite
        candidate sets used for self-dual characters in the root-number stage)."""
        table, cl, local = self.table, self.cl, self.local
        e, E, ell = self.e, self.E, self.ell
        K = local.K
        Iel = local.I.elements()
        tau = local.tau
        powers = {}
        x = identity(len(tau))
        for a in range(e):
            powers[x] = a
            x = mul(x, tau)
        assert len(powers) == e
        eps = 1 + 0j
        h = 0
        orbits_done = set()
        details = []
        for j in range(1, e):
            if j in orbits_done:
                continue
            orb, t = [], j
            while t not in orb:
                orb.append(t)
                t = (t * ell) % e
            orbits_done |= set(orb)
            s = Cyc.zero(e)
            for i in Iel:
                s = s + table.values[nu][cl.class_of(self.to_global(i))] * Cyc.root_of_unity(e, e, (-j * powers[i]) % e)
            q = s.rational() / e
            assert q.denominator == 1
            n_theta = int(q)
            if n_theta == 0:
                continue
            e_theta = e // gcd(e, j)
            f_theta = mult_order(ell, e_theta)
            qF = ell ** f_theta
            b = 0 if ((qF - 1) // e_theta) % 2 == 0 else e_theta // 2
            b_tau = b * (e // e_theta)
            Phi = identity(len(tau))
            for _ in range(f_theta):
                Phi = mul(local.phi, Phi)
            for _ in range(b_tau):
                Phi = mul(tau, Phi)
            def trace_on_Vtheta(g):
                s2 = Cyc.zero(e)
                for i in Iel:
                    s2 = s2 + table.values[nu][cl.class_of(self.to_global(mul(g, i)))] * Cyc.root_of_unity(e, e, (-j * powers[i]) % e)
                return s2 / e
            pw = [None]
            g = identity(len(tau))
            for k in range(1, n_theta + 1):
                g = mul(g, Phi)
                pw.append(trace_on_Vtheta(g))
            el = [Cyc.one(e)]
            for k in range(1, n_theta + 1):
                acc = Cyc.zero(e)
                for i in range(1, k + 1):
                    term = el[k - i] * pw[i]
                    acc = acc + (term if i % 2 == 1 else -term)
                el.append(acc / k)
            det = el[n_theta].to_complex()
            # theta_F on k_F^x and the Gauss sum, numerically
            Kr = self.res
            gF = Kr.power(self.gen, (Kr.q - 1) // (qF - 1))
            zres = Kr.power([c % ell for c in K.zeta], E // e_theta)
            dlog = {}
            z = Kr.one()
            for m in range(e_theta):
                dlog[tuple(z)] = m
                z = Kr.mul(z, zres)
            periods = [0j] * e_theta
            u = Kr.one()
            for _ in range(qF - 1):
                w = Kr.power(u, (qF - 1) // e_theta)
                m = dlog[tuple(w)]
                tr = Kr.zero(); y = u
                for _ in range(f_theta):
                    tr = Kr.add(tr, y); y = Kr.frob(y)
                periods[m] += cmath.exp(2j * cmath.pi * (tr[0] % ell) / ell)
                u = Kr.mul(u, gF)
            # g(conj theta_F) = sum_m conj(zeta_e^{j(-m)}) eta_m
            g_bar = sum(cmath.exp(-2j * cmath.pi * ((j * (-m)) % e) / e) * periods[m] for m in range(e_theta))
            if abs(abs(g_bar) ** 2 - qF) > 1e-6 * qF:
                raise HardFailure(f"Gauss sum modulus check fails at {ell}: |g|^2 = {abs(g_bar)**2} vs {qF}")
            eps *= det * (g_bar / math.sqrt(qF)) ** n_theta
            h += f_theta * n_theta
            details.append({"theta_exponent": j, "orbit": orb, "e_theta": e_theta, "f_theta": f_theta, "n_theta": n_theta, "b": b})
        return eps, h, details

# ------------------------------------------------------------------ the global root-number construction global

def global_root_number(eps_inf_b, eps_locals, M):
    """W = i^{-b} prod_ell num_ell / ell^{h_ell/2}, exactly in Q(zeta_M)."""
    W = Cyc.root_of_unity(M, 4, (-eps_inf_b) % 4)
    for ell, (num, h) in eps_locals.items():
        W = W * num.embed(M)
        if h:
            s = sqrt_prime(ell, M)
            W = W * (s ** h).inverse()
    return W