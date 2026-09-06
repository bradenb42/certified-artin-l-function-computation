"""Euler factors at ramified primes, and the identities that check them.

Provides the determinant character, the Euler factor at a ramified prime by
two independent routes (isotypic multiplicities over the Frobenius coset, and
power traces with Newton's identities), the factor from a matrix model, and
the identity relating the product of the Euler factors over the characters
appearing in a subgroup's permutation character to the local factor of the
Dedekind zeta function of the corresponding subfield.
"""
from __future__ import annotations
from fractions import Fraction
from math import gcd

from .cyclo import Cyc, phi
from .perm import PermGroup, mul, inverse, identity
from .schur import eigen_multiplicities, multiplicity, Pair

# ------------------------------------------------------------------ polynomials with Cyc coefficients

def cpoly_mul(a, b):
    e = max(a[0].e, b[0].e)
    out = [Cyc.zero(e) for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = out[i + j] + x * y
    return out

def cpoly_pow(a, n):
    r = [Cyc.one(a[0].e)]
    for _ in range(n):
        r = cpoly_mul(r, a)
    return r

def cpoly_eq(a, b):
    n = max(len(a), len(b))
    for i in range(n):
        x = a[i] if i < len(a) else Cyc.zero(a[0].e)
        y = b[i] if i < len(b) else Cyc.zero(b[0].e)
        if not (x == y):
            return False
    return True

def cpoly_json(a):
    return [c.to_json() for c in a]

# ------------------------------------------------------------------ determinant character

def det_character(table, nu):
    """det rho_chi as a class function: on g_k, zeta_o^{sum_j j m_j}."""
    cl, e = table.cl, table.e
    vals = []
    for k in range(cl.r):
        o = cl.orders[k]
        m = eigen_multiplicities(table, nu, k)
        s = sum(j * mj for j, mj in enumerate(m)) % o
        vals.append(Cyc.root_of_unity(e, o, s))
    return vals

def conductor_exponent_of_class_function(vals, deg, delta, N_by_class, e_inertia):
    """f_ell for a class function given by its values."""
    val = vals[0] * 0 + deg * delta
    for k, Nk in N_by_class.items():
        val = val - vals[k] * Nk
    q = val.rational() / e_inertia
    assert q.denominator == 1
    return int(q)

# ------------------------------------------------------------------ the ramified-Euler-factor construction: Euler factor at a ramified prime

def euler_factor_C(table, nu, D_elems_global, I_elems_global, Fr_global, cl, f):
    """Route (C): P(T) = prod_j (1 - zeta_f^j T)^{<chi|_D, psi_j>}, psi_j(Fr^k i) = zeta_f^{jk}."""
    e = table.e
    E = e * f // gcd(e, f)
    # coset exponent k of each sigma in D: sigma I = Fr^k I
    Iset = set(I_elems_global)
    k_of = {}
    x = identity(len(Fr_global))
    for k in range(f):
        for i in I_elems_global:
            k_of[mul(x, i)] = k
        x = mul(x, Fr_global)
    assert set(k_of) == set(D_elems_global)
    P = [Cyc.one(E)]
    mults = []
    for j in range(f):
        s = Cyc.zero(E)
        for sg in D_elems_global:
            s = s + table.values[nu][cl.class_of(sg)].embed(E) * Cyc.root_of_unity(E, f, (-j * k_of[sg]) % f)
        q = s.rational() / len(D_elems_global)
        assert q.denominator == 1 and q >= 0
        m = int(q)
        mults.append(m)
        for _ in range(m):
            P = cpoly_mul(P, [Cyc.one(E), -Cyc.root_of_unity(E, f, j)])
    return P, mults

def euler_factor_B(table, nu, I_elems_global, Fr_global, cl, f, dimVI):
    """Route (B): power sums s_k = tr(rho(Fr)^k | V^I) = (1/|I|) sum_i chi(Fr^k i), Newton identities."""
    e = table.e
    E = e * f // gcd(e, f)
    pw = [Cyc.zero(E)]
    x = identity(len(Fr_global))
    for k in range(1, dimVI + 1):
        x = mul(x, Fr_global)
        s = Cyc.zero(E)
        for i in I_elems_global:
            s = s + table.values[nu][cl.class_of(mul(x, i))].embed(E)
        pw.append(s / len(I_elems_global))
    # Newton: e_0 = 1, k e_k = sum_{i=1}^k (-1)^{i-1} e_{k-i} s_i
    el = [Cyc.one(E)]
    for k in range(1, dimVI + 1):
        acc = Cyc.zero(E)
        for i in range(1, k + 1):
            term = el[k - i] * pw[i]
            acc = acc + (term if i % 2 == 1 else -term)
        el.append(acc / k)
    # P(T) = sum (-1)^k e_k T^k
    return [el[k] if k % 2 == 0 else -el[k] for k in range(dimVI + 1)]

def euler_factor_A(rho, dim_model, field, I_elems_global, Fr_global):
    """Route (A): with a matrix model rho (multiplicity a), e_I = sum_{i in I} rho(i), the
    characteristic polynomial of rho(Fr) on the image of e_I (equals P^a)."""
    F = field
    n = dim_model
    Z = [[Cyc.zero(F) for _ in range(n)] for _ in range(n)]
    for i in I_elems_global:
        M = rho(i)
        Z = [[Z[r][c] + M[r][c] for c in range(n)] for r in range(n)]
    # column space basis of Z (over Q(zeta_F))
    basis = []
    for col in range(n):
        v = [Z[r][col] for r in range(n)]
        for piv, b in basis:
            if not v[piv].is_zero():
                fct = v[piv]
                v = [x - fct * y for x, y in zip(v, b)]
        piv = next((r for r in range(n) if not v[r].is_zero()), None)
        if piv is None:
            continue
        inv = v[piv].inverse()
        v = [x * inv for x in v]
        basis = [(pv, [x - b[piv] * y for x, y in zip(b, v)]) for pv, b in basis]
        basis.append((piv, v))
    d = len(basis)
    if d == 0:
        return [Cyc.one(F)], 0
    Mf = rho(Fr_global)
    pivots = [pv for pv, _ in basis]
    vecs = [b for _, b in basis]
    # action matrix A: image of basis vector j = Mf b_j, coordinates at pivots
    A = []
    for b in vecs:
        img = [sum((Mf[r][c] * b[c] for c in range(n) if not b[c].is_zero()), Cyc.zero(F)) for r in range(n)]
        co = [img[pv] for pv in pivots]
        # check in span
        w = [Cyc.zero(F) for _ in range(n)]
        for ci, bb in zip(co, vecs):
            w = [x + ci * y for x, y in zip(w, bb)]
        assert all(x == y for x, y in zip(img, w)), "V^I not Fr-stable?"
        A.append(co)   # A[j][i]: coefficient of b_i in image of b_j
    # det(1 - T A^T) via Faddeev-LeVerrier over Q(zeta_F) (dimension small)
    At = [[A[j][i] for j in range(d)] for i in range(d)]
    # characteristic polynomial coefficients c_k of det(xI - At): Faddeev-LeVerrier
    Mk = [[Cyc.zero(F) for _ in range(d)] for _ in range(d)]
    c = [Cyc.one(F)]
    I_d = [[Cyc.one(F) if i == j else Cyc.zero(F) for j in range(d)] for i in range(d)]
    for k in range(1, d + 1):
        # M_k = A M_{k-1} + c_{k-1} I
        Mk = [[sum((At[i][t] * Mk[t][j] for t in range(d)), Cyc.zero(F)) + (c[-1] if i == j else Cyc.zero(F)) for j in range(d)] for i in range(d)]
        tr = sum((sum((At[i][t] * Mk[t][i] for t in range(d)), Cyc.zero(F)) for i in range(d)), Cyc.zero(F))
        c.append(-(tr / k))
    # det(xI - A) = x^d + c_1 x^{d-1} + ... ; det(1 - T A) = T^d det(1/T - A) = sum c_k T^k
    return c, d

def check_euler(table, nu, P, mults, dimVI, f_ell, swan, fdeg):
    d = len(P) - 1
    ok = (d == dimVI == table.degrees[nu] - f_ell + swan)
    return ok

# ------------------------------------------------------------------ the Euler-factor identity system: identities

def subgroup_zeta_factor(H_pair: Pair, table, cl, euler_polys, degrees):
    """prod_chi P_ell(chi;T)^{n_chi(H)}, n_chi(H) = <chi|_H, 1>."""
    E = max(P[0].e for P in euler_polys)
    prod = [Cyc.one(E)]
    ns = []
    for nu, P in enumerate(euler_polys):
        n = multiplicity(table, nu, H_pair)
        ns.append(n)
        for _ in range(n):
            prod = cpoly_mul(prod, [c.embed(E) for c in P])
    return prod, ns

def direct_zeta_factor(residue_decomposition, E):
    """Z_{H,ell}(T) = prod_{p | ell} (1 - T^{f(p/ell)}) from the (e_j, f_j) of R_H at ell."""
    prod = [Cyc.one(E)]
    for e_j, f_j in residue_decomposition:
        fac = [Cyc.zero(E) for _ in range(f_j + 1)]
        fac[0] = Cyc.one(E); fac[f_j] = -Cyc.one(E)
        prod = cpoly_mul(prod, fac)
    return prod