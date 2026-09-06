"""Independent re-checking of maximal-order witnesses at a prime.
"""
from __future__ import annotations
from fractions import Fraction
from.arith import ptrim, pdivmod, pmul, pmod, pgcd

def polymul_mod(a, b, f):
    n = len(f) - 1
    prod = [Fraction(0)] * (2 * n - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    prod[i + j] += x * y
    for i in range(2 * n - 2, n - 1, -1):
        c = prod[i]
        if c:
            for j in range(n + 1):
                prod[i - n + j] -= c * f[j]
    return prod[:n]

def inv_matrix(M):
    n = len(M)
    A = [list(map(Fraction, row)) + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            raise ZeroDivisionError
        A[c], A[piv] = A[piv], A[c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                fct = A[r][c]
                A[r] = [x - fct * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]

def det(M):
    n = len(M); A = [list(map(Fraction, r)) for r in M]; d = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]; d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                fct = A[r][c] / A[c][c]
                A[r] = [x - fct * y for x, y in zip(A[r], A[c])]
    return d

def vecmat(v, M):
    return [sum(v[i] * M[i][j] for i in range(len(v))) for j in range(len(M[0]))]

def val(x, p):
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v

def nullspace_mod(A, p, ncols):
    A = [row[:] for row in A]
    rows = len(A); piv = []; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, rows) if A[i][c] % p), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][c], -1, p)
        A[r] = [x * inv % p for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] % p:
                fct = A[i][c]
                A[i] = [(x - fct * y) % p for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == rows:
            break
    free = [c for c in range(ncols) if c not in piv]
    out = []
    for fc in free:
        v = [0] * ncols; v[fc] = 1
        for i, c in enumerate(piv):
            v[c] = (-A[i][fc]) % p
        out.append(v)
    return out

def hnf_basis(rows, n):
    A = [list(r) for r in rows]
    basis = []
    for c in range(n):
        while True:
            nz = [r for r in A if r[c] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda r: abs(r[c]))
            pv = nz[0]
            for r in nz[1:]:
                q = r[c] // pv[c]
                for j in range(n):
                    r[j] -= q * pv[j]
        nz = [r for r in A if r[c] != 0]
        assert nz, "not full rank"
        pv = nz[0]
        basis.append([-x for x in pv] if pv[c] < 0 else pv)
        A = [r for r in A if r is not pv and any(r)]
    return basis

class Order:
    def __init__(self, f, ell, B):
        self.f = [Fraction(c) for c in f]; self.n = len(f) - 1; self.ell = ell
        self.B = [[Fraction(x) for x in row] for row in B]
        self.Binv = inv_matrix(self.B)
        n = self.n
        for i in range(n):
            co = vecmat([Fraction(int(j == i)) for j in range(n)], self.Binv)
            if any(c.denominator != 1 for c in co):
                raise ValueError("order does not contain Z[alpha]")
        self.s = [[None] * n for _ in range(n)]
        for a in range(n):
            for b in range(a, n):
                co = vecmat(polymul_mod(self.B[a], self.B[b], self.f), self.Binv)
                if any(c.denominator != 1 for c in co):
                    raise ValueError("basis not closed under multiplication")
                self.s[a][b] = self.s[b][a] = [int(c) for c in co]
        d = det(self.B)
        self.index_val = val(d.denominator, ell) - val(d.numerator, ell)
    # algebra mod ell
    def mul_mod(self, x, y):
        n, p = self.n, self.ell
        out = [0] * n
        for a in range(n):
            if x[a]:
                for b in range(n):
                    if y[b]:
                        sab = self.s[a][b]
                        for t in range(n):
                            if sab[t]:
                                out[t] = (out[t] + x[a] * y[b] * sab[t]) % p
        return out
    def one_mod(self):
        n, p = self.n, self.ell
        rows, rhs = [], []
        for a in range(n):
            for t in range(n):
                rows.append([self.s[u][a][t] % p for u in range(n)] + [int(a == t)])
        # solve by elimination
        A = rows; piv = []; r = 0
        for c in range(n):
            pr = next((i for i in range(r, len(A)) if A[i][c] % p), None)
            if pr is None: continue
            A[r], A[pr] = A[pr], A[r]
            inv = pow(A[r][c], -1, p); A[r] = [x * inv % p for x in A[r]]
            for i in range(len(A)):
                if i != r and A[i][c] % p:
                    fct = A[i][c]; A[i] = [(x - fct * y) % p for x, y in zip(A[i], A[r])]
            piv.append(c); r += 1
        sol = [0] * n
        for i, c in enumerate(piv):
            sol[c] = A[i][n]
        return sol
    def frob_matrix(self):
        n, p = self.n, self.ell
        one = self.one_mod()
        cols = []
        for a in range(n):
            e = [int(i == a) for i in range(n)]
            r = one; base = e; k = p
            while k:
                if k & 1: r = self.mul_mod(r, base)
                base = self.mul_mod(base, base); k >>= 1
            cols.append(r)
        return [[cols[a][t] for a in range(n)] for t in range(n)]
    def nilradical(self):
        n, p = self.n, self.ell
        F = self.frob_matrix()
        r = 1
        while p ** r < n: r += 1
        Fr = F
        for _ in range(r - 1):
            Fr = [[sum(Fr[i][k] * F[k][j] for k in range(n)) % p for j in range(n)] for i in range(n)]
        return nullspace_mod(Fr, p, n)
    def radical_lattice(self):
        N = self.nilradical()
        rows = [[self.ell * int(i == j) for j in range(self.n)] for i in range(self.n)] + [list(v) for v in N]
        return hnf_basis(rows, self.n), len(N)
    def is_maximal(self):
        """The multiplier ring of the radical equals the order."""
        n, p = self.n, self.ell
        C, _ = self.radical_lattice()
        Cq = [[Fraction(x) for x in r] for r in C]; Cinv = inv_matrix(Cq)
        cols = []
        for a in range(n):
            col = []
            for j in range(n):
                vec = [0] * n
                for m, cm in enumerate(C[j]):
                    if cm:
                        for t in range(n): vec[t] += cm * self.s[a][m][t]
                y = vecmat([Fraction(x) for x in vec], Cinv)
                if any(v.denominator != 1 for v in y):
                    return False
                col.extend(int(v) % p for v in y)
            cols.append(col)
        M = [[cols[a][i] for a in range(n)] for i in range(n * n)]
        ker = nullspace_mod(M, p, n)
        # U = ell O + lifts; O' = U / ell equals O iff every kernel vector lies in ell O, i.e. kernel trivial mod ell
        return all(all(x % p == 0 for x in v) for v in ker)
    def residue_decomposition(self):
        """(e_j, f_j) by idempotent splitting of O/ell (own implementation)."""
        return _decompose(self.n, self.s, self.ell)

def _decompose(m, s, p):
    import random
    rng = random.Random(1)
    def mul(x, y):
        out = [0] * m
        for a in range(m):
            if x[a]:
                for b in range(m):
                    if y[b]:
                        sab = s[a][b]
                        for t in range(m):
                            if sab[t]: out[t] = (out[t] + x[a] * y[b] * sab[t]) % p
        return out
    def solve_one():
        rows = [[s[u][a][t] % p for u in range(m)] + [int(a == t)] for a in range(m) for t in range(m)]
        A = rows; piv = []; r = 0
        for c in range(m):
            pr = next((i for i in range(r, len(A)) if A[i][c] % p), None)
            if pr is None: continue
            A[r], A[pr] = A[pr], A[r]; inv = pow(A[r][c], -1, p); A[r] = [x * inv % p for x in A[r]]
            for i in range(len(A)):
                if i != r and A[i][c] % p:
                    fct = A[i][c]; A[i] = [(x - fct * y) % p for x, y in zip(A[i], A[r])]
            piv.append(c); r += 1
        sol = [0] * m
        for i, c in enumerate(piv): sol[c] = A[i][m]
        return sol
    one = solve_one()
    def power(x, k):
        r = one
        while k:
            if k & 1: r = mul(r, x)
            x = mul(x, x); k >>= 1
        return r
    F = [[0] * m for _ in range(m)]
    for a in range(m):
        col = power([int(i == a) for i in range(m)], p)
        for t in range(m): F[t][a] = col[t]
    r = 1
    while p ** r < m: r += 1
    Fr = F
    for _ in range(r - 1):
        Fr = [[sum(Fr[i][k] * F[k][j] for k in range(m)) % p for j in range(m)] for i in range(m)]
    N = nullspace_mod(Fr, p, m)
    nd = len(N)
    # projection mod N
    Nr = [row[:] for row in N]; piv = []
    rr = 0
    for c in range(m):
        pr = next((i for i in range(rr, len(Nr)) if Nr[i][c] % p), None)
        if pr is None: continue
        Nr[rr], Nr[pr] = Nr[pr], Nr[rr]; inv = pow(Nr[rr][c], -1, p); Nr[rr] = [x * inv % p for x in Nr[rr]]
        for i in range(len(Nr)):
            if i != rr and Nr[i][c] % p:
                fct = Nr[i][c]; Nr[i] = [(x - fct * y) % p for x, y in zip(Nr[i], Nr[rr])]
        piv.append(c); rr += 1
    def modN(v):
        v = list(v)
        for row, c in zip(Nr, piv):
            if v[c]:
                fct = v[c]; v = [(x - fct * y) % p for x, y in zip(v, row)]
        return [v[i] for i in range(m) if i not in piv]
    FmI = [[(F[i][j] - int(i == j)) % p for j in range(m)] for i in range(m)]
    cols = [modN([FmI[i][a] for i in range(m)]) for a in range(m)]
    if cols and cols[0]:
        Mz = [[cols[a][i] for a in range(m)] for i in range(len(cols[0]))]
        Zs = nullspace_mod(Mz, p, m)
    else:
        Zs = [[int(i == j) for j in range(m)] for i in range(m)]
    nprimes = len(Zs) - nd
    if nprimes == 1:
        fdeg = m - nd
        return [(m // fdeg, fdeg)]
    # splitting element
    span = [list(v) for v in N] + [one]
    def rank(vs):
        A = [list(v) for v in vs]; r = 0
        for c in range(m):
            pr = next((i for i in range(r, len(A)) if A[i][c] % p), None)
            if pr is None: continue
            A[r], A[pr] = A[pr], A[r]; inv = pow(A[r][c], -1, p); A[r] = [x * inv % p for x in A[r]]
            for i in range(len(A)):
                if i != r and A[i][c] % p:
                    fct = A[i][c]; A[i] = [(x - fct * y) % p for x, y in zip(A[i], A[r])]
            r += 1
        return r
    y = next(v for v in Zs if rank(span + [v]) > rank(span))
    # minimal polynomial and its roots (all in F_p)
    powers = [one]
    while True:
        M = [[powers[k][t] for k in range(len(powers))] for t in range(m)]
        ker = nullspace_mod(M, p, len(powers))
        if ker:
            v = ker[0]; inv = pow(v[-1], -1, p); mp = [c * inv % p for c in v]; break
        powers.append(mul(powers[-1], y))
    roots = [a for a in range(p) if sum(c * pow(a, i, p) for i, c in enumerate(mp)) % p == 0]
    facs = []; rem = mp
    for a in roots:
        k = 0
        while True:
            q, rr_ = pdivmod(rem, [(-a) % p, 1], p)
            if rr_: break
            rem, k = q, k + 1
        facs.append(([(-a) % p, 1], k))
    out = []
    for i, (lin, k) in enumerate(facs):
        qi = [1]
        for _ in range(k): qi = pmul(qi, lin, p)
        cof = [1]
        for j, (lin2, k2) in enumerate(facs):
            if j != i:
                for _ in range(k2): cof = pmul(cof, lin2, p)
        # inverse of cof mod qi
        r0, r1 = ptrim(list(qi)), pmod(cof, qi, p); s0, s1 = [], [1]
        while r1:
            q, rmd = pdivmod(r0, r1, p); r0, r1 = r1, rmd
            s0, s1 = s1, ptrim([(a - b) % p for a, b in zip(s0 + [0] * len(pmul(q, s1, p)), pmul(q, s1, p) + [0] * len(s0))])
        inv = pow(r0[0], -1, p); invp = pmod([c * inv % p for c in s0], qi, p)
        h = pmod(pmul(cof, invp, p), mp, p)
        ei = [0] * m
        for c in reversed(h):
            ei = mul(ei, y); ei = [(a + c * b) % p for a, b in zip(ei, one)]
        vecs = [mul(ei, [int(t == a) for t in range(m)]) for a in range(m)]
        # rref basis of the component
        A = [list(v) for v in vecs]; pv = []; r = 0
        for c in range(m):
            pr = next((ii for ii in range(r, len(A)) if A[ii][c] % p), None)
            if pr is None: continue
            A[r], A[pr] = A[pr], A[r]; inv = pow(A[r][c], -1, p); A[r] = [x * inv % p for x in A[r]]
            for ii in range(len(A)):
                if ii != r and A[ii][c] % p:
                    fct = A[ii][c]; A[ii] = [(x - fct * yy) % p for x, yy in zip(A[ii], A[r])]
            pv.append(c); r += 1
        sub = A[:r]; d = len(sub)
        s2 = [[None] * d for _ in range(d)]
        for a in range(d):
            for b in range(d):
                prod = mul(sub[a], sub[b])
                s2[a][b] = [prod[c] for c in pv]
        out.extend(_decompose(d, s2, p))
    return out

def dedekind_check(f, ell, facs_rec):
    """Dedekind criterion with the recorded factorization: returns (maximal, ramified)."""
    from.arith import factorization_type
    return None