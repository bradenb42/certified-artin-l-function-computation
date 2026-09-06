"""Linear algebra and univariate polynomials over a prime field F_p.
"""
from __future__ import annotations
import random

# ------------------------------------------------------------- matrices mod p

def mat_mul(A, B, p):
    n, m, k = len(A), len(B), len(B[0]) if B else 0
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(row, col)) % p for col in Bt] for row in A]

def mat_vec(A, v, p):
    return [sum(a * b for a, b in zip(row, v)) % p for row in A]

def rref(A, p):
    """Reduced row echelon form; returns (R, pivots)."""
    A = [row[:] for row in A]
    rows, cols = len(A), len(A[0]) if A else 0
    pivots = []
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] % p), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], -1, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] % p:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return A[:r], pivots

def nullspace(A, p, ncols=None):
    """Basis of {x: A x = 0} over F_p (list of vectors)."""
    if ncols is None:
        ncols = len(A[0])
    if not A:
        return [[1 if i == j else 0 for i in range(ncols)] for j in range(ncols)]
    R, piv = rref(A, p)
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fcol in free:
        v = [0] * ncols
        v[fcol] = 1
        for i, c in enumerate(piv):
            v[c] = (-R[i][fcol]) % p
        basis.append(v)
    return basis

def rank(A, p):
    return len(rref(A, p)[1])

def charpoly(A, p):
    """Characteristic polynomial (ascending coefficients) via Hessenberg reduction."""
    n = len(A)
    H = [row[:] for row in A]
    # Hessenberg form
    for m in range(1, n - 1):
        piv = next((i for i in range(m, n) if H[i][m - 1] % p), None)
        if piv is None:
            continue
        if piv != m:
            H[piv], H[m] = H[m], H[piv]
            for row in H:
                row[piv], row[m] = row[m], row[piv]
        inv = pow(H[m][m - 1], -1, p)
        for i in range(m + 1, n):
            if H[i][m - 1] % p:
                f = (H[i][m - 1] * inv) % p
                H[i] = [(x - f * y) % p for x, y in zip(H[i], H[m])]
                for row in H:
                    row[m] = (row[m] + f * row[i]) % p
    # charpoly of Hessenberg matrix by recurrence
    polys = [[1]]  # p_0 = 1
    for m in range(1, n + 1):
        # p_m(x) = (x - H[m-1][m-1]) p_{m-1}(x) - sum_{i=1}^{m-1} H[m-i-1][m-1] * (prod_{j=m-i}^{m-1} H[j][j-1]) p_{m-i-1}(x)
        pm = poly_sub(poly_shift(polys[m - 1]), poly_scale(polys[m - 1], H[m - 1][m - 1], p), p)
        prod = 1
        for i in range(1, m):
            prod = (prod * H[m - i][m - i - 1]) % p
            coef = (H[m - i - 1][m - 1] * prod) % p
            if coef:
                pm = poly_sub(pm, poly_scale(polys[m - i - 1], coef, p), p)
        polys.append(pm)
    return polys[n]

# ---------------------------------------------------------------- polynomials

def poly_trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a

def poly_shift(a):
    return [0] + list(a)

def poly_scale(a, c, p):
    return [(x * c) % p for x in a]

def poly_add(a, b, p):
    n = max(len(a), len(b))
    return poly_trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])

def poly_sub(a, b, p):
    n = max(len(a), len(b))
    return poly_trim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(n)])

def poly_mul(a, b, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return poly_trim(r)

def poly_divmod(a, b, p):
    a = list(a)
    b = poly_trim(list(b))
    if not b:
        raise ZeroDivisionError
    inv = pow(b[-1], -1, p)
    q = [0] * max(0, len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = (a[i + len(b) - 1] * inv) % p
        q[i] = c
        if c:
            for j, y in enumerate(b):
                a[i + j] = (a[i + j] - c * y) % p
    return poly_trim(q), poly_trim(a[:len(b) - 1] if len(b) > 1 else [])

def poly_mod(a, b, p):
    return poly_divmod(a, b, p)[1]

def poly_monic(a, p):
    a = poly_trim(list(a))
    if not a:
        return a
    inv = pow(a[-1], -1, p)
    return poly_scale(a, inv, p)

def poly_gcd(a, b, p):
    a, b = poly_trim(list(a)), poly_trim(list(b))
    while b:
        a, b = b, poly_mod(a, b, p)
    return poly_monic(a, p)

def poly_powmod(base, k, mod, p):
    r = [1]
    base = poly_mod(base, mod, p)
    while k:
        if k & 1:
            r = poly_mod(poly_mul(r, base, p), mod, p)
        base = poly_mod(poly_mul(base, base, p), mod, p)
        k >>= 1
    return r

def poly_eval(a, x, p):
    r = 0
    for c in reversed(a):
        r = (r * x + c) % p
    return r

def roots_mod_p(f, p, rng=random):
    """All roots in F_p of f (with multiplicity ignored)."""
    f = poly_monic(f, p)
    if len(f) <= 1:
        return []
    if p <= 5000:
        return [a for a in range(p) if poly_eval(f, a, p) == 0]
    # g = gcd(f, x^p - x)
    xp = poly_powmod([0, 1], p, f, p)
    g = poly_gcd(f, poly_sub(xp, [0, 1], p), p)
    out = []
    _split_linear(g, p, out, rng)
    return sorted(out)

def _split_linear(g, p, out, rng):
    g = poly_monic(g, p)
    d = len(g) - 1
    if d <= 0:
        return
    if d == 1:
        out.append((-g[0]) % p)
        return
    if d == 2:
        # solve directly via quadratic formula when possible
        pass
    while True:
        a = rng.randrange(p)
        h = poly_powmod([a, 1], (p - 1) // 2, g, p)
        h = poly_sub(h, [1], p)
        w = poly_gcd(g, h, p)
        if 0 < len(w) - 1 < d:
            _split_linear(w, p, out, rng)
            _split_linear(poly_divmod(g, w, p)[0], p, out, rng)
            return
