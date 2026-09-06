"""Resolvent polynomials and the alignment of local roots.

Builds the resolvent of an invariant over a transversal from roots computed
modulo a power of a prime, recovering exact integer coefficients by centring;
evaluates invariants at a permuted root vector; and decides which
permutations align a set of local roots with the global numbering.
"""
from __future__ import annotations
import random

from .padic import lifted_roots, frobenius_perm, GaloisRing, irreducible_poly, roots_in_GF, hensel_lift
from .perm import mul, inverse, identity, PermGroup, symmetric
from .ramified import discriminant

def coset_transversal(big: PermGroup, small: PermGroup):
    """Representatives tau of the left cosets tau*small in big (small <= big)."""
    reps = []
    seen = set()
    for x in big.elements():
        if x in seen:
            continue
        reps.append(x)
        for h in small.elements():
            seen.add(mul(x, h))
    return reps

def sorted_roots(ring, roots):
    return sorted(roots, key=lambda x: tuple(ring.residue(x)))

_ROOT_CACHE = {}

NUMBERING = {}   # (tuple(f), p) -> list of residues mod p fixing the order (from the descent certificate)

def roots_at(f, ell, k, seed=0, r=None):
    """Roots of f in Z_{ell^r}/ell^k, sorted by residue (deterministic for a given seed);
    memoized, and served by truncation from a cached higher precision when available."""
    given = NUMBERING.get((tuple(f), ell))
    key = (tuple(f), ell, seed, r, tuple(given) if given else None)
    for (kk, ring, roots, facs) in _ROOT_CACHE.get(key, []):
        if kk == k:
            return ring, roots, facs
        if kk > k:
            ring2 = GaloisRing(ell, k, ring.g)
            return ring2, [ring2.truncate(x, k) for x in roots], facs
    rng = random.Random(seed)
    if r is not None:
        # canonical numbering lives in the minimal ring; embed it into the larger one so that
        # the numbering is independent of the ambient ring
        ring0, roots0, facs = roots_at(f, ell, k, seed=seed, r=None)
        if r % ring0.r:
            raise ValueError("requested degree not a multiple of the root field degree")
        if r == ring0.r:
            return ring0, roots0, facs
        ring = GaloisRing(ell, k, irreducible_poly(ell, r, rng))
        gsmall = [c % ring.mod for c in ring0.g]
        w = roots_in_GF(gsmall, ring.res, rng)[0]
        w = hensel_lift(ring, gsmall, w)
        def emb(a):
            return ring.eval_poly(list(a), w)
        roots = [emb(a) for a in roots0]
        assert all(ring.valuation(ring.eval_poly(f, x)) >= k for x in roots)
        _ROOT_CACHE.setdefault(key, []).append((k, ring, roots, facs))
        return ring, roots, facs
    ring, roots, facs = lifted_roots(f, ell, k, r=r, rng=rng)
    roots = sorted_roots(ring, roots)
    if given is not None:
        # order the roots as in the descent certificate (residues mod ell must be distinct and in Z/ell)
        res = [ring.residue(x) for x in roots]
        order = []
        for g in given:
            j = next(i for i, a in enumerate(res) if a == ring.res.from_int(g % ell))
            order.append(j)
        assert sorted(order) == list(range(len(roots))), "descent numbering does not match the roots"
        roots = [roots[j] for j in order]
    _ROOT_CACHE.setdefault(key, []).append((k, ring, roots, facs))
    return ring, roots, facs

def evaluate_permuted(F, roots, tau, ring):
    """F(tau alpha) = F(alpha_{tau(1)}, ..., alpha_{tau(n)})."""
    return F.evaluate([roots[tau[i]] for i in range(len(roots))], ring)

def resolvent(F, transversal, roots, ring):
    """prod_{tau}(x - F(tau alpha)) with integer coefficients recovered by centering (T1).
    Returns (coefficients ascending, list of the values)."""
    vals = [evaluate_permuted(F, roots, tau, ring) for tau in transversal]
    poly = [ring.one()]
    for v in vals:
        new = [ring.zero() for _ in range(len(poly) + 1)]
        for i, c in enumerate(poly):
            new[i + 1] = ring.add(new[i + 1], c)
            new[i] = ring.sub(new[i], ring.mul(c, v))
        poly = new
    coeffs = []
    for c in poly:
        z = ring.centered_int(c)
        if z is None:
            raise ValueError("resolvent coefficient not in Z/ell^k: not a G-set of values")
        coeffs.append(z)
    return coeffs, vals

def invariant_value(F, roots, ring):
    z = ring.centered_int(F.evaluate(roots, ring))
    if z is None:
        raise ValueError("invariant value not rational at this numbering")
    return z

def alignment_single_level(F, c, transversal, roots, ring):
    """Survivors tau' in the transversal of S_n/G with F(beta^{tau'}) = c mod ell^k (T2)."""
    target = ring.from_int(c)
    out = []
    for idx, tau in enumerate(transversal):
        v = evaluate_permuted(F, roots, tau, ring)
        if all((x - y) % ring.mod == 0 for x, y in zip(v, target)):
            out.append(idx)
    return out


def hadamard_disc_bound(R):
    """|disc R| <= m^m ||R||_2^{2m-2} for monic R of degree m (Mahler/Hadamard)."""
    m = len(R) - 1
    norm2sq = sum(c * c for c in R)
    return m ** m * norm2sq ** (m - 1)

def exact_discriminant(R, vals_at, f, p, seed=0):
    """disc R = prod_{i<j}(r_i - r_j)^2 computed from the resolvent values at the numbering
    prime, at a precision exceeding twice the Hadamard bound (exact by T1).
    vals_at(ring, roots) must return the resolvent values in the given ring."""
    bound = 2 * hadamard_disc_bound(R) + 1
    K, y = 1, p
    while y <= bound:
        y *= p
        K += 1
    ring, roots, _ = roots_at(f, p, K, seed=seed)
    vals = vals_at(ring, roots)
    d = ring.one()
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            t = ring.sub(vals[i], vals[j])
            d = ring.mul(d, ring.mul(t, t))
    z = ring.centered_int(d)
    if z is None:
        raise ValueError("discriminant not rational: values are not the roots of an integer polynomial")
    return z


def squarefree_certificate(R, max_primes=60):
    """A prime ell with gcd(R, R') = 1 modulo ell (hence disc R != 0), or None."""
    from .fpoly import poly_gcd, poly_trim
    from .chartable import is_prime
    q = 2
    tried = 0
    while tried < max_primes:
        while not is_prime(q):
            q += 1
        Rb = [c % q for c in R]
        dRb = poly_trim([(i * c) % q for i, c in enumerate(Rb)][1:])
        if dRb and len(poly_gcd(Rb, dRb, q)) == 1:
            return q
        tried += 1
        q += 1
    return None

def divides_disc(R, ell):
    """ell | disc R  <=>  R mod ell is not squarefree (R monic)."""
    from .fpoly import poly_gcd, poly_trim
    Rb = [c % ell for c in R]
    dRb = poly_trim([(i * c) % ell for i, c in enumerate(Rb)][1:])
    return (not dRb) or len(poly_gcd(Rb, dRb, ell)) > 1


def alignment_chain(chain, roots, ring):
    """Survivors of the multi-level alignment test: chain = [(T_i, F_i, c_i)] with T_i a
    transversal of G_{i-1}/G_i (G_0 = S_n).  Survivors after level i are the products
    tau_1 ... tau_i with F_j(beta^{tau_1...tau_j}) = c_j for all j <= i.  Returns the
    survivors (as permutations) after the last level (exactly one at the policy precision,
    exactly one survivor remains at the policy precision."""
    n = len(roots)
    from .perm import identity, mul
    surv = [identity(n)]
    for T, F, c in chain:
        target = ring.from_int(c)
        new = []
        for s in surv:
            for t in T:
                tau = mul(s, t)
                v = evaluate_permuted(F, roots, tau, ring)
                if all((x - y) % ring.mod == 0 for x, y in zip(v, target)):
                    new.append(tau)
        surv = new
        if not surv:
            break
    return surv