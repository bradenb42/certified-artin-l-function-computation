"""Multivariate invariants used to separate subgroups.

An invariant is a sparse integer polynomial in the roots, together with an
optional Tschirnhaus substitution.  The module provides the permutation
action, the height data that drives the precision policy, and the standard
constructions: products over a subset of the roots and orbit sums, whose
stabiliser is a prescribed subgroup.
"""
from __future__ import annotations
from math import comb

from .perm import identity, mul, inverse

class Invariant:
    """F = sum c_m x^m stored as dict exponent-tuple -> int.  Optional Tschirnhaus
    transform T (ascending integer coefficients): the invariant is F(T(x_1),..,T(x_n))
    and its the precision policy height data are those of F o T."""
    def __init__(self, n, terms, T=None, label=None):
        self.n = n
        self.terms = {tuple(m): int(c) for m, c in terms.items() if c}
        self.T = list(T) if T else None
        self.label = label

    # -- the precision policy data
    def degree(self):
        d = max((sum(m) for m in self.terms), default=0)
        if self.T:
            d *= len(self.T) - 1
        return d

    def norm1(self):
        s = sum(abs(c) for c in self.terms.values())
        if self.T:
            s *= sum(abs(c) for c in self.T) ** max((sum(m) for m in self.terms), default=0)
        return s

    def height_bound(self, R):
        """B_F = ||F||_1 R^{deg F}."""
        return self.norm1() * R ** self.degree()

    # -- action: (rho F)(x) = F(x_{rho(1)},...,x_{rho(n)}), i.e. x_i -> x_{rho(i)}
    def act(self, rho):
        out = {}
        for m, c in self.terms.items():
            m2 = [0] * self.n
            for i, e in enumerate(m):
                if e:
                    m2[rho[i]] += e
            m2 = tuple(m2)
            out[m2] = out.get(m2, 0) + c
        return Invariant(self.n, out, self.T, self.label)

    def __eq__(self, other):
        return self.terms == other.terms and self.T == other.T

    def __hash__(self):
        return hash(frozenset(self.terms.items()))

    def evaluate(self, vals, ring):
        """Evaluate at vals (list of ring elements); ring provides add, mul, one, zero, from_int."""
        if self.T:
            vals = [_horner(self.T, v, ring) for v in vals]
        total = ring.zero()
        powcache = {}
        for m, c in self.terms.items():
            term = ring.from_int(c)
            for i, e in enumerate(m):
                if e:
                    key = (i, e)
                    if key not in powcache:
                        powcache[key] = ring.power(vals[i], e)
                    term = ring.mul(term, powcache[key])
            total = ring.add(total, term)
        return total

    def to_json(self):
        return {"n": self.n, "terms": [[list(m), c] for m, c in self.terms.items()],
                "tschirnhaus": self.T, "label": self.label,
                "degree": self.degree(), "norm1": self.norm1()}

    @staticmethod
    def from_json(d):
        return Invariant(d["n"], {tuple(m): c for m, c in d["terms"]}, d.get("tschirnhaus"), d.get("label"))

def _horner(T, v, ring):
    r = ring.zero()
    for c in reversed(T):
        r = ring.add(ring.mul(r, v), ring.from_int(c))
    return r

def set_product(n, S, label=None):
    """prod_{i in S} x_i."""
    m = [0] * n
    for i in S:
        m[i] = 1
    return Invariant(n, {tuple(m): 1}, label=label)

def orbit_sum(n, H_elements, base, label=None):
    """F_H = sum_{h in H} h. m_t with m_t = prod_j x_{t_j}^{j+1}."""
    terms = {}
    for h in H_elements:
        m = [0] * n
        for j, t in enumerate(base):
            m[h[t]] += j + 1
        m = tuple(m)
        terms[m] = terms.get(m, 0) + 1
    return Invariant(n, terms, label=label)

def stabilizer_in(F, group_elements):
    return [g for g in group_elements if F.act(g) == F]
