"""Permutation groups.

Permutations are tuples of images on the points 0..n-1, composed by
(p*q)(i) = p[q[i]].  On disk they are 1-indexed image lists.  PermGroup is a
deterministic Schreier-Sims implementation providing order, membership,
element enumeration and conjugacy classes with power maps and class matrices;
these are the only group operations the rest of the package uses.
"""
from __future__ import annotations
import random
from math import gcd
from typing import Iterable

# ---------------------------------------------------------------- permutations

def identity(n):
    return tuple(range(n))

def mul(p, q):
    """(p*q)(i) = p[q[i]]."""
    return tuple(p[i] for i in q)

def inverse(p):
    r = [0] * len(p)
    for i, j in enumerate(p):
        r[j] = i
    return tuple(r)

def power(p, k):
    n = len(p)
    if k < 0:
        p, k = inverse(p), -k
    r = identity(n)
    while k:
        if k & 1:
            r = mul(r, p)
        p = mul(p, p)
        k >>= 1
    return r

def conj(g, x):
    """g x g^-1."""
    return mul(mul(g, x), inverse(g))

def cycles(p):
    n = len(p)
    seen = [False] * n
    out = []
    for i in range(n):
        if not seen[i]:
            c = []
            j = i
            while not seen[j]:
                seen[j] = True
                c.append(j)
                j = p[j]
            out.append(c)
    return out

def cycle_type(p):
    return tuple(sorted((len(c) for c in cycles(p)), reverse=True))

def order(p):
    o = 1
    for c in cycles(p):
        o = o * len(c) // gcd(o, len(c))
    return o

def from_cycles(n, cyc):
    """cyc: list of cycles on 0-indexed points."""
    p = list(range(n))
    for c in cyc:
        for a, b in zip(c, c[1:] + c[:1]):
            p[a] = b
    return tuple(p)

def to_json(p):
    return [i + 1 for i in p]

def from_json(lst):
    return tuple(i - 1 for i in lst)

# ------------------------------------------------------------ Schreier-Sims

class PermGroup:
    """A permutation group with a base and strong generating set.

    Interface used by the rest of the pipeline:
      n, generators, order(), contains(g), elements() (enumeration, guarded
      by a size limit), classes()  -> ConjugacyClasses, random_element().
    """

    def __init__(self, generators: Iterable, n=None):
        gens = [tuple(g) for g in generators]
        if n is None:
            if not gens:
                raise ValueError("need n for the trivial group")
            n = len(gens[0])
        self.n = n
        self.generators = [g for g in gens if g != identity(n)]
        self._base = []
        self._transversals = []   # list of dict point -> perm mapping base point to point
        self._sgs_levels = []     # list of list of strong generators fixing earlier base points
        self._schreier_sims()
        self._elements = None
        self._classes = None

    # -- Schreier-Sims (deterministic, adequate for n <= ~30, |G| <= ~10^7)
    def _orbit_transversal(self, gens, b):
        n = self.n
        t = {b: identity(n)}
        queue = [b]
        while queue:
            x = queue.pop()
            for g in gens:
                y = g[x]
                if y not in t:
                    t[y] = mul(g, t[x])
                    queue.append(y)
        return t

    def _sift(self, g, start=0):
        """Return (residue, level) after sifting g through levels >= start."""
        for lvl in range(start, len(self._base)):
            b = self._base[lvl]
            y = g[b]
            t = self._transversals[lvl]
            if y not in t:
                return g, lvl
            g = mul(inverse(t[y]), g)
        return g, len(self._base)

    def _level_gens(self, lvl):
        """Strong generators fixing base[0..lvl-1]: those stored at levels >= lvl."""
        out = []
        for l in range(lvl, len(self._sgs_levels)):
            out.extend(self._sgs_levels[l])
        return out

    def _rebuild_transversals(self):
        for lvl in range(len(self._base)):
            self._transversals[lvl] = self._orbit_transversal(self._level_gens(lvl), self._base[lvl])

    def _schreier_sims(self):
        n = self.n
        gens = list(self.generators)
        if not gens:
            return
        first = next(i for g in gens for i in range(n) if g[i] != i)
        self._base = [first]
        self._sgs_levels = [list(gens)]
        self._transversals = [{}]
        self._rebuild_transversals()
        # fixed-point iteration: sift all Schreier generators at all levels
        changed = True
        while changed:
            changed = False
            for lvl in range(len(self._base)):
                t = self._transversals[lvl]
                gens_l = self._level_gens(lvl)
                for y, u in list(t.items()):
                    for s in gens_l:
                        w = t[s[y]]
                        sch = mul(inverse(w), mul(s, u))
                        if sch == identity(n):
                            continue
                        h, l2 = self._sift(sch, lvl + 1)
                        if h != identity(n):
                            if l2 == len(self._base):
                                pt = next(i for i in range(n) if h[i] != i)
                                self._base.append(pt)
                                self._sgs_levels.append([])
                                self._transversals.append({})
                            self._sgs_levels[l2].append(h)
                            self._rebuild_transversals()
                            changed = True
                            break
                    if changed:
                        break
                if changed:
                    break

    # -- basic queries
    def order(self):
        o = 1
        for t in self._transversals:
            o *= len(t)
        return o

    def contains(self, g):
        g = tuple(g)
        if len(g) != self.n:
            return False
        h, _ = self._sift(g)
        return h == identity(self.n)

    def base(self):
        return list(self._base)

    def strong_generators(self):
        out = []
        for lvl in self._sgs_levels:
            for g in lvl:
                if g not in out:
                    out.append(g)
        return out

    def random_element(self, rng=random):
        n = self.n
        g = identity(n)
        for t in self._transversals:
            g = mul(g, rng.choice(list(t.values())))
        return g

    def elements(self, limit=2_000_000):
        """Enumerate all elements (as a list) via the transversal product."""
        if self._elements is not None:
            return self._elements
        if self.order() > limit:
            raise ValueError(f"group order {self.order()} exceeds enumeration limit {limit}")
        elts = [identity(self.n)]
        for t in reversed(self._transversals):
            tv = list(t.values())
            elts = [mul(u, e) for u in tv for e in elts]
        assert len(elts) == self.order()
        self._elements = elts
        return elts

    def classes(self, limit=2_000_000):
        if self._classes is None:
            self._classes = ConjugacyClasses(self, limit=limit)
        return self._classes

    def is_subgroup_of(self, other: "PermGroup"):
        return all(other.contains(g) for g in self.generators)

    def to_json(self):
        return {"n": self.n,
                "generators": [to_json(g) for g in self.generators],
                "order": self.order(),
                "base": [b + 1 for b in self._base],
                "strong_generators": [to_json(g) for g in self.strong_generators()]}

    @staticmethod
    def from_json(d):
        return PermGroup([from_json(g) for g in d["generators"]], n=d["n"])


def symmetric(n):
    if n <= 1:
        return PermGroup([], n=n)
    return PermGroup([from_cycles(n, [list(range(n))]), from_cycles(n, [[0, 1]])], n=n)

def alternating(n):
    if n <= 2:
        return PermGroup([], n=n)
    gens = [from_cycles(n, [[0, 1, 2]])]
    if n > 3:
        if n % 2:
            gens.append(from_cycles(n, [list(range(n))]))
        else:
            gens.append(from_cycles(n, [list(range(1, n))]))
    return PermGroup(gens, n=n)

# ------------------------------------------------------------ conjugacy classes

class ConjugacyClasses:
    """Classes C_1={1},...,C_r with representatives, sizes, power maps.

    Class identification is by dictionary lookup on the enumerated group
    (the pipeline's class labels always come with this lookup; conjugators
    are recorded where the certificate needs them, see the verifier specification header).
    """

    def __init__(self, G: PermGroup, limit=2_000_000):
        self.G = G
        n = G.n
        elts = G.elements(limit=limit)
        index = {g: i for i, g in enumerate(elts)}
        cls = [-1] * len(elts)
        reps, sizes = [], []
        gens = G.generators
        # identity first
        order_ = [index[identity(n)]] + [i for i in range(len(elts)) if i != index[identity(n)]]
        for start in order_:
            if cls[start] != -1:
                continue
            k = len(reps)
            reps.append(elts[start])
            cls[start] = k
            queue = [elts[start]]
            cnt = 1
            while queue:
                x = queue.pop()
                for g in gens:
                    y = conj(g, x)
                    j = index[y]
                    if cls[j] == -1:
                        cls[j] = k
                        cnt += 1
                        queue.append(y)
            sizes.append(cnt)
        self.reps = reps
        self.sizes = sizes
        self.r = len(reps)
        self._index = index
        self._cls = cls
        self.orders = [order(g) for g in reps]
        self.cycle_types = [cycle_type(g) for g in reps]
        self.exponent = 1
        for o in self.orders:
            self.exponent = self.exponent * o // gcd(self.exponent, o)
        self.inverse_class = [self.class_of(inverse(g)) for g in reps]
        self._power_cache = {}
        # members by class (lazy)
        self._members = None

    def class_of(self, g):
        return self._cls[self._index[tuple(g)]]

    def members(self, k):
        if self._members is None:
            m = [[] for _ in range(self.r)]
            for g, i in self._index.items():
                m[self._cls[i]].append(g)
            self._members = m
        return self._members[k]

    def power_map(self, k, t):
        """class of g_k^t."""
        key = (k, t % self.orders[k])
        if key not in self._power_cache:
            self._power_cache[key] = self.class_of(power(self.reps[k], key[1]))
        return self._power_cache[key]

    def centralizer_order(self, k):
        return self.G.order() // self.sizes[k]

    def class_matrix(self, i):
        """M_i = (a_{ijk})_{j,k}, a_{ijk} = #{x in C_i: x^{-1} g_k in C_j}."""
        r = self.r
        M = [[0] * r for _ in range(r)]
        for x in self.members(i):
            xi = inverse(x)
            for k, gk in enumerate(self.reps):
                j = self.class_of(mul(xi, gk))
                M[j][k] += 1
        return M

    def to_json(self):
        return {"r": self.r,
                "representatives": [to_json(g) for g in self.reps],
                "sizes": self.sizes,
                "orders": self.orders,
                "cycle_types": [list(c) for c in self.cycle_types],
                "exponent": self.exponent,
                "inverse_class": self.inverse_class,
                "power_maps": [[self.power_map(k, t) for t in range(self.orders[k])] for k in range(self.r)]}