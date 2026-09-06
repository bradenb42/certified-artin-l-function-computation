"""The certificate format shared with the descent step.

The input certificate records f, the numbering of the roots, generators of the
Galois group, and the descent chain with its invariants and resolvents.  The
output certificate records those verbatim together with the results of every
stage.  Permutations are 1-indexed image lists; composition is
(p*q)(i) = p(q(i)), matching the action on polynomials in the roots.
"""
from __future__ import annotations
import json
import os
from fractions import Fraction

from .perm import PermGroup, from_json, to_json

CERT_VERSION = "artin-cert-1"

class DescentCertificate:
    def __init__(self, data: dict, source: str):
        self.data = data
        self.source = source
        self.f = [int(c) for c in data["f"]]
        if self.f[-1] != 1:
            raise ValueError("f must be monic")
        self.n = len(self.f) - 1
        gens = [from_json(g) for g in data["group"]["generators"]]
        for g in gens:
            if len(g) != self.n or sorted(g) != list(range(self.n)):
                raise ValueError("generator is not a permutation of 1..n")
        self.G = PermGroup(gens, n=self.n)
        if "order" in data["group"] and int(data["group"]["order"]) != self.G.order():
            raise ValueError("recorded group order disagrees with the generators")
        self.verified = bool(data.get("verified", False))
        self.numbering = data.get("numbering")
        self.chain = data.get("chain")

    @staticmethod
    def load(path):
        with open(path) as fh:
            return DescentCertificate(json.load(fh), source=os.path.abspath(path))

    @staticmethod
    def from_generators(f, generators_1indexed, numbering=None):
        data = {"version": "inline", "f": list(f), "n": len(f) - 1,
                "group": {"generators": [list(g) for g in generators_1indexed]},
                "verified": False, "numbering": numbering}
        return DescentCertificate(data, source="inline (no descent certificate: G is an unverified input)")


def _json_default(o):
    if isinstance(o, Fraction):
        return str(o)
    if isinstance(o, tuple):
        return list(o)
    raise TypeError(f"not serializable: {type(o)}")

def dump_json(obj, path):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=_json_default)

def load_json(path):
    with open(path) as fh:
        return json.load(fh)