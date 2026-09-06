"""The direct route to a Frobenius class.

Computes the roots of f in a residue field, reads off the Frobenius
permutation, aligns the numbering, and conjugates: this always determines the
class exactly, at a cost that grows with the degree of the residue field.  It
is the fallback when the cheaper tests leave several candidates.
"""
from __future__ import annotations
from .resolvent import roots_at, alignment_single_level, coset_transversal
from .padic import frobenius_perm
from .perm import mul, inverse, symmetric, cycle_type

class DirectRoute:
    def __init__(self, G, cl, f, F_chain, c_chain, R_chain, policy, seed=0):
        self.G, self.cl, self.f = G, cl, f
        self.F, self.c, self.R = F_chain, c_chain, R_chain
        self.T = coset_transversal(symmetric(G.n), G)
        self.policy = policy
        self.seed = seed

    def excluded(self, ell):
        from .resolvent import divides_disc
        return self.policy.Delta % ell == 0 or divides_disc(self.R, ell)

    def frobenius_class(self, ell):
        k = self.policy.k(ell) if self.excluded(ell) else 1
        ring, roots, _ = roots_at(self.f, ell, k, seed=self.seed)
        assert len(roots) == self.G.n
        phi = frobenius_perm(ring, roots)
        surv = alignment_single_level(self.F, self.c, self.T, roots, ring)
        if len(surv) != 1:
            raise RuntimeError(f"alignment at {ell} left {len(surv)} survivors")
        rho = self.T[surv[0]]
        sigma = mul(inverse(rho), mul(phi, rho))
        assert self.G.contains(sigma)
        return self.cl.class_of(sigma), {"k": k, "survivor": surv[0], "cycle_type": list(cycle_type(sigma))}