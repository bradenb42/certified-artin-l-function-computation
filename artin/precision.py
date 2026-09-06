"""The precision policy: the single authority on working precision.

One Policy object is created per run.  Every stage registers the pairs
(invariant, index) it will use, which can only raise precisions, and obtains
its working precision at each prime from the policy; working below it raises
PolicyViolation.  Every consultation is logged so that a verifier can confirm
that no stage worked below the precision its own data require.

Two modes are available: a conservative one, whose bound covers every test in
advance, and a sharp one, which applies the bound test by test and needs
about a factor of the index fewer digits.
"""
from __future__ import annotations
import math

class PolicyViolation(Exception):
    pass

def mu(m):
    return 1 if m <= 1 else m * (m - 1)

def ilog(base, x):
    """floor(log_base x) for integers x >= 1."""
    k, y = 0, base
    while y <= x:
        y *= base
        k += 1
    return k

class Policy:
    """mode "conservative": the precision policy literally, M_F = (2B_F+2)^{mu(m)} covers T1-T5 for
    every registered pair without knowing any resolvent.
    mode "sharp": the same theorem applied test by test (the precision policy and (4)):
    T1/T2 need ell^k > (2B_F+2)^m; T3 (pairwise comparison of resolvent roots) is not used
    by the pipeline, every comparison being against a recorded integer; T4 needs k > v_ell(Delta);
    T5 (root counts of a resolvent R) needs k > v_ell(disc R), read from the exact disc R
    registered by register_disc once R is known.  Both modes are rigorous; sharp gives
    O(m) instead of O(m^2) digits."""
    def __init__(self, f, disc, root_bound=None, refined=True, mode="conservative"):
        self.f = list(f)
        self.n = len(f) - 1
        self.Delta = abs(int(disc))
        # Cauchy bound R = 1 + max |a_i|
        self.R = root_bound if root_bound is not None else 1 + max(abs(int(c)) for c in f[:-1])
        self.refined = refined     # use v_ell(Delta) when known
        self.pairs = []            # dicts: label, degree, norm1, m, B_F, M_F
        self.Mstar = 1
        self.consultations = {}    # ell -> k issued
        self.history = []          # (label, M* before, M* after)
        self.mode = mode
        self.discs = []            # (label, |disc R|) for T5 in sharp mode
        self.extra = {}            # T5 precisions raised beyond k(ell) by count_Zl_roots_policy
        self.log = []              # every consultation with the number of pairs registered at the time

    # ---- registration (may happen at any time; only increases precisions, the precision policy
    def register(self, F, m, label=None):
        B = F.height_bound(self.R)
        MF = (2 * B + 2) ** (mu(m) if self.mode == "conservative" else m)
        rec = {"label": label or F.label, "degree": F.degree(), "norm1": F.norm1(), "m": m,
               "B_F": B, "log2_M_F": MF.bit_length() - 1}
        self.pairs.append(rec)
        before = self.Mstar
        self.Mstar = max(self.Mstar, MF)
        self.history.append((rec["label"], before.bit_length(), self.Mstar.bit_length()))
        return rec

    def register_disc(self, label, disc):
        """Sharp mode, T5: record the exact discriminant of a resolvent whose local root
        counts will be taken (the precision policy: precision > v_ell(disc R))."""
        d = abs(int(disc))
        if d == 0:
            raise ValueError("resolvent not squarefree")
        self.discs.append((label, d))

    @property
    def M(self):
        return self.Mstar * self.Delta

    # ---- the authority
    def k(self, ell, v_disc=None):
        """Precision at ell.  If refined and v_ell(Delta) is supplied (or ell does not divide
        Delta), the sharper bound is used; otherwise k = floor(log_ell M) + 1."""
        if self.refined and (v_disc is not None or self.Delta % ell):
            v = 0 if v_disc is None else v_disc
            k = max(ilog(ell, self.Mstar) + 1, v + 1)
        else:
            k = ilog(ell, self.M) + 1
        if self.mode == "sharp":
            for _, d in self.discs:
                vd = 0
                while d % ell == 0:
                    d //= ell
                    vd += 1
                k = max(k, vd + 1)
        if v_disc is None and self.Delta % ell == 0 and self.refined:
            assert ell ** k > self.M
        else:
            assert ell ** k > self.Mstar and k > (v_disc or 0)
        prev = self.consultations.get(ell)
        self.consultations[ell] = max(prev, k) if prev else k
        self.log.append({"ell": ell, "k": k, "pairs": len(self.pairs), "v_disc": v_disc})
        return k

    def check(self, ell, k_used, v_disc=None):
        need = self.k(ell, v_disc)
        if k_used < need:
            raise PolicyViolation(f"precision {k_used} at {ell} below policy {need}")
        return True

    def bits(self, ell, v_disc=None):
        return self.k(ell, v_disc) * math.log2(ell)

    def to_json(self):
        return {"f": self.f, "Delta_bits": self.Delta.bit_length(), "root_bound_R": self.R,
                "refined": self.refined, "mode": self.mode,
                "registered_discs": [{"label": l, "bits": d.bit_length()} for l, d in self.discs],
                "pairs": self.pairs, "log2_Mstar": self.Mstar.bit_length() - 1, "log2_M": self.M.bit_length() - 1,
                "consultations": {str(l): k for l, k in sorted(self.consultations.items())},
                "T5_extra": self.extra,
                "consultation_log": self.log,
                "history": self.history}
