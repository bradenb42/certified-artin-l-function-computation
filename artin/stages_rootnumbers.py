"""Root numbers: local factors, the global product and the self-dual cases.

Each local factor is computed at every tamely ramified prime and multiplied
with the archimedean factor.  For self-dual characters the local factors lie
in a finite set (fourth roots of unity in the orthogonal case, signs in the
symplectic one), which makes the global root number exact; the corresponding
predictions are checked and a mismatch stops the run.
"""
from __future__ import annotations
from math import gcd
import cmath

from .cyclo import Cyc
from .rootnumbers import TameRootNumbers, lcm
from .filtration import HardFailure
from .perm import from_json

def fourth_root_index(x):
    """k with x = i^k, or None."""
    for k in range(4):
        if x == Cyc.root_of_unity(x.e, 4, k) if x.e % 4 == 0 else (x == 1 if k == 0 else x == -1 if k == 2 else False):
            return k
    return None

def run_root_numbers(G, cl, table, objects, local_json, AR, log=print):
    wild = [ell for ell, r in local_json["ramified"].items() if r["status"] != "ok"]
    tame = {}
    for ell, (rec, local, D, I, Fr) in objects.items():
        if I.order() == 1:
            continue
        tame[ell] = TameRootNumbers(local, rec, table, cl, from_json(rec["matching"]), log)
    b_of = {d["chi"]: d["b"] for d in AR["characters"]}
    out = {"wild_primes": wild, "characters": []}
    conj_row = table.galois_action[table.e - 1] if table.e > 1 else list(range(table.r))
    Wnum = {}
    for nu in range(table.r):
        real, fs = table.is_real[nu], table.indicator[nu]
        entry = {"chi": nu + 1, "degree": table.degrees[nu], "real": real, "frobenius_schur": fs, "b": b_of[nu + 1],
                 "epsilon_infinity": f"i^{-b_of[nu+1] % 4}"}
        if wild:
            entry["W"] = None
            entry["status"] = f"unknown: wild primes {wild} (wild root numbers are outside the implemented scope)"
            out["characters"].append(entry)
            continue
        Wc = cmath.exp(-1j * cmath.pi / 2 * b_of[nu + 1])
        k4 = (-b_of[nu + 1]) % 4          # exponent of i in W when all local factors are in mu_4
        exact_mu4 = True
        loc = {}
        for ell, TR in tame.items():
            val, h, det = TR.epsilon(nu, log)
            if abs(abs(val) - 1) > 1e-8:
                raise HardFailure(f"|epsilon_{ell}(chi_{nu+1})| != 1 (numerically {abs(val)})")
            # finite candidate set for self-dual characters: mu_4 (orthogonal, Deligne) or +-1 (symplectic)
            k = None
            for kk in range(4):
                if abs(val - 1j ** kk) < 1e-7:
                    k = kk
                    break
            loc[str(ell)] = {"complex": [round(val.real, 12), round(val.imag, 12)], "sqrt_ell_power": h,
                             "in_mu4": k is not None, "i_exponent": k, "orbits": det}
            Wc *= val
            if k is None:
                exact_mu4 = False
            else:
                k4 = (k4 + k) % 4
        entry["local"] = loc
        entry["W_complex"] = [round(Wc.real, 12), round(Wc.imag, 12)]
        Wnum[nu] = Wc
        entry["modulus_one"] = abs(abs(Wc) - 1) < 1e-9
        sc = {}
        if real and fs == 1:
            if not exact_mu4:
                raise HardFailure(f"orthogonal chi_{nu+1}: a local root number is not a fourth root of unity")
            entry["W_exact"] = f"i^{k4}"
            sc["orthogonal_W=1"] = (k4 == 0)
            if table.degrees[nu] == 1:
                sc["quadratic_or_trivial_W=1"] = (k4 == 0)
        if real and fs == -1:
            if not exact_mu4 or any(v["i_exponent"] % 2 for v in loc.values()):
                raise HardFailure(f"symplectic chi_{nu+1}: a local root number is not +-1")
            entry["W_exact"] = f"i^{k4}"
            sc["symplectic_W=+-1"] = (k4 % 2 == 0)
        for name, ok in sc.items():
            if not ok:
                raise HardFailure(f"shortcut {name} fails for chi_{nu+1}: W = i^{k4}")
        entry["shortcuts"] = sc
        entry["status"] = "ok" if (sc or True) else ""
        entry["exact"] = exact_mu4
        out["characters"].append(entry)
    # W(conj chi) = conj W(chi) (numerically for non-self-dual characters)
    for nu, W in Wnum.items():
        nb = conj_row[nu]
        if nb in Wnum and abs(Wnum[nb] - W.conjugate()) > 1e-9:
            raise HardFailure(f"W(conj chi_{nu+1}) != conj W(chi_{nu+1})")
    done = [c for c in out["characters"] if c.get("W_complex") is not None]
    if done:
        log(f"root numbers: {len(done)} characters; exact (mu_4) for {[c['chi'] for c in done if c['exact']]}; "
            f"shortcuts on {[c['chi'] for c in done if c['shortcuts']]}; values {[(c['chi'], c['W_complex']) for c in done]}")
    else:
        log(f"root numbers: not computed (wild primes {wild})")
    return out