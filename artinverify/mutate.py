"""Single-datum alterations of a certificate, for testing the checker.
"""
import copy, json, sys

def mutations(cert):
    out = {}
    S = cert["sections"]
    # 1. filtration: raise one i_L value by 1 (delta adjusted consistently)
    c = copy.deepcopy(cert)
    for ell, rec in c["sections"]["conductors"]["primes"].items():
        fil = rec["filtration"]
        if fil["i_L"]:
            k = next(iter(fil["i_L"]))
            fil["i_L"][k] += 1
            fil["delta"] += 1
            out["filtration"] = c
            break
    # 2. conductor exponent of one nontrivial character at one prime, +1 (partial conductor adjusted)
    c = copy.deepcopy(cert)
    CJ = c["sections"]["conductors"]
    for ell, rec in CJ["primes"].items():
        ce = rec["conductor_exponents"][1]
        ce["f_ell"] += 1
        for cc in CJ["conductors"]:
            if cc["chi"] == ce["chi"]:
                cc["exponents"][ell] += 1
                cc["partial_conductor"] *= int(ell)
        out["conductor"] = c
        break
    # 3. Frobenius class at one unramified prime: another class (same cycle type if possible)
    c = copy.deepcopy(cert)
    recs = c["sections"]["classes"]["primes"]
    cts = c["sections"]["group"]["classes"]["cycle_types"]
    for rec in recs:
        k = rec["class"] - 1
        same = [j for j in range(len(cts)) if cts[j] == cts[k] and j != k]
        newk = same[0] if same else (k + 1) % len(cts)
        rec["class"] = newk + 1
        rec["mutated_from"] = k + 1
        out["frobenius"] = c
        break
    # 4. Euler factor: negate the T-coefficient of P for a character with deg P >= 1 at a ramified prime
    c = copy.deepcopy(cert)
    done = False
    for ell, rec in c["sections"]["euler"]["ramified"].items():
        for chi, cc in rec["characters"].items():
            if len(cc["P"]) >= 2:
                cc["P"][1]["c"] = [(-x if isinstance(x, int) else str(-int(x))) for x in cc["P"][1]["c"]]
                done = True
                break
        if done:
            break
    out["euler"] = c
    # 5. root number: negate one local numerator (and W) for the first character with a local factor
    c = copy.deepcopy(cert)
    for cc in c["sections"]["rootnumbers"]["characters"]:
        if cc.get("local"):
            ell, loc = next(iter(cc["local"].items()))
            loc["complex"] = [-loc["complex"][0], -loc["complex"][1]]
            if loc["i_exponent"] is not None:
                loc["i_exponent"] = (loc["i_exponent"] + 2) % 4
            cc["W_complex"] = [-cc["W_complex"][0], -cc["W_complex"][1]]
            if "W_exact" in cc:
                cc["W_exact"] = "mutated"
            for a in c["sections"].get("analytic", {}).get("characters", []):
                if a.get("chi") == cc["chi"] and "W" in a:
                    a["W"] = [-a["W"][0], -a["W"][1]]
            out["rootnumber"] = c
            break
    return out

if __name__ == "__main__":
    cert = json.load(open(sys.argv[1]))
    for name, c in mutations(cert).items():
        p = sys.argv[1].replace(".json", f".mut_{name}.json")
        json.dump(c, open(p, "w"))
        print(p)
