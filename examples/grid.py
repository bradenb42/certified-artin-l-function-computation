"""Cost decomposition over a grid of inputs: runs the pipeline, records the per-stage
timings and the statistics (rational classes needing cyclotomic refinement, undetermined
filtrations, largest wild inertia (lcm of the ramification indices at wild primes), size of the largest Brauer
induction datum), and packages certificates, conductor and Frobenius
tables and verification defects.  Regenerates every table:

    python examples/grid.py <outdir> [names...]
"""
import sys, os, json, time, shutil
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from artin.run import run_pipeline
from artin.certificate import load_json
from artin.ramified import discriminant, factorint

GRID = {
 # name: (f ascending, generators 1-indexed, options)
 "d3_S3_x3-x-1": ([-1, -1, 0, 1], [[2, 3, 1], [2, 1, 3]], {}),
 "d4_S4_x4-x-1": ([-1, -1, 0, 0, 1], [[2, 3, 4, 1], [2, 1, 3, 4]], {}),
 "d4_D4_x4-2": ([-2, 0, 0, 0, 1], [[2, 3, 4, 1], [1, 4, 3, 2]], {}),
 "d5_C5_cond11": ([1, 3, -3, -4, 1, 1], [[2, 3, 4, 5, 1]], {}),
 "d5_S5_x5-x-1": ([-1, -1, 0, 0, 0, 1], [[2, 3, 4, 5, 1], [2, 1, 3, 4, 5]], {"fe_test": False}),
 "d5_A5_x5+20x+16": ([16, 20, 0, 0, 0, 1], [[2, 3, 1, 4, 5], [2, 3, 4, 5, 1]], {}),
 "d6_D6_x6-2": ([-2, 0, 0, 0, 0, 0, 1], [[2, 3, 4, 5, 6, 1], [1, 6, 5, 4, 3, 2]], {"class_confirmation": False}),
 "d6_S6_x6-x-1": ([-1, -1, 0, 0, 0, 0, 1], [[2, 3, 4, 5, 6, 1], [2, 1, 3, 4, 5, 6]], {"class_confirmation": False}),
 "d7_S7_x7-x-1": ([-1, -1, 0, 0, 0, 0, 0, 1], [[2, 3, 4, 5, 6, 7, 1], [2, 1, 3, 4, 5, 6, 7]], {"class_confirmation": False}),
 "d8_S8_x8-x-1": ([-1, -1, 0, 0, 0, 0, 0, 0, 1], [[2, 3, 4, 5, 6, 7, 8, 1], [2, 1, 3, 4, 5, 6, 7, 8]], {"class_confirmation": False}),
}

def stats(d):
    C = load_json(os.path.join(d, "classes.json")); L = load_json(os.path.join(d, "local.json"))
    CJ = load_json(os.path.join(d, "conductors.json")); G = load_json(os.path.join(d, "group.json"))
    methods = {m: sum(1 for r in C["primes"] if r["method"] == m) for m in ("A5", "A6", "A7")}
    refine = sum(1 for orb in C["rational_classes"] if len(orb) >= 2)
    undetermined = sum(1 for p in CJ["primes"].values() if not p["filtration_candidates"]["unique"])
    wild = [r for r in L["ramified"].values() if r["status"].startswith("wild")]
    largest_wild_E = max([r["e_lcm"] for r in wild], default=0)
    return {"|G|": G["order"], "r": G["classes"]["r"], "rational_classes_needing_refinement": refine,
            "primes_by_method": methods, "undetermined_filtrations": undetermined,
            "largest_wild_inertia_lcm_e": largest_wild_E, "wild_primes": [r["ell"] for r in wild],
            "largest_brauer_datum": "n/a"}

def main():
    out = sys.argv[1]; names = sys.argv[2:] or list(GRID)
    os.makedirs(out, exist_ok=True)
    results = load_json(os.path.join(out, "grid_results.json")) if os.path.exists(os.path.join(out, "grid_results.json")) else {}
    for name in names:
        f, gens, extra = GRID[name]
        d = os.path.join(out, "runs", name)
        shutil.rmtree(d, ignore_errors=True)
        D = abs(discriminant(f))
        big = max(factorint(D)) if D > 1 else 1
        opts = {"class_bound_X": 200, "precision_check": False, "local_unramified_check": 2, "quiet": True,
                "root_numbers": big <= 5000, "fe_test": big <= 5000 and len(f) - 1 <= 5, "fe_cap": 2000}
        opts.update(extra)
        t = time.time()
        try:
            run_pipeline({"f": f, "generators": gens, "run_dir": d, "options": opts})
            rec = {"f": f, "n": len(f) - 1, "disc": D, "largest_ramified_prime": big, "options": opts,
                   "timings": load_json(os.path.join(d, "timings.json")), "stats": stats(d), "total": round(time.time() - t, 1), "status": "ok"}
            # tables
            CJ = load_json(os.path.join(d, "conductors.json")); C = load_json(os.path.join(d, "classes.json"))
            rec["conductor_table"] = [{"chi": c["chi"], "exponents": c["exponents"], "partial": c["partial_conductor"], "unknown_at": c["unknown_at"]} for c in CJ["conductors"]]
            rec["frobenius_table"] = [[r["ell"], r["class"], r["method"]] for r in C["primes"]]
            an = os.path.join(d, "analytic.json")
            if os.path.exists(an):
                AN = load_json(an)
                rec["defects"] = [{"chi": c["chi"], "tests": [{"t": t_["t"], "defect": t_["defect"], "bound": t_["bound"]} for t_ in c["tests"]]} for c in AN["characters"] if c.get("tests")]
            shutil.copy(os.path.join(d, "CERT.json"), os.path.join(out, f"CERT_{name}.json"))
        except Exception as ex:
            rec = {"f": f, "n": len(f) - 1, "status": f"failed: {ex}", "total": round(time.time() - t, 1)}
        results[name] = rec
        print(name, rec["status"], rec["total"], rec.get("timings"), flush=True)
        json.dump(results, open(os.path.join(out, "grid_results.json"), "w"), indent=1)
    render(results, out)

def render(results, out):
    stages = ["numbering_and_group_placement", "character_table", "ramified_primes", "precision_check",
              "class_assignment", "class_confirmation", "local_descents",
              "filtrations_and_conductors", "euler_factors", "archimedean", "root_numbers", "functional_equation"]
    lines = ["# Cost decomposition over the grid (seconds per stage)\n", "| input | n | |G| | log10 disc | " + " | ".join(s.replace("_", " ") for s in stages) + " | total |", "|" + "---|" * (len(stages) + 5)]
    import math
    for name, r in results.items():
        if r["status"] != "ok":
            lines.append(f"| {name} | {r['n']} | | | {r['status']} |"); continue
        T = r["timings"]
        lines.append(f"| {name} | {r['n']} | {r['stats']['|G|']} | {math.log10(r['disc']):.1f} | " + " | ".join(str(T.get(s, "")) for s in stages) + f" | {r['total']} |")
    lines += ["\n# Statistics\n", "| input | |G| | r | rational classes needing refinement | primes by method | undetermined filtrations | largest wild inertia (lcm e) | wild primes | largest Brauer datum |", "|---|---|---|---|---|---|---|---|---|"]
    for name, r in results.items():
        if r["status"] != "ok": continue
        s = r["stats"]
        lines.append(f"| {name} | {s['|G|']} | {s['r']} | {s['rational_classes_needing_refinement']} | {s['primes_by_method']} | {s['undetermined_filtrations']} | {s['largest_wild_inertia_lcm_e']} | {s['wild_primes']} | {s['largest_brauer_datum']} |")
    lines += ["\n# Conductor tables (exponents at the tame ramified primes; 'unknown_at' lists wild primes)\n"]
    for name, r in results.items():
        if r["status"] != "ok": continue
        lines.append(f"\n## {name}\n")
        for c in r["conductor_table"]:
            lines.append(f"- chi_{c['chi']}: {c['exponents']} partial conductor {c['partial']} unknown at {c['unknown_at']}")
    lines += ["\n# Verification defects\n"]
    for name, r in results.items():
        if r.get("defects"):
            lines.append(f"\n## {name}\n")
            for c in r["defects"]:
                lines.append(f"- chi_{c['chi']}: " + ", ".join(f"t={t['t']}: defect {t['defect']:.2e} bound {t['bound']:.2e}" for t in c["tests"]))
    open(os.path.join(out, "GRID.md"), "w").write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main()
