"""The archimedean stage.
"""
from .archimedean import conjugation_class, archimedean_data, permutation_check, isolate_real_roots, real_root_count
from .filtration import HardFailure

def run_archimedean(G, cl, table, f, assigner, log=print):
    k_c, r, used = conjugation_class(G, cl, f, assigner.stage, log)
    iso = isolate_real_roots(f)
    if len(iso) != r:
        raise HardFailure("isolating intervals do not match the Sturm count")
    data = archimedean_data(table, cl, k_c, r, G.n)
    if not all(d["a+b=chi(1)"] for d in data):
        raise HardFailure("a + b != chi(1)")
    pc = permutation_check(table, cl, data, G, r)
    if not pc["ok"]:
        raise HardFailure(f"parity classification vs real root count fails: {pc}")
    if pc["regular_sum_chi1_chi(c)"] != (G.order() if r == G.n else 0):
        raise HardFailure("regular representation sum fails")
    log(f"archimedean: a+b = chi(1) on all {table.r} characters; sum n_chi b_chi = {pc['sum_n_chi_b']} = complex pairs; parities {[d['parity'] for d in data]}; odd 2-dim: {[d['chi'] for d in data if d['odd_two_dimensional']]}")
    return {"real_roots": r, "isolating_intervals": [[str(a), str(b)] for a, b in iso], "class_of_c": k_c + 1,
            "separation": used, "characters": data, "permutation_check": pc}