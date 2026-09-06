"""Assembly of the self-contained certificate.

Collects the per-stage records of a run directory into a single file that
embeds everything the standalone verifier needs, with no references to other
files.
"""
from __future__ import annotations
import json, os

from .certificate import load_json, dump_json

SECTIONS = ["config", "descent", "group", "chartable", "ramified", "precision", "classes", "local",
            "conductors", "euler", "archimedean", "rootnumbers", "analytic", "coefficients", "falsifier"]

def write_certificate(run_dir, out_path=None):
    cert = {"format": "artin-CERT-1", "sections": {}}
    for name in SECTIONS:
        p = os.path.join(run_dir, name + ".json")
        if os.path.exists(p):
            cert["sections"][name] = load_json(p)
    base = load_json(os.path.join(run_dir, "certificate.json"))
    cert["summary"] = {k: v for k, v in base.items() if k in ("numbering", "group", "character_table", "precision", "ramified_primes") or k in ("frobenius_classes", "local_data", "conductors", "euler_factors", "archimedean", "root_numbers", "functional_equation", "falsifier", "schur")}
    cert["input"] = base["input"]
    out_path = out_path or os.path.join(run_dir, "CERT.json")
    dump_json(cert, out_path)
    return out_path