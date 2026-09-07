# Repository structure

The repository separates computation, independent verification, documentation, examples, and tests so that the certificate checker is not merely a second entry point into the prover.

## Top-level map

| path | role |
| --- | --- |
| `artin/` | Artin L-function computation pipeline and supporting mathematics |
| `artinverify/` | independent certificate checker |
| `tests/` | unit, regression, certificate-rejection, and end-to-end tests |
| `examples/` | runnable examples, sample inputs/certificates, and committed recorded results |
| `docs/` | design, verification, validation, performance, scope, and mathematical derivations |
| `.github/workflows/` | continuous-integration configuration |
| `pyproject.toml` | package metadata, dependencies, and pytest configuration |
| `requirements.txt` | direct runtime dependency list for simple installs |
| `runs/`, `out/` | generated working/output directories; ignored by Git |

## `artin/`: computation side

The `artin` package owns the prover/computation side. Its modules cover the mathematical objects and stages used to construct a run and ultimately `CERT.json`. The `stages_*.py` modules are the per-stage orchestration layer used by `artin.run`.

Validation helpers such as the abelian, splitting-field, and newform cross-checks live on this side because they compare computed results against independent mathematical constructions; they are not part of the standalone certificate checker.

## `artinverify/`: checker side

`artinverify` is the independent checker for the serialized certificate. Keep this package independent of implementation details in `artin/` so an error in the computation is not automatically reproduced in verification.

If certificate fields or semantics change, update the checker, `docs/verification.md`, and the relevant tests together.

## Documentation

The documentation tree is organized by purpose:

- `design.md` describes the computational stages and architecture.
- `verification.md` describes the certificate and independent checker.
- `validation.md` records external and internal cross-checks.
- `performance.md` documents runtime behavior and costs.
- `scope.md` records supported and unsupported cases.
- `theory/` contains derivations and proofs underlying the implementation.

Repository-maintenance conventions belong in this file and `CONTRIBUTING.md`, rather than being mixed into the mathematical documents.

## Generated versus committed artifacts

A normal run writes stage JSON files, models, logs, timings, and `CERT.json` beneath a selected run directory. These should normally live under `runs/` or `out/`, both of which are ignored.

Files under `examples/` are different: they are committed intentionally when they serve as reproducible fixtures, sample certificates, or recorded validation results. Changes to such files should be treated as evidence changes and reviewed alongside the code or configuration that produced them.
