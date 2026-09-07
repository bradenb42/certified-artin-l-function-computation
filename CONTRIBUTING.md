# Contributing

This repository contains both the Artin L-function computation pipeline and an independent certificate checker. Changes should preserve that separation and keep generated run output out of the source tree.

## Development setup

Use Python 3.10 or later and install the project with its test dependencies:

```bash
python -m pip install -e '.[test]'
python -m pytest -q
```

The GitHub Actions workflow runs the test suite on pull requests and pushes.

## Repository boundaries

- `artin/` contains the computation pipeline, mathematical routines, validation helpers, and stage drivers.
- `artinverify/` contains the independent certificate checker. Keep it independent of `artin/`; do not import prover-side implementation details into the checker.
- `tests/` contains regression and end-to-end tests.
- `examples/` contains committed examples and recorded fixtures intended to be reproducible.
- `docs/` contains design, verification, validation, performance, scope, and theory documentation.
- `runs/` and `out/` are generated working directories and are intentionally ignored.

See `docs/repository_structure.md` for a fuller map.

## Change discipline

Keep pull requests focused. Organization-only changes should not silently alter mathematical behavior, certificate semantics, or recorded validation evidence.

When changing a computational stage, update the corresponding tests and any affected documentation. When changing the certificate format or checker logic, update both `docs/verification.md` and tests that exercise rejection of altered certificates.

Avoid committing local caches, virtual environments, coverage output, editor metadata, or generated run directories.

## Style and files

Follow `.editorconfig` and `.gitattributes`. Python uses four-space indentation and UTF-8 text files use LF line endings.

Prefer descriptive module and test names that match the mathematical stage or responsibility they implement. Keep generated data distinguishable from source code and documentation.
