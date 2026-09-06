# Certified Artin L-function computation

Compute Artin L-functions from a polynomial.

Given a monic separable polynomial `f` in `Z[x]` and generators of its Galois
group as permutations of the roots, this package computes, for every
irreducible character of that group:

- the conductor,
- the Euler factors at the ramified primes and at every prime up to a bound,
- the Gamma factor and the sign of the functional equation,
- the root number, exactly for self-dual characters,

together with a certificate that a separate program checks without redoing the
work.

Everything is exact where exactness is possible: character values live in
cyclotomic fields as normal forms, resolvents are recovered as integer
polynomials, ramification data come from maximal orders, and real root counts
come from Sturm sequences. The two places where floating point enters — the
Gauss sums behind the root numbers and the approximate functional equation —
are the two places where the result is a check rather than a value, and each
carries an explicit error bound.

## Install

Python 3.10 or later.

```bash
pip install -r requirements.txt      # sympy, mpmath, numpy
```

## Use

```bash
python -m artin.run --f "x^5-x-1" --gens "[[2,3,4,5,1],[2,1,3,4,5]]" --run-dir runs/x5
python -m artin.verify runs/x5                  # re-check the run in place
python -m artinverify.verify runs/x5/CERT.json  # check the certificate independently
```

Generators are 1-indexed image lists: `[2,3,4,5,1]` is the cycle sending root 1
to root 2, and so on. Instead of `--f`/`--gens` you can pass a certificate
from a descent computation:

```bash
python -m artin.run --descent examples/descent_x5-x-1.json --run-dir runs/x5
python -m artin.run --config examples/config_example.json
```

A run writes its configuration first, then one file per stage:

| file | contents |
| --- | --- |
| `group.json` | generators, order, classes with power maps |
| `chartable.json` | the character table and its certificate |
| `schur.json`, `models/` | Schur index data and matrix models |
| `ramified.json` | ramified primes with witnesses and decompositions |
| `precision.json` | the precision policy and every consultation |
| `classes.json` | the Frobenius class at each unramified prime |
| `local.json` | decomposition and inertia groups, matchings |
| `conductors.json` | filtrations and conductor exponents |
| `euler.json` | Euler factors and the identity system |
| `archimedean.json` | Gamma factors, parities, complex conjugation |
| `rootnumbers.json` | local and global root numbers |
| `analytic.json`, `falsifier.json` | functional-equation defects and the perturbation experiment |
| `CERT.json` | all of the above in one self-contained file |
| `timings.json`, `log.txt` | per-stage times and a log |

## Example

```python
from artin.run import run_pipeline

run_pipeline({
    "f": [-1, -1, 0, 0, 0, 1],  # x^5 - x - 1, ascending
    "generators": [[2, 3, 4, 5, 1], [2, 1, 3, 4, 5]],
    "run_dir": "runs/x5",
    "options": {"class_bound_X": 200},
})
```

## Options

Passed under `"options"`; the defaults are in `artin/run.py`.

| option | default | meaning |
| --- | --- | --- |
| `class_bound_X` | 200 | assign Frobenius classes at unramified primes up to this bound |
| `precision_mode` | `"sharp"` | `"sharp"` or `"conservative"` precision policy |
| `precision_check` | `True` | recompute everything at double precision and compare |
| `class_confirmation` | `True` | check the classes against an explicit splitting field when the group is small |
| `root_numbers` | `True` | compute root numbers (Gauss sums; cost grows with the ramified primes) |
| `fe_test` | `True` | run the functional-equation test and the perturbation experiment |
| `compute_models` | `True` | build matrix models of the characters |
| `class_enumeration_limit` | 2000000 | refuse to enumerate groups larger than this |

## How it is checked

Beyond the internal identities enforced at every prime, the output is compared
with independent computations:

- abelian fields of conductor up to 500 against class field theory — 530 fields, every character agreeing;
- a sample of degree at most 5 against explicitly constructed Galois closures;
- odd two-dimensional characters against weight-one newforms computed from binary quadratic forms — 5100 coefficients, no mismatch;
- the residue of a subfield's Dedekind zeta function computed from `f` alone against the product of L-values from the character data;
- the certificate checker's rejection of certificates in which a single datum has been altered.

See [docs/validation.md](docs/validation.md).

## Documentation

- [docs/design.md](docs/design.md) — what each stage does and why
- [docs/verification.md](docs/verification.md) — the certificate and the checker
- [docs/validation.md](docs/validation.md) — the cross-checks and their results
- [docs/performance.md](docs/performance.md) — where the time goes
- [docs/scope.md](docs/scope.md) — what is covered and what is not
- [docs/theory/](docs/theory/README.md) — the derivations and proofs implemented by the code

## Tests

```bash
python -m pip install -e '.[test]'
python -m pytest -q
```

Forty-nine tests covering the group and field arithmetic, each stage, the
end-to-end runs, the certificate checker and its rejection of altered
certificates. The two analytic tests take about two minutes; the rest run in
well under one.

## Layout

```text
artin/            the computation, one module per piece of mathematics
  stages_*.py     the per-stage drivers used by artin/run.py
artinverify/      the independent certificate checker
examples/         runnable scripts, sample certificates, recorded results
docs/             design, verification, validation, performance, scope
  theory/         derivations and proofs underlying each computational stage
tests/            the test suite
```

## Licence

[MIT](LICENSE).
