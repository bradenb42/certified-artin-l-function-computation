# Validation

Recorded independent validation data for the Artin L-function pipeline.

Raw data:

- [`examples/results/weight_one_newforms.json`](../examples/results/weight_one_newforms.json)
- [`examples/results/abelian_class_field_checks.json`](../examples/results/abelian_class_field_checks.json)
- [`examples/results/explicit_galois_closure_checks.json`](../examples/results/explicit_galois_closure_checks.json)
- [`examples/CERT_example_C4.json`](../examples/CERT_example_C4.json) — sample self-contained certificate

## Weight-one newforms

The recorded dataset contains 17 odd two-dimensional character comparisons. Each comparison checks 300 coefficients, for 5100 coefficients total. The recorded results contain zero coefficient mismatches; every level, nebentypus, and root-number comparison is marked successful.

| D | |G| | chi | level | class number | coefficients compared | mismatches | level | nebentypus | root number |
|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|:---:|
| -23 | 6 | 3 | 23 | 3 | 300 | 0 | ok | ok | ok |
| -31 | 6 | 3 | 31 | 3 | 300 | 0 | ok | ok | ok |
| -59 | 6 | 3 | 59 | 3 | 300 | 0 | ok | ok | ok |
| -211 | 6 | 3 | 211 | 3 | 300 | 0 | ok | ok | ok |
| -239 | 6 | 3 | 239 | 15 | 300 | 0 | ok | ok | ok |
| -283 | 6 | 3 | 283 | 3 | 300 | 0 | ok | ok | ok |
| -283 | 24 | 3 | 283 | 3 | 300 | 0 | ok | ok | ok |
| -1291 | 6 | 3 | 1291 | 9 | 300 | 0 | ok | ok | ok |
| -1319 | 6 | 3 | 1319 | 45 | 300 | 0 | ok | ok | ok |
| -1327 | 6 | 3 | 1327 | 15 | 300 | 0 | ok | ok | ok |
| -1579 | 6 | 3 | 1579 | 9 | 300 | 0 | ok | ok | ok |
| -1823 | 6 | 3 | 1823 | 45 | 300 | 0 | ok | ok | ok |
| -419 | 6 | 3 | 419 | 9 | 300 | 0 | ok | ok | ok |
| -499 | 6 | 3 | 499 | 3 | 300 | 0 | ok | ok | ok |
| -643 | 6 | 3 | 643 | 3 | 300 | 0 | ok | ok | ok |
| -743 | 6 | 3 | 743 | 21 | 300 | 0 | ok | ok | ok |
| -823 | 6 | 3 | 823 | 9 | 300 | 0 | ok | ok | ok |

## Abelian class-field checks

The recorded class-field dataset contains 530 cases. Every record is marked `status: "ok"`. These are the data underlying the README statement that abelian fields of conductor up to 500 were checked against class field theory.

## Explicit Galois-closure checks

The recorded Galois-closure dataset contains 15 polynomial examples. At primes where the explicit reference construction is available and the pipeline supplies comparable local data, the recorded decomposition data are compared directly. Wild primes for which the tame local construction is not applicable, and cases whose explicit reference computation is over budget, are retained explicitly rather than reported as comparisons.

## Certificate example

`examples/CERT_example_C4.json` is a complete `artin-CERT-1` sample certificate for the quartic cyclotomic example. It is kept as a concrete input for the independent checker and as a reference for the serialized certificate format.
