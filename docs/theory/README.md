# Derivations

These are the arguments the code implements: what each stage computes, why the
quantity it computes is the right one, and what a checker has to verify.  They
were written before the code and are the specification it was built against;
where an argument fixes a normalisation or a bound, the corresponding module
follows it exactly.

Notation is uniform across the pages: `f` is a monic separable polynomial in
`Z[x]`, `G` its Galois group as a permutation group on the roots, `N` the
splitting field, `chi` an irreducible character of `G`, and `ell` a prime.

## Setting up

| derivation subject code                                                                       |                                                                                                                                                    |                                        |
| --------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| [Character table, Schur indices and matrix models](https://claude.ai/chat/character-table.md) | the table with its certificate, Schur index bounds, models over cyclotomic fields, and why none of the L-function data depends on the model chosen | `artin/chartable.py`, `artin/schur.py` |
| [The ramified primes](https://claude.ai/chat/ramified-primes.md)                              | which primes ramify, decided from the discriminants of the factors, with witnesses                                                                 | `artin/ramified.py`                    |
| [The precision policy](https://claude.ai/chat/precision-policy.md)                            | one bound covering every test, and the sharper bound test by test                                                                                  | `artin/precision.py`                   |

## Frobenius classes at the unramified primes

| derivation subject code                                                                               |                                                                                         |                    |
| ----------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | ------------------ |
| [Rational classes and what coset actions determine](https://claude.ai/chat/rational-classes.md)       | coset actions determine the rational class of the Frobenius and nothing finer           | `artin/classes.py` |
| [Choosing and realising the separating coset actions](https://claude.ai/chat/separating-subgroups.md) | the set-cover choice of separating subgroups, their resolvents, and the excluded primes | `artin/classes.py` |
| [Refinement inside a rational class](https://claude.ai/chat/cyclotomic-refinement.md)                 | twisted resolvents over `Z[zeta_o]` and what they resolve                               | `artin/twisted.py` |
| [The direct route to a Frobenius class](https://claude.ai/chat/direct-route.md)                       | residue-field arithmetic when the cheaper tests are expensive                           | `artin/direct.py`  |

## The ramified primes

| derivation subject code                                                                                  |                                                                                 |                       |
| -------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- | --------------------- |
| [Local Galois groups and compositum assembly](https://claude.ai/chat/local-galois-groups.md)             | the decomposition group at a ramified prime                                     | `artin/local.py`      |
| [Matching the local roots to the global numbering](https://claude.ai/chat/matching.md)                   | the matching problem and its solution                                           | `artin/matching.py`   |
| [Independence of the local data from the choices](https://claude.ai/chat/well-definedness.md)            | the matchings form one coset, so the recorded groups are well defined           | `artin/matching.py`   |
| [The ramification filtration](https://claude.ai/chat/ramification-filtration.md)                         | lower and upper numbering, Herbrand, the polygons                               | `artin/filtration.py` |
| [When the polygons do not determine the filtration](https://claude.ai/chat/filtration-uniqueness.md)     | the enumeration of consistent filtrations and the resolvents that separate them | `artin/filtration.py` |
| [Artin conductor exponents](https://claude.ai/chat/conductor-exponents.md)                               | the exponents and the identity fixing their weighted sum                        | `artin/filtration.py` |
| [Global conductors, determinant characters and cross-checks](https://claude.ai/chat/global-conductor.md) | the determinant character and the conductor-discriminant identities             | `artin/euler.py`      |

## L-function data

| derivation subject code                                                                                  |                                                                 |                               |
| -------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- | ----------------------------- |
| [Euler factors at a ramified prime](https://claude.ai/chat/euler-factors.md)                             | three routes to the same factor and the dimension checks        | `artin/euler.py`              |
| [Identities among the Euler factors](https://claude.ai/chat/euler-identities.md)                         | what the subgroup identities pin down, and what they leave free | `artin/euler.py`              |
| [Complex conjugation, Gamma factors and archimedean root numbers](https://claude.ai/chat/archimedean.md) | the class of complex conjugation from real root counts          | `artin/archimedean.py`        |
| [Local root numbers at tame primes](https://claude.ai/chat/tame-root-numbers.md)                         | Gauss sums and the Lubin-Tate normalisation                     | `artin/rootnumber.py`         |
| [Local root numbers at wild primes](https://claude.ai/chat/wild-root-numbers.md)                         | the Brauer induction argument                                   | derivation only               |
| [The global root number](https://claude.ai/chat/global-root-numbers.md)                                  | the product formula and the self-dual shortcuts                 | `artin/stages_rootnumbers.py` |

## Testing the result

| derivation subject code                                                                |                                                                      |                                                                                       |
| -------------------------------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| [The smoothed functional-equation test](https://claude.ai/chat/functional-equation.md) | the test function, the tail bound, and which errors the test detects | `artin/analytic.py`, `artin/falsifier.py`                                             |
| [Dedekind zeta identities for the subfields](https://claude.ai/chat/subfield-zeta.md)  | the subfield identities and what they add                            | `artin/stages_analytic.py`                                                            |
| [The certificate](https://claude.ai/chat/certificate.md)                               | what the run records                                                 | `artin/certwriter.py`                                                                 |
| [The verifier](https://claude.ai/chat/verifier.md)                                     | the claims a checker walks, and what acceptance proves               | `artinverify/`                                                                        |
| [External correctness checks](https://claude.ai/chat/external-checks.md)               | the three families checked against independent computations          | `artin/check_abelian.py`, `artin/check_splitting_field.py`, `artin/check_newforms.py` |
| [Cost](https://claude.ai/chat/cost.md)                                                 | how the cost decomposes, and what is feasible at a given budget      | `examples/grid.py`                                                                    |

## Reading them

The pages are self-contained but cross-referenced.  A reader following the computation from the start can go in the order above; a reader interested in one stage can start at its page and follow the links backwards.  Each numbered result is stated and proved in place, and the code comments name the result they implement wherever a choice would otherwise look arbitrary.
