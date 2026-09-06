# Scope

What the package computes.

## Covered

- Any monic separable f in Z[x] whose Galois group is supplied as a permutation group on the roots.
- Every irreducible character: conductor, Euler factors at all primes up to a bound and at the ramified primes, Gamma factor, and, at tamely ramified primes, the root number.
- Certificates for all of it, checkable by a separate program.

## Wild ramification

The local Galois group is constructed at primes whose ramification indices are prime to the prime itself. At a wildly ramified prime the run records the prime and the ramification data known from the factors, and marks the conductor exponent, the Euler factor and the root number there as not computed; the characters involved carry that mark into their global conductor and root number. Data at the remaining primes are complete.
Reference values for several wild primes (the ramification index, residue degree, number of primes and discriminant valuation) are recorded in `examples/results_splitting_field.json`, read off explicitly constructed Galois closures.

## Degree and group size

When the input certificate contains a descent chain, the group is placed in the root numbering along that chain and the cost is small. Otherwise the placement uses a resolvent of degree n!/|G|, which is cheap for the symmetric groups and grows fast otherwise; runs of that kind reach degree 6 for general groups and degree 7 to 8 for the symmetric ones.

## Root numbers

Root numbers are computed at tamely ramified primes. The Gauss sums run over the residue field, so the cost grows with the ramified prime; the stage is configurable and is normally left on up to a few thousand, above which the root numbers are marked as not computed. A character ramified wildly is marked the same way.

## Truncation in the analytic test

The functional-equation test compares truncated sums against an explicit tail bound, so it resolves errors above that bound at primes below the truncation point. Each verdict records whether the bound was tight enough to be informative, and the perturbation experiment in `examples/` shows which errors are caught at a given truncation point.
