# The certificate and the checker

Every run writes `CERT.json`: one self-contained file holding the group and
its classes, the character table with its certificate, the ramified-prime
witnesses, the precision log, the Frobenius classes with the resolvents that
determined them, the local groups and matchings, the filtrations and conductor
exponents, the Euler factors, the archimedean data, the root numbers, the
Dirichlet coefficients and the functional-equation defects.

`artinverify` reads that file and rejects it at the first claim that fails.
The package imports only the permutation group module from the main package.
Its cyclotomic arithmetic, its polynomial arithmetic over finite fields, its
Sturm sequences, its discriminants and its maximal-order checks are separate
implementations, so the two agree only when both are right.

```bash
python -m artinverify.verify runs/example/CERT.json
```

The claims are checked in this order.

1. **Group and table.** The recorded class representatives form a system with the recorded sizes and orders; the power maps are right; the central characters match the table; the eigenvector identities hold for the recorded classes, with the class matrices recomputed; the joint eigenvalue tuples are distinct; the rows are orthonormal.
2. **Ramified primes.** The recorded factors multiply to f; the discriminants and their prime factorisations are recomputed; each prime's witness is rechecked: an odd valuation, the Dedekind criterion re-run from the recorded factorisation, or a maximal-order witness (containment, closure under multiplication, index, radical, maximality, decomposition).
3. **Precision.** The bound is recomputed from the recorded pairs and every logged consultation is checked against the pairs registered at that time.
4. **Frobenius classes.** At each prime the cycle type of f mod ell must be the cycle type of the recorded class, and the factorisation type of every non-excluded resolvent must be that class's signature on the cosets.
5. **Local data.** The decomposition group lies in G, has order e f, contains the Frobenius and normalises the inertia group; the two matchings differ by an element of G; the orbit sizes match the factor data.
6. **Filtrations and conductors.** The recorded ramification numbers are transported by the matching, the wild inertia subgroup is identified, and every conductor exponent is recomputed from those numbers and the table, together with both identities.
7. **Euler factors.** Each is recomputed from the local groups and the table; the degree, the conductor and Swan exponents and the product identity are checked.
8. **Archimedean data.** The real root count, the class of complex conjugation, the Gamma exponents and the parity are recomputed.
9. **Root numbers.** Moduli, the identification of the local factors as roots of unity, the product, and the predictions for self-dual characters.
10. **Functional equation.** The recorded coefficients are checked at every prime against the recorded classes and Euler factors, and the defects are recomputed with the checker's own test function and its own tail bound.

## Testing the checker

`artinverify.mutate` produces copies of an accepted certificate in which one
datum has been altered: one ramification number, one conductor exponent, one
Frobenius class, one Euler-factor coefficient, one local root number. Each is
rejected by the claim that governs it. The test suite runs this end to end.

## What acceptance means

The checker confirms that the recorded data are consistent and correctly derived from f and the given group. That the given group is the Galois group of f is established by the descent step that produces the input certificate. As a cross-check on it, the driver verifies that the Frobenius classes observed at the small primes generate G and warns when they do not, since a group that is too large passes every later stage unnoticed.
