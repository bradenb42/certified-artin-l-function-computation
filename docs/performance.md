# Cost

`examples/grid.py` runs a grid of inputs, records the time of each stage, and
regenerates every table below into `examples/grid/`.

```bash
python examples/grid.py out            # all grid entries
python examples/grid.py out d6_S6_x6-x-1
```

## Where the time goes

Seconds, on one core:

| input | n | |G| | placement | table | classes | local | conductors | Euler | root numbers | functional eq. | total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| x^3-x-1 | 3 | 6 | 0.01 | 0.01 | 0.03 | 0.03 | 0.00 | 0.01 | 0.01 | 26.8 | 27 |
| x^4-x-1 | 4 | 24 | 0.03 | 0.02 | 0.02 | 0.07 | 0.00 | 0.36 | 1.9 | 101.5 | 104 |
| x^5+x^4-4x^3-3x^2+3x+1 | 5 | 5 | 0.17 | 0.02 | 2.2 | 0.61 | 0.00 | 0.07 | 0.03 | 49.1 | 55 |
| x^5-x-1 | 5 | 120 | 0.07 | 0.15 | 0.04 | 1.3 | 0.02 | 10.5 | 16.4 | off | 31 |
| x^6-2 | 6 | 12 | 37.1 | 0.02 | 1.9 | 2.5 | 0.00 | 0.03 | wild | off | 41 |
| x^6-x-1 | 6 | 720 | 0.09 | 3.8 | 0.15 | 3.6 | 0.02 | 5.1 | 1.3 | off | 14 |
| x^7-x-1 | 7 | 5040 | 0.77 | 5.8 | 1.7 | 27.6 | 0.01 | 6.2 | off | off | 43 |

Notes on the columns:

- **Placement** is the resolvent of degree n!/|G| that locates the group in the numbering. It is negligible for the symmetric groups and 37 s for the dihedral sextic (index 60). With a descent chain in the input certificate the term disappears; the figures above are for runs without one, where degree 7 and above is reached for the symmetric groups.
- **Table** grows with |G|: 3.8 s at order 720, 5.8 s at order 5040.
- **Local** grows with the order of the decomposition group and the residue degree of the local field, which drives the precision.
- **Euler** is dominated by the identity system, and by the maximal-order computations at primes dividing the discriminants of the subgroup resolvents.
- **Functional equation** is the largest stage whenever it runs, from building the test function and assigning Frobenius classes out to the truncation point.
- **Root numbers** cost seconds at primes up to a few thousand, through Gauss sums over the residue field, and the stage is configurable for larger primes.

## Statistics over the grid

| input | |G| | classes | rational classes needing refinement | primes by route | non-unique filtrations | largest wild inertia | wild primes |
| --- | ---: | ---: | ---: | --- | ---: | --- | --- |
| x^3-x-1 | 6 | 3 | 0 | 45 cheap | 0 | none | none |
| x^4-x-1 | 24 | 5 | 0 | 46 cheap | 0 | none | none |
| x^4-2 | 8 | 5 | 0 | 45 cheap | 0 | 4 | 2 |
| C5 quintic | 5 | 5 | 1 | 8 cheap, 36 twisted, 1 direct | 0 | none | none |
| x^5-x-1 | 120 | 7 | 0 | 44 cheap | 0 | none | none |
| x^5+20x+16 | 60 | 5 | 1 | 25 cheap, 19 twisted | 0 | 5 | 2, 5 |
| x^6-2 | 12 | 6 | 0 | 44 cheap | 0 | 6 | 2, 3 |
| x^6-x-1 | 720 | 11 | 0 | 45 cheap | 0 | none | none |
| x^7-x-1 | 5040 | 15 | 0 | 46 cheap | 0 | none | none |

"Non-unique filtrations" counts primes where several filtrations are
consistent with the ramification polygons; there are none in this grid, all
of its built primes being tame. "Largest wild inertia" is the least common
multiple of the ramification indices at the wild primes.

## Harder inputs

Quartics with inertia of order 3 and up to three ramified primes take one to four seconds. Sextics with Galois groups of order 18 to 72 take two to thirty-five seconds. The septic with group PSL(2,7) takes about a minute: the placement resolvent has index 30, and the pair of classes of order 7 is resolved by the direct route because the twisted resolvents would need a Galois ring of degree 42 at several thousand bits.
