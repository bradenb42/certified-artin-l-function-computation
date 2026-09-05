# Performance

Recorded cost decomposition and validation statistics for the current benchmark grid.

Raw data: [`examples/results/performance_grid.json`](../examples/results/performance_grid.json).

## Cost decomposition over the grid (seconds per stage)

| input | n | |G| | log10 disc | A3 numbering and group placement | A1 table schur models | A2 ramified primes | A3 precision check | A4 A7 class assignment | A4 A7 class confirmation | A8 A10 local descents matching | A11 A13 filtration conductors | A14 A16 euler identities | A17 archimedean | A18 A20 root numbers | A21 A22 functional equation | total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| d3_S3_x3-x-1 | 3 | 6 | 1.4 | 0.007 | 0.009 | 0.005 | 0.0 | 0.029 | 0.165 | 0.026 | 0.001 | 0.011 | 0.001 | 0.014 | 26.782 | 27.1 |
| d4_S4_x4-x-1 | 4 | 24 | 2.5 | 0.025 | 0.02 | 0.002 | 0.0 | 0.021 | 0.331 | 0.07 | 0.001 | 0.361 | 0.003 | 1.876 | 101.544 | 104.3 |
| d4_D4_x4-2 | 4 | 8 | 3.3 | 0.011 | 0.01 | 0.001 | 0.0 | 0.023 | 0.181 | 0.026 | 0.0 | 0.02 | 0.002 | 0.0 | 7.428 | 7.7 |
| d5_C5_cond11 | 5 | 5 | 4.2 | 0.166 | 0.017 | 0.001 | 0.0 | 2.193 | 2.565 | 0.614 | 0.004 | 0.07 | 0.006 | 0.029 | 49.076 | 54.8 |
| d5_S5_x5-x-1 | 5 | 120 | 3.5 | 0.069 | 0.154 | 0.002 | 0.0 | 0.044 | 2.51 | 1.327 | 0.017 | 10.471 | 0.014 | 16.416 | 0.001 | 31.0 |
| d5_A5_x5+20x+16 | 5 | 60 | 9.0 | 0.021 | 0.06 | 0.014 | 0.0 | 4.139 | 4.326 | 0.812 | 0.001 | 8.503 | 0.005 | 0.0 | 59.565 | 77.5 |
| d6_S6_x6-x-1 | 6 | 720 | 4.7 | 0.094 | 3.834 | 0.003 | 0.0 | 0.153 | 0.0 | 3.599 | 0.019 | 5.13 | 0.036 | 1.283 | 0.0 | 14.2 |
| d6_D6_x6-2 | 6 | 12 | 6.2 | 37.067 | 0.017 | 0.002 | 0.0 | 1.868 | 0.0 | 2.45 | 0.004 | 0.027 | 0.004 | 0.0 | 0.0 | 41.4 |
| d7_S7_x7-x-1 | 7 | 5040 | 5.9 | 0.872 | 5.807 | 0.003 | 0.0 | 1.705 | 0.0 | 27.569 | 0.016 | 6.183 | 0.609 | 0.0 | 0.0 | 42.8 |

## Statistics

| input | |G| | r | rational classes needing refinement | primes by method | undetermined filtrations | largest wild inertia (lcm e) | wild primes | largest Brauer datum |
|---|---|---|---|---|---|---|---|---|
| d3_S3_x3-x-1 | 6 | 3 | 0 | {'A5': 45, 'A6': 0, 'A7': 0} | 0 | 0 | [] | n/a (A19 not built) |
| d4_S4_x4-x-1 | 24 | 5 | 0 | {'A5': 46, 'A6': 0, 'A7': 0} | 0 | 0 | [] | n/a (A19 not built) |
| d4_D4_x4-2 | 8 | 5 | 0 | {'A5': 45, 'A6': 0, 'A7': 0} | 0 | 4 | [2] | n/a (A19 not built) |
| d5_C5_cond11 | 5 | 5 | 1 | {'A5': 8, 'A6': 36, 'A7': 1} | 0 | 0 | [] | n/a (A19 not built) |
| d5_S5_x5-x-1 | 120 | 7 | 0 | {'A5': 44, 'A6': 0, 'A7': 0} | 0 | 0 | [] | n/a (A19 not built) |
| d5_A5_x5+20x+16 | 60 | 5 | 1 | {'A5': 25, 'A6': 19, 'A7': 0} | 0 | 5 | [2, 5] | n/a (A19 not built) |
| d6_S6_x6-x-1 | 720 | 11 | 0 | {'A5': 45, 'A6': 0, 'A7': 0} | 0 | 0 | [] | n/a (A19 not built) |
| d6_D6_x6-2 | 12 | 6 | 0 | {'A5': 44, 'A6': 0, 'A7': 0} | 0 | 6 | [2, 3] | n/a (A19 not built) |
| d7_S7_x7-x-1 | 5040 | 15 | 0 | {'A5': 46, 'A6': 0, 'A7': 0} | 0 | 0 | [] | n/a (A19 not built) |

## Conductor tables

Exponents at the tame ramified primes; `unknown_at` lists wild primes.

### d3_S3_x3-x-1

- chi_1: {'23': 0} partial conductor 1 unknown at []
- chi_2: {'23': 1} partial conductor 23 unknown at []
- chi_3: {'23': 1} partial conductor 23 unknown at []

### d4_S4_x4-x-1

- chi_1: {'283': 0} partial conductor 1 unknown at []
- chi_2: {'283': 1} partial conductor 283 unknown at []
- chi_3: {'283': 1} partial conductor 283 unknown at []
- chi_4: {'283': 2} partial conductor 80089 unknown at []
- chi_5: {'283': 1} partial conductor 283 unknown at []

### d4_D4_x4-2

- chi_1: {} partial conductor 1 unknown at ['2']
- chi_2: {} partial conductor 1 unknown at ['2']
- chi_3: {} partial conductor 1 unknown at ['2']
- chi_4: {} partial conductor 1 unknown at ['2']
- chi_5: {} partial conductor 1 unknown at ['2']

### d5_C5_cond11

- chi_1: {'11': 0} partial conductor 1 unknown at []
- chi_2: {'11': 1} partial conductor 11 unknown at []
- chi_3: {'11': 1} partial conductor 11 unknown at []
- chi_4: {'11': 1} partial conductor 11 unknown at []
- chi_5: {'11': 1} partial conductor 11 unknown at []

### d5_S5_x5-x-1

- chi_1: {'19': 0, '151': 0} partial conductor 1 unknown at []
- chi_2: {'19': 1, '151': 1} partial conductor 2869 unknown at []
- chi_3: {'19': 3, '151': 3} partial conductor 23615200909 unknown at []
- chi_4: {'19': 1, '151': 1} partial conductor 2869 unknown at []
- chi_5: {'19': 3, '151': 3} partial conductor 23615200909 unknown at []
- chi_6: {'19': 2, '151': 2} partial conductor 8231161 unknown at []
- chi_7: {'19': 3, '151': 3} partial conductor 23615200909 unknown at []

### d5_A5_x5+20x+16

- chi_1: {} partial conductor 1 unknown at ['2', '5']
- chi_2: {} partial conductor 1 unknown at ['2', '5']
- chi_3: {} partial conductor 1 unknown at ['2', '5']
- chi_4: {} partial conductor 1 unknown at ['2', '5']
- chi_5: {} partial conductor 1 unknown at ['2', '5']

### d6_S6_x6-x-1

- chi_1: {'67': 0, '743': 0} partial conductor 1 unknown at []
- chi_2: {'67': 1, '743': 1} partial conductor 49781 unknown at []
- chi_3: {'67': 3, '743': 3} partial conductor 123364683646541 unknown at []
- chi_4: {'67': 4, '743': 4} partial conductor 6141217316608457521 unknown at []
- chi_5: {'67': 2, '743': 2} partial conductor 2478147961 unknown at []
- chi_6: {'67': 1, '743': 1} partial conductor 49781 unknown at []
- chi_7: {'67': 6, '743': 6} partial conductor 15218845171211140441021264681 unknown at []
- chi_8: {'67': 3, '743': 3} partial conductor 123364683646541 unknown at []
- chi_9: {'67': 6, '743': 6} partial conductor 15218845171211140441021264681 unknown at []
- chi_10: {'67': 4, '743': 4} partial conductor 6141217316608457521 unknown at []
- chi_11: {'67': 8, '743': 8} partial conductor 37714550129811583584401487826861465441 unknown at []

### d6_D6_x6-2

- chi_1: {} partial conductor 1 unknown at ['2', '3']
- chi_2: {} partial conductor 1 unknown at ['2', '3']
- chi_3: {} partial conductor 1 unknown at ['2', '3']
- chi_4: {} partial conductor 1 unknown at ['2', '3']
- chi_5: {} partial conductor 1 unknown at ['2', '3']
- chi_6: {} partial conductor 1 unknown at ['2', '3']

### d7_S7_x7-x-1

- chi_1: {'776887': 0} partial conductor 1 unknown at []
- chi_2: {'776887': 1} partial conductor 776887 unknown at []
- chi_3: {'776887': 5} partial conductor 283001847899423500561443133207 unknown at []
- chi_4: {'776887': 1} partial conductor 776887 unknown at []
- chi_5: {'776887': 9} partial conductor 103090984807943025944764838037316061168780943464484727 unknown at []
- chi_6: {'776887': 10} partial conductor 80090045914488433597150520728296362813230720825293146104849 unknown at []
- chi_7: {'776887': 5} partial conductor 283001847899423500561443133207 unknown at []
- chi_8: {'776887': 4} partial conductor 364276719650893245171361 unknown at []
- chi_9: {'776887': 10} partial conductor 80090045914488433597150520728296362813230720825293146104849 unknown at []
- chi_10: {'776887': 5} partial conductor 283001847899423500561443133207 unknown at []
- chi_11: {'776887': 10} partial conductor 80090045914488433597150520728296362813230720825293146104849 unknown at []
- chi_12: {'776887': 11} partial conductor 62220915500369175711989476597043976416882375009799516397957825063 unknown at []
- chi_13: {'776887': 10} partial conductor 80090045914488433597150520728296362813230720825293146104849 unknown at []
- chi_14: {'776887': 20} partial conductor 6414415454584865433839689336775354230592213635977501595897821879115089669868586820377874035607872530199638140901312801 unknown at []
- chi_15: {'776887': 15} partial conductor 22665630992149900224666345275509242975139027694331455377443469506077343354868850295620743 unknown at []

## Verification defects (A21, relative to the bound)

### d3_S3_x3-x-1

- chi_2: t=1.1: defect 2.89e-15 bound 7.69e-119, t=1.3: defect 8.88e-16 bound 1.28e-111
- chi_3: t=1.1: defect 2.22e-16 bound 2.23e-80, t=1.3: defect 6.11e-16 bound 3.72e-73

### d4_S4_x4-x-1

- chi_2: t=1.1: defect 8.88e-16 bound 3.09e-163, t=1.3: defect 2.11e-15 bound 5.15e-156
- chi_3: t=1.1: defect 2.66e-15 bound 1.59e-124, t=1.3: defect 1.33e-15 bound 2.65e-117
- chi_4: t=1.1: defect 7.46e-14 bound 1.99e-10, t=1.3: defect 2.68e-12 bound 5.02e-09
- chi_5: t=1.1: defect 1.11e-15 bound 5.25e-86, t=1.3: defect 2.33e-15 bound 8.75e-79

### d5_C5_cond11

- chi_2: t=1.1: defect 7.70e-13 bound 1.85e-135, t=1.3: defect 6.60e-13 bound 3.08e-128
- chi_3: t=1.1: defect 3.24e-14 bound 1.85e-135, t=1.3: defect 2.82e-14 bound 3.08e-128
- chi_4: t=1.1: defect 3.26e-14 bound 1.85e-135, t=1.3: defect 2.82e-14 bound 3.08e-128
- chi_5: t=1.1: defect 7.69e-13 bound 1.85e-135, t=1.3: defect 6.60e-13 bound 3.08e-128
