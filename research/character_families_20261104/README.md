# Coupled character-family discovery campaign

No new kissing-number record was found. The substantive engine change is
joint assignment/assembly of complete sign-code families using exact Walsh
kernels, including sparse stabilizers and removable base families. Continuous
proposals can deform every supplied orbit shape together. See the
[operator documentation](../../discovery/docs/CHARACTER_FAMILIES.md).

| Experiment | Actual search | Strongest result |
|---|---|---|
|17D mixed amplitudes|1200 new integer norm72 sign orbits generated;192 enter the joint model together with all406 incumbent families|5730; gain0, all old families removable; numerical solver optimum only|
|19D extra heads|All112 cap shapes can change jointly; one or two additional missing-head orbits proposed, six numerical restarts|All proposals violate the retained base; worst cosine about0.53033|
|21D nonlinear labels|12 joint radius-two/larger neighborhoods,232–254 cosets per solve; exact Walsh table independently agrees with the line rule|960 of1024 per component,3840 additions in four components; gain0|
|21D continuous shapes|16/32-coset full-code proposals and a61-of64-coset proposal; every selected amplitude shape changes jointly|Best61-coset proposal has max cosine0.51423868, above0.5; not feasible|
|21D stochastic labels|8 restarts of3m updates, all labels and holes change|No gain over960 per component|
|21D hole proposals|16 restarts of1m swaps with63 holes|At least13 uncovered vertices in the best proposal; no admissible labelling or configuration|

The17D candidates have more than eight nonzero head coordinates and at least
two distinct nonzero magnitudes. Their heights are sampled from0,±1,±2,±3.
This avoids repeating the parent's full flat-coset/norm-six tier catalogue.
Whole RM(1,4) orbits, including all roots, octads, caps and poles of the5730
control, remain removable. No point antipodality assumption is added by the
optimizer, though complete chosen sign orbits enforce their actual symmetry.
Half the192 new families are cheapest by exact incumbent-blocker price;
half are sampled from the remaining generated families. The full generated
price census is saved. Even the cheapest generated32-point orbit blocks856
old points. This diagnoses the weakness of this amplitude grammar; it does
not rule out coordinated trades or other shapes.

## Current19D seed with rational coordinates

The primary source is Kevin Li,
[Improved Kissing Numbers in Nineteen and Twenty-One Dimensions](https://hexagonmath.org/2610.00248v1),
14 October2026. `build19.py` reconstructs the130 octads from the printed binary
generators, the10668-point hyperplane body, the112 printed coset representatives,
their1792 selected words, and the forced marks. It uses diagonal tail factor
817/1000 and axis mark747/1000, replacing the source'ssqrt(2/3) and373/500.
These are altered coordinates at the existing count, not a new lower bound.

The existing unequal-norm exact backend checks all77619570 pairs and proves
rank19. All raw coordinates are integers; the unused radical coefficients
are zero. A separate exact block calculation obtains maxima for all five
norm-class pair types involving added points. Python-integer comparisons
then prove cosine≤14999/30000, i.e. a gap of1/30000 below1/2, for every such
pair. The body retains its original contacts. This gives a less fragile seed
for future coupled searches. `control19_source.json` records the exact maxima.

The19D numerical searches used mark746/1000, explicitly fixed in their script,
and tail817/1000. Their zero-extra-family control has loss0; adding one or two
families collapses toward a positive-loss stationary point. This is numerical
failure, not an impossibility certificate. The best shapes and full constraint
counts are saved. No infeasible point collection is promoted to exact status.

The related GitHub repository `Felpix-Studios/D19-D21-Kissing`, read at commit
5cec83161d2ed017ad4f43e359ea691f3f39ebcf, contains older12276/30779 constructions.
It was not mistaken for the current12460/31560 recipe. The current formula
source above supplied the choices used here.

## Scope and reproduction

The21D line-label control uses the primary source's punctured Golay polynomial
and parity map. Its3840 additions are a reconstruction. Neither a replay of
that count nor an exact kernel is credited as held-out discovery. The earlier
1–24 benchmark remains a distinct, unchanged measurement. A61-of64-coset
continuous proposal would add3904 points if feasible; its numerical failure
does not establish any bound on arbitrary nonlinear assignments.

```sh
OPENBLAS_NUM_THREADS=1 python3 research/character_families_20261104/coupled17.py
OPENBLAS_NUM_THREADS=1 python3 research/character_families_20261104/build19.py
OPENBLAS_NUM_THREADS=1 python3 research/character_families_20261104/continuous19.py
OPENBLAS_NUM_THREADS=1 python3 research/character_families_20261104/discrete21.py
OPENBLAS_NUM_THREADS=1 python3 research/character_families_20261104/continuous21.py --extra 2 --fill-one --restarts 1 --iterations 1000
OPENBLAS_NUM_THREADS=1 python3 -m unittest discovery.tests.test_characters discovery.tests.test_numerical_shapes
```

Stochastic C++ proposal commands and exact input hashes are recorded in
`manifest.json`. The first three61-coset numerical runs completed optimization
but hit an output-path error; their diagnostics are preserved in a log. The
corrected script reran the first restart and saved its complete best shapes.
No result from that failed export is used as a geometric certificate.

The mathematical blocker remains substantive: no certified1155+,5731+,
12461+ or31561+ construction has been found. The negative outcomes concern
specified finite catalogues or numerical searches, not unrestricted dimensions.
