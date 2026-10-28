# Continuous tiny-point completion in 43 dimensions

The 28 October 2026 certificate proves **tau(43) >= 2,919,782**, improving
this repository's preceding **2,919,743** by **39 points**. It deletes 46
points from the previous 3,408-point tiny layer and inserts 85 rational-head
points, leaving **3,447** explicit tiny points. The 118 complete families
and five pure heads are the same inputs as before.

The count is

```text
1,084,590 + 1,388,940 + 413,100 + 29,700 + 3,447 + 5 = 2,919,782.
```

The complete-family counts retain the published P48n minimum/design inputs
and population identities in `../POPULATIONS.txt`. The explicit tiny-body
pool, new coordinates and every required angle are checked directly. This
update does not freshly enumerate the large body shells or claim optimality.

## Replay

From this directory:

```bash
python reproduce.py
python ../partial43/verify.py
python verify_coordinates.py
python negative_controls.py
```

All four require only Python 3.9+ and NumPy. The root package's `../verify.py`
also runs both mathematical representations. Both check **6,370,783**
angle conditions, including all **5,939,181** distinct tiny-point pairs.

The ambient replay uses integer products of 48-coordinate body vectors and
squares positive expressions. The independent replay reconstructs coefficient
coordinates in a positive-definite rank-40 metric, reconstructs the linear
residue map from eight family residues, and compares integer products with
integer square-root cutoffs. Both check the same **22,001** short-difference
refinements. Neither accepts a floating-point angle or an optimizer score.

`negative_controls.py` checks that both implementations reject duplicate
points and a deliberately invalid angle to a pure head. Its receipts and the
two successful mathematical replays are supplied below.

## Why the continuous search helped

For a tiny body, write its normalized point as `(a L/sqrt(3), h)` divided by
`sqrt(10 + |h|^2)`, with `L L^T = q`. Introduce the unit four-vector

```text
t = (sqrt(10), h) / sqrt(10 + |h|^2).
```

If a retained family has body-radius index `r`, head `g`, and a bound `beta`
on twice the body inner product, its angle constraint becomes linear in `t`:

```text
(beta / sqrt(10*(r+|g|^2)), 2*g / sqrt(r+|g|^2)) . t <= 1.
```

Retained individual tiny points give constraints of the same form, using
their actual body products. Pure heads give linear constraints on the last
three entries. Together with `0 <= t[0] <= 1`, these define a polytope.
A vertex `v` with `|v| > 1` and `v[0] > 0` gives a feasible unit direction
by radial normalization. Its head is `sqrt(10)*v[1:]/v[0]`.

The search removes small groups of old constraints and explores this
continuous polytope. It rounds proposed heads to rationals, then checks them
with arbitrary-size integers. Joint packing of the proposals is essential:
individually admissible points can conflict with one another. Three rounds
increased the tiny count by 32, 5 and 2; a fourth sampled round found no gain.
The final selection from 594 rational proposals deletes 46 old points and
inserts 85. These finite searches establish no construction ceiling.

This transfers PackingStar's fill/correct principle to the existing body/head
model. The new step is continuous head placement after shared deletions;
no reinforcement-learning implementation or new Gram-matrix rank assumption
is used.

## Reproduction and optional further search

`reproduce.py` assembles the exact registered arrays from `previous_points.npz`
and `delta.npz`. It compares them with `../partial43/explicit_points_43.npz`;
`--output /tmp/reproduced43.npz` also writes a separate copy. This reproduction
uses no floating point or optimizer.

The discovery tools additionally require SciPy:

```bash
python search.py --body 14496 --delete 152 1348 1690
python select.py --output /tmp/proposed43-delta.npz
```

`search.py` explores continuous head placements for one body after specified
old-point deletions. Every printed head passes exact checks against the
retained configuration; multiple printed heads still need mutual checks.
`select.py` rebuilds all old/new and new/new conflicts of the shipped finite
proposal pool and optimizes the net count. On the recorded environment it
recovers the registered 85-insertion, 46-deletion delta. Solver output alone
is not a proof. Numerical polytope discovery can vary across software versions;
the supplied rational certificate and its exact assembly are portable.

## Files

- `previous_points.npz`: retained 3,408-point predecessor, with its original denominator.
- `delta.npz`: deleted indices and all 85 inserted rational-head coordinates.
- `candidate_pool.npz`: 594 rational proposals from the continuous searches.
- `reproduce.py`: exact coordinate assembly and comparison.
- `verify_coordinates.py`: independent coefficient-coordinate proof.
- `search.py`, `select.py`: optional discovery tools.
- `negative_controls.py`: invalid-certificate rejection checks.
- `ambient_verification.json`, `coordinate_verification.json`: full exact replay receipts.
- `reproduction.json`, `negative_controls.json`: assembly and rejection receipts.
- `search_receipt.json`: finite search scope and fingerprints.
