# SUPERSEDED — dimension 29 by the norm-8 frame layer

**τ(29) ≥ 209 594**, against the published 209 496 — a +98 improvement, held for part of
2026-09-15.

**Superseded by this project's own
[`../../improved/dim29-30-frame-layer/`](../../improved/dim29-30-frame-layer/), which gives
209 968.** Same template, a different layer: instead of four norm-6 heads at height √(5/2) on the
32 deep holes of D₅ — 128 points, and an axis cut from 40 to 10 — a whole Leech *frame* of
norm-8 heads at height √2 on the cross-polytope of ℝ⁵, 48 × 10 = 480 points, with 32 of the 40
axis points surviving. Better by +374. Everything below is still true; the height-√(5/2)
coincidences are real, and the layer is still the largest one that *that* height admits. It is
simply not the largest layer.

---

## The original account

# Dimension 29 — a deletion-free layer at height √(5/2) in the D₅ deep holes

**K(29) ≥ 209 594**, against the published 209 496 (Ma et al. 2025, arXiv:2511.13391, the
PackingStar configurations; the value in Cohn's table). A gain of **+98**, the largest of the three
layer gains in the 25–31 range (dimension 30 is +50, dimension 31 is +4).

    189 616 equator  +  19 840 caps  +  10 axis  +  128 layer  =  209 594

    python verify29.py     # exact, ~3 minutes; the whole claim
    python axbound29.py    # exact, seconds; the axis is at its ceiling of 10

Joint work in progress with H. Cohn and B. Lindow.

## The idea: the height is a free parameter

Every earlier layer in this repository fixed the head at `x = f/√3` (height √2 — dimensions 28,
30) or at `x = u/2` (height √3 — dimension 31), and the criterion was written as if it were a
statement about the direction system alone. It is not. A head at height `h` sits at `|x|² = 4 −
h²`, and three separate thresholds move with `h`:

| what | threshold | why |
|---|---|---|
| does the head delete equator points? | `max ⟨x,u⟩ ≤ 2` | a deleted equator point costs 1 |
| may an axis point survive? | `cos(a,w) ≤ 1/h` | the layer blocks the axis otherwise |
| how many directions per head? | `cos(w,w′) ≤ 1 − 2/h²` | two directions of one head |
| how many heads per direction? | `cos ≤ (2−h²)/(4−h²)`, so `k ≤ 2/(2−|x|²)` | the simplex bound |

For `D₅` the right height is `h² = 5/2`, where three of them coincide at once:

* `|x|² = 4 − h² = 3/2` is **exactly the norm-6 scale**, so the head is `v/2` with `v` a norm-6
  Leech vector, and `⟨v/2, u⟩ ≤ 3/2 < 2` for every minimal `u`: the layer **deletes nothing**;
* `1/h = 2/√10` is **exactly D₅'s covering cosine on S⁴**, attained at its 32 deep holes
  `(±1,±1,±1,±1,±1)/√5` — so deep holes are legal directions at all, and the axis survives;
* `1 − 2/h² = 1/5` is **exactly the cosine of two sign patterns at Hamming distance 2**, so one
  head carries **16 directions instead of 4**.

That last coincidence is where the points come from: heads are expensive (each forbids a
552-vector set of owners) and directions are free.

`|x|² = 3/2` is a knife edge in both directions. The simplex bound `k ≤ 2/(2−|x|²)` is an integer
exactly there (4 heads per direction); and `1 − 2/(4−|x|²) ≥ 1/5` makes it the smallest `|x|²` at
which Hamming-distance-2 codes are legal. Raising it to 1.6 buys a fifth head per direction but
drops each head to 4 directions, so the head count — and with it the forbidden set — rises
fivefold, and the probability that a class image is clean falls from 1.8·10⁻⁵ to 6·10⁻²⁵.

## The construction

**Heads.** Four Golay **dodecads** covering each of the 24 coordinates exactly twice, with
opposite signs on the 4 coordinates each pair shares. The resulting `v₁…v₄` have norm 48 (Cohn
units), Gram 48 on the diagonal and −16 off it, and **sum to zero** — which is precisely the
simplex bound `k ≤ 4` at `|x|² = 3/2`, attained and tight.

**Layer.** `v_i` on the 16 **even** sign patterns, `−v_i` on the 16 **odd** ones:
`4 × 32 = 128` points `(±v_i/√32, √(5/2) w_j)`. Three different families of pairs sit at exactly
2 — the same head at Hamming distance 2, two heads on one direction, and `v_i` against `−v_j` at
distance 1 — so nothing here has slack.

**Axis.** The published 40-point axis loses 32 of its points to the layer. Replacing it with
the ten points `±Eᵢ` loses none, and `axbound29.py` proves 10 is the ceiling: against all 32 deep holes
the admissible region is exactly `‖a‖₁ ≤ √2`, which together with the cap constraint
`|a_i| + |a_j| ≤ √(3/2)` forces `max|a_i| ≥ 0.922857…` — every admissible point within 22.65° of a
coordinate axis, hence at most one per `±E_i`.

**Classes.** All 14 direction groups are dirty, so every class has to avoid
`U = {u : |⟨v_i,u⟩| = 24}`, `|U| = 4320` (8 blocks of 552 overlapping in 96). A pool of 17 million
monomial images per 47 s yielded 12 185 distinct clean classes, and greedy assembly found 14
pairwise disjoint ones **on the first seed**: no owner had to be given up, so the layer is worth
its full 128.

## Why 29 and not 31

Owner density. Dimension 29 places 14 classes, `6944 = 3.53%` of the Leech shell, so a clean image
is disjoint from the other 13 with probability about `e^{−16}` and a 12 000-image pool reaches a
perfect packing at once. Dimension 31 places 42 classes at 10.60%, where the same step costs
`e^{−52}` and greedy assembly stalls at 386 shared owners. The dimension with the *best* layer
geometry has the *worst* room for it.

## Files

    verify29.py          the claim, exact in every decision (~3 min)
    verify29.log         its output
    axbound29.py         the axis ceiling of 10, exact, sympy only, no data files
    axbound29.log        its output
    data/classes.npy     14 x 496 x 24, int8: the owner classes in Cohn coordinates (norm 32)
    data/heads.npy       4 x 24, int8: the heads, norm 48, Gram 48/-16, summing to zero
    data/directions.json unit vectors of R^5 as sympy strings, in the orthonormal D5 frame of
                         Cohn's dimension-29 block: 14 direction groups (12 zero-sum triangles
                         and 2 antipodal pairs, covering the 40 cap directions once each), the
                         10 axis directions, the 32 deep holes and their sign parity
    lib/golay.py         the Golay code in the coordinates of Cohn's file (shared with dim30/31)
    lib/leech.py         the 196560 minimal vectors and the lattice membership test

The y-side geometry is written in the orthonormal `D₅` frame recovered from the published cap
directions (`research/collab2531/dim29/build29b.py`, which reads the dimension-29 block of
`dimensions25-31.txt`), in which the 40 cap directions are exactly the 40 vectors `(±1,±1,0,0,0)/√2` and the deep holes exactly `(±1,…,±1)/√5`. Everything on that side is
therefore exact, and the verifier never touches a floating-point direction. Which partition of the
40 roots into 12 triangles and 2 antipodal pairs is used does not matter — every shipped class
avoids `U` entirely, so no group is cleaner than another.

## What would improve it

Only the class. The layer is at its ceiling: 4 heads per direction is the simplex bound, 16
directions per head is `A(5,2) = 16`, so `32 × 4 = 128` is the most this height can hold; the axis
is at its exact ceiling of 10; a second layer at any other height on the same directions would
need `⟨v_i,w⟩ ≤ −3` for all four heads, which `Σv_i = 0` forbids; and free heads die for the
reason recorded in FINDINGS §49. The published cap height is forced too: `c² = 2/3` is the *smallest*
value that keeps the 496-class, hence the largest directions-per-head.

What is left is the direction-side counting slack. D₅'s triangle graph is vertex-transitive with
`n = 80` and `α = 6`, so `χ_f(G₅) = 80/6 = 40/3` exactly, and the fractional ceiling for the cap
side is `2·496·40/3 = 39680/3 = 13 226.67` against the `26·496 = 12 896` an integral partition
realises — **330.67 short**. That is the `13⅓ against 13` gap recorded in FINDINGS §§45–47, and it
is a question about the Leech class, not about dimension 29.

## History of the bound in this dimension

| date | value | source |
|---|---|---|
| 2025 | 209 496 | Ma et al., arXiv:2511.13391 (PackingStar), the value in Cohn's table: 12 zero-sum triangles + 2 antipodal pairs of cap directions, 40 axis points. The coordinates are read from `dimensions25-31.txt` of H. Cohn, *Table of kissing number bounds*, <https://hdl.handle.net/1721.1/153312> |
| 2026-09-15 | **209 594** | **this package** — the height-√(5/2) layer, 128 points, the axis replaced by `±Eᵢ` |
