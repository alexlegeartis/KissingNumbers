# SUPERSEDED — dimension 30 by the norm-8 frame layer

**τ(30) ≥ 220 494**, against the published 220 440 — a +54 improvement, held from 2026-09-13 to
2026-09-15.

**Superseded by this project's own
[`../../improved/dim29-30-frame-layer/`](../../improved/dim29-30-frame-layer/), which gives
221 012.** The same height √2, but the head is taken from a lattice vector of squared length 8
instead of 6, which raises the cap threshold from `⟨f,u⟩ ≤ 2.12` to `⟨v,u⟩ ≤ 2` *with equality*
and lets a whole Leech frame carry the whole cross-polytope of ℝ⁶: 48 × 12 = 576 points against
54. Better by +454. Everything below is still true; the E₆* deep-hole geometry is real, and the
norm-6 head is simply the more expensive one.

---

## The original account

# Dimension 30: a deletion-free layer at height √2 in the E6* directions

**τ(30) ≥ 220 494**, against the published 220 440 (Cohn's table; Ma et al. 2025,
arXiv:2511.13391, the PackingStar configurations). An improvement of **+54**: eleven *head
lines* worth six points each. The layer itself still deletes nothing; what it
gives back is 6 owners lost in the packing of the 24 class images — 12 points — not to the layer.

```
python verify30.py          # ~1 min, exact arithmetic in every decision; exits non-zero on failure
```

It reads only `data/` and rebuilds the 196 560 Leech minimal vectors from the Golay code
(`lib/golay.py`, `lib/leech.py`, the coordinates of Cohn's file), so nothing depends on the code
that produced the configuration. Every pair of points is decided exactly: integer inner products
for the ℝ²⁴ parts, sympy in ℚ(√2, √3) for the ℝ⁶ parts.

## The configuration

Norm-4 units: a point of ℝ³⁰ is `(x, y)` with `y ∈ ℝ⁶`, two points are compatible iff their inner
product is at most 2.

| | points |
|---|---|
| equator `(u, 0)`, `u` a minimal vector that is not an owner | 184 662 |
| caps `(√(2/3) u, (2/√3) z)`, `u` an owner, `z` a direction of its triangle | 3 × 11 898 |
| axis `(0, 2a)`, a rotated copy of E6 at 30° from the cap roots | 72 |
| **layer** `(±f/√3, √2 w)`, `f` one of 11 head lines, `w` one of the three directions it carries | **6 × 11** |
| **total** | **220 494** |

The first three rows are the published template: 24 classes of up to 496 owners, one per zero-sum
triangle of the cap E6. A class is *not* merely 60°-free: its owners share all three directions
of their triangle, so their inner products are at most 1 — an angle of at least arccos(1/4) =
75.52°. The classes here are
monomial images of the published 248-line class, chosen so that the layer costs nothing; 6 of
the 11 904 owner slots are empty, because two images shared that minimal vector and only one of
them may own it. Each empty slot costs two points (three caps, one equator point back), so the
layer's 66 points net **+54**.

## Why the layer is free, and why it has six points per head

A head is a Leech vector `f` of norm 6. The point `(f/√3, √2 w)` has norm 4 for any unit `w`, and

| partner | inner product | bound |
|---|---|---|
| equator `(v, 0)` | `⟨f, v⟩/√3 ≤ 3/√3 = √3` | free: a norm-6 vector meets a minimal vector at at most 3 |
| axis `(0, 2a)` | `2√2 ⟨w, a⟩ ≤ 2` iff `⟨w, a⟩ ≤ 1/√2` | the E6* deep holes meet the rotated axis copy at exactly `±1/√2`: touching, nothing deleted |
| cap of owner `u` at direction `z` | `(√2/3)⟨f, u⟩ + (2√6/3)⟨z, w⟩`, at most `2√2/3 + 1 = 1.9428` | free unless `⟨f, u⟩ = 3` **and** `⟨z, w⟩ = √6/4`, the only case above 2 |
| another layer point `(±f′/√3, √2 w′)` | `⟨f, f′⟩/3 + 2⟨w, w′⟩` | ≤ 2 for every pair used here |

so the layer deletes nothing at all; the only obstruction is an owner at `⟨f, u⟩ = 3` — 552
minimal vectors per head — sitting in a class whose triangle carries a root at cosine `√6/4` to
`w`. Each deep hole makes exactly **16 of the 24 triangles dirty**, and those 16 classes are the
ones that have to be moved off `A(f) = {u : ⟨f, u⟩ = 3}`.

Two facts multiply the yield of one such avoidance:

* **Three directions, one dirty set.** The 54 deep holes fall into **9 groups of 6** (three lines
  each) with *identical* dirty sets, and each group contains exactly two zero-sum triples of
  pairwise non-acute directions. A head vector may carry any set of directions that is pairwise
  non-acute (`|x|² + 2⟨w, w′⟩ ≤ 2` with `|x|² = 2`), so it carries a whole triple — three points
  for the same 16 dirty classes. (The largest non-acute set of deep holes has four directions, but
  a 4-set dirties 22 triangles: 4/22 is worse than 3/16.)
* **Both signs.** Every class is antipodal, so `C ∩ A(f) = ∅` is the same condition as
  `C ∩ A(−f) = ∅`: the head `−f` is free once `f` is placed. At one direction the two points of a
  line meet at `⟨f, −f⟩/3 + 2 = −2 + 2 = 0`.

So one avoidance buys **six points**.

### The head system: `⟨f, f′⟩ = 4`, not 3

The layer allows `|⟨f, f′⟩| ≤ 6 − 6M`, where `M` is the largest cosine between the two heads'
direction triples — so `M = 1/2` permits 3 and `M = 1/4` permits 4. **Each group of six deep
holes carries two zero-sum triples**, and choosing triple 1 for groups 0–5 and triple 0 for
groups 6–8 puts *all 36 group pairs* at cosine 1/4 at once, where `⟨f, f′⟩ = 4` is legal. Such
heads exist in quantity: the Leech lattice contains a 22-line clique at inner product 4.

That is worth much more than the extra room suggests, because it is not the union of the conflict
sets that matters but how their avoidances correlate. Measured over 8.2 million random monomial
images, `P(dodge both)/P(dodge a)P(dodge b)` is 1.003, 0.992, 0.987, **1.039**, **1.290** at
`⟨f, f′⟩ = 0, 1, 2, 3, 4`. The per-class supply of usable images rises accordingly —
1.04·10⁻⁵ → 7.6·10⁻⁵ at nine head lines, 6.5·10⁻⁷ → 2.1·10⁻⁶ at eleven — and eleven head lines
now cost about what nine used to. Two lines in the *same* group are the exception: they must take
opposite triples, whose cosine is 1/2, so they are capped at 3.

## What that costs, and where it stops

Each class lies in exactly 6 of the 9 dirty sets, so `N` head lines spread over the groups give
every class a *burden* of `2N/3` conflict sets to dodge, and the layer is `9b` points when `b` is
the largest burden a single class carries (the linear program `max Σ n_g` subject to
`max_class burden ≤ b` has optimum `1.5b`; the dual assigns 1/16 to every class). A random
monomial image of the class avoids one `A(f)` with probability 0.135, so a pool of depth `b`
costs roughly `0.135^-b` images per class. The integer optimum, with at most two lines per group
(a third would have to be orthogonal to one of the other two), is

| burden `b` | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|
| head lines | 9 | 10 | 11 | 12 | 14 |
| layer points | 54 | 60 | **66** | 72 | 84 |

and `b = 8` — eleven head lines, two groups doubled — is what the `⟨f, f′⟩ = 4` head system makes
affordable. Its pool here is 4.8 million images, about 40 000 usable per class.

What stops the count there is the *packing*. Twenty-four pairwise disjoint images of one 248-line
class is the hard part, and it is a maximum-clique problem in a graph of edge density 0.535 (the
chance two images are disjoint): with `N` images per class a greedy construction reaches a clique
of about `log(N)/log(1/0.535)` — 17 at `N = 50 000` — while the maximum is about twice that and we
need exactly 24. Perfect packings therefore exist in enormous numbers (`Π N_p · 0.535²⁷⁶ ≈ 10³⁷`)
but greedy cannot reach them until the pool is ~1.2 million per class; everything between is local
search. Exhaustive per-position sweeps and ruin-and-recreate get to **six** shared owners, which
is where the 12 points go.

## Files

```
verify30.py          the verification (exact); prints ALL CHECKS PASS  K(30) >= 220494
verify30.log         its output
data/heads.npy       9 x 24 int8: the head lines (Leech vectors of norm 48); both signs are used
data/classes.npy     24 x 496 x 24 int8; a row of zeros is an owner slot given up to another class
data/directions.json unit vectors of R^6 as sympy strings: the 24 zero-sum triangles of cap
                     directions, the 72 axis directions, the 27 layer directions, and
                     head_dirs -- the three directions each head line carries
lib/golay.py         the Golay code in the coordinates of Cohn's file
lib/leech.py         the 196560 minimal vectors from it
```

## The ceiling

Everything except the class and the layer is tight: the equator is the whole Leech kissing
configuration, the axis is `τ(6) = 72`, and the cap system is exactly optimal *for these
classes*: an owner is worth exactly +2 (one equator point out, three caps in) and can serve only
one triangle, so equator + caps + axis is exactly `196 560 + 2·NOWN + 72`, which at the realised
`NOWN = 24 × 496` is exactly **220 440**. The published value is therefore the exact optimum of
the template, and every point above it has to come from the layer. What is *not* proved is
`|S_z| ≤ 496`: the Delsarte bound on a class is 425 lines, not 248. The two levers are therefore

* **the class**: +96 for every extra line, against 248 realised and 425.45 by the Delsarte LP —
  the eigenvalues of the forbidden (`|ip| = 2`) graph on the 98 280 Leech lines are exactly
  4600, 1000, 76, −20, and Hoffman's ratio bound `98280·20/4620 = 9360/22` is attained by the LP.
  The same number follows with no reference to the Leech lattice, from the fact that a class is a
  two-angle line system (`|cos| ∈ {0, 1/4}`): see `research/collab2531/FRAMEWORK_2531.md`;
* **the direction set**: 72 is `τ(6)`, but the bound that actually applies is the *fractional*
  independence number of the "cosine > 1/2" graph on the directions (at most 78 by the
  dimension-6 LP). E6 ∪ E6* gives exactly 72 again, by König's theorem on its bipartite
  cosine-`√6/4` graph.

Working notes: `research/collab2531/BRAINSTORM_30.md`, FINDINGS §46, KNOWLEDGE.md §129.

## History of the bound in dimension 30

| date | bound | who |
|---|---|---|
| 2025 | 220 440 | Ma et al., arXiv:2511.13391 (PackingStar), the value in Cohn's table |
| 2026-09-13 | 220 450 | this package, five heads at two points each |
| 2026-09-14 | 220 452 | six heads |
| 2026-09-14 | 220 458 | the same six heads on a whole triple of directions (three points each) |
| 2026-09-15 | 220 490 | 9 head lines, six points each, 2 owners given up |
| 2026-09-15 | **220 494** | 11 head lines: the heads moved to `⟨f, f′⟩ = 4` — this package |
