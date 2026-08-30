# Dimensions 73–95: the cap construction over Γ₇₂, with zero-sum triples

Twenty-three dimensions, gains from **+4 214** to **+184 807 648**.

```
python scripts/calibrate.py      # ~3 min  THE KEY CHECK: reproduces Cohn's table, dims 26-31
python scripts/triples.py        # instant: zero-sum triple partitions of A_1 ... E_8
python final73-80.py             # dimensions 73-80, from the 120 shipped classes
python final81-95.py             # dimensions 81-95, from the 32000-class size table
python scripts/verify_classes.py data/disjoint_classes.npz     # exact check of the 120
```

| dim | k | W | scheme | gain | **this work** | previously |
|---|---|---|---|---|---|---|
| 73 | 1 | A₁ | 1 pair + 2 poles, t=3/4 | +4 214 | **6 218 179 814** | 6 218 175 602 |
| 74 | 2 | A₂ | 2 triples + 6 poles, t=2/3 | +16 854 | **6 218 192 454** | 6 218 175 606 |
| 75 | 3 | A₃ | 4 triples + 12 poles, t=2/3 | +33 708 | **6 218 209 308** | 6 218 175 612 |
| 76 | 4 | D₄ | 8 triples + 24 poles, t=2/3 | +67 416 | **6 218 243 016** | 6 218 175 624 |
| 77 | 5 | D₅ | 12 triples + 2 pairs + 40 poles, t=2/3 | +109 552 | **6 218 285 152** | 6 218 175 640 |
| 78 | 6 | E₆ | 24 triples + 72 poles, t=2/3 | +202 240 | **6 218 377 840** | 6 218 175 672 |
| 79 | 7 | E₇ | 42 triples + 126 poles, t=2/3 | +353 854 | **6 218 529 454** | 6 218 175 726 |
| 80 | 8 | E₈ | 80 triples + 240 poles, t=2/3 | **+673 628** | **6 218 849 228** | 6 218 175 840 |
| 81 | 9 | record | 102 T +  0 P +306 poles | +858 670 | 6 219 034 270 | 6 218 175 906 |
| 82 | 10 | record | 170 T +  0 P +510 poles | +1 429 978 | 6 219 605 578 | 6 218 176 110 |
| 83 | 11 | Z[sqrt2] | 84 T + 176 P +604 poles | +1 446 090 | 6 219 621 690 | 6 218 176 204 |
| 84 | 12 | K_k | 252 T +  0 P +756 poles | +2 117 704 | 6 220 293 304 | 6 218 176 441 |
| 85 | 13 | K_k | 306 T +  0 P +918 poles | +2 569 814 | 6 220 745 414 | 6 218 176 754 |
| 86 | 14 | record | 644 T +  0 P +1932 poles | +5 392 244 | 6 223 567 844 | 6 218 177 532 |
| 87 | 15 | Lambda_k | 780 T +  0 P +2562 poles | +6 527 998 | 6 224 703 598 | 6 218 178 164 |
| 88 | 16 | Lambda_k | 1440 T +  0 P +4320 poles | +12 030 892 | 6 230 206 492 | 6 218 179 920 |
| 89 | 17 | Lambda_k | 1782 T +  0 P +5704 poles | +14 879 248 | 6 233 054 848 | 6 218 181 330 |
| 90 | 18 | Lambda_k | 2466 T +  0 P +7624 poles | +20 569 188 | 6 238 744 788 | 6 218 183 254 |
| 91 | 19 | Lambda_k | 3556 T +  0 P +11934 poles | +29 625 434 | 6 247 801 034 | 6 218 187 548 |
| 92 | 20 | Lambda_k | 5800 T +  0 P +19415 poles | +48 230 067 | 6 266 405 667 | 6 218 195 048 |
| 93 | 21 | Lambda_k | 9240 T +  0 P +29712 poles | +76 674 760 | 6 294 850 360 | 6 218 205 368 |
| 94 | 22 | Lambda_k | 16632 T +  0 P +49852 poles | +137 564 000 | 6 355 739 600 | 6 218 225 496 |
| **95** | 23 | Lambda_k | 31050 T +  0 P +93074 poles | +254 889 446 | **6 473 065 046** | 6 218 268 750 |

Every dimension 81-95 now has a construction of its own: dimension 81 clears dimension 80 by
**184 904**, so no row here is carried by monotonicity from the row below it. The cap directions are the better,
per k, of the best **lattice** and the **record** kissing configuration of R^k; the record wins
exactly where it splits perfectly into zero-sum triples — 306 = 3·102, 510 = 3·170 and
1932 = 3·644 — and loses where it is a layered non-lattice packing whose triple structure is
too thin for the 4/3 factor. Dimension 83 (k = 11) used to be carried from dimension 82, because Λ₁₁'s 146
triples give only 292 units against dimension 82's 340; it now has a construction of its own
from the 604-point record configuration of ℝ¹¹ realised over ℤ[√2], whose 84 triples and 176
pairs give 344.

All fifteen configurations and partitions are re-checked exactly by
[`scripts/verify_parts.py`](scripts/verify_parts.py): each spans exactly k dimensions, is a
60-degree code, and its parts are disjoint with every triple summing to zero at 120 degrees.
Since 2026-08-30 **every one of the fifteen is perfect** — 3T = |W|, no pair parts and no
leftover vectors anywhere in k = 9…23. Dimension 83's ℤ[√2] configuration, checked separately
by [`scripts/verify11.py`](scripts/verify11.py), is the one exception, and provably so: 84
triples is its exact maximum, and 84 T + 176 P its exact optimum on 2T + P, both settled by
CP-SAT in [`scripts/k11_optimal.py`](scripts/k11_optimal.py).

## The two class families

Dimension 72+k reads the `T+P` largest classes of ONE family of pairwise-disjoint Γ₇₂ classes,
and two are shipped. They have opposite shapes, which is why both are kept:

| family | classes | sizes | total lines | how |
|---|---|---|---|---|
| `class_sizes_32000.npy` | 32 000 | 867–2106 | 47 150 230 | automorphism images of one greedy 2106-line class |
| `class_sizes_gpu.npy` | 32 000 | 1932–2118 | 65 570 602 | greedy on a 80 000 000-line pool of conflict density 0.004827 (`scripts/gpu/`) |

The image family's first few thousand classes are excellent and only its tail collapses; the
fresh family is nearly flat. So which one wins depends on how deep the dimension reads, and the
crossover is at **dimension 86**:

| dim | classes read | images | fresh | fresh wins by |
|---|---|---|---|---|
| 81 | 102 | 214 591 | 214 323 | -0.1% |
| 82 | 170 | 357 367 | 356 826 | -0.2% |
| 83 | 260 | 545 983 | 545 227 | -0.1% |
| 84 | 252 | 529 237 | 528 491 | -0.1% |
| 85 | 306 | 642 224 | 641 446 | -0.1% |
| 86 | 644 | 1 346 134 | 1 347 578 | +0.1% |
| 87 | 780 | 1 627 715 | 1 631 359 | +0.2% |
| 88 | 1 440 | 2 981 444 | 3 006 643 | +0.8% |
| 89 | 1 782 | 3 674 580 | 3 718 386 | +1.2% |
| 90 | 2 466 | 5 043 564 | 5 140 391 | +1.9% |
| 91 | 3 556 | 7 179 645 | 7 403 375 | +3.1% |
| 92 | 5 800 | 11 402 771 | 12 052 663 | +5.7% |
| 93 | 9 240 | 17 447 448 | 19 161 262 | +9.8% |
| 94 | 16 632 | 28 881 727 | 34 378 537 | +19.0% |
| 95 | 31 050 | 46 221 022 | 63 699 093 | +37.8% |

`final81-95.py` takes the better of the two per dimension, by the rule stated there. They are
never mixed, since disjointness ACROSS families is not claimed. Declining the GPU family
entirely costs no RECORD: on the image family alone, which is regenerable on a CPU, every
dimension 81–95 still clears its published value, the smallest margin being 858 364 at
dimension 81.

**What it costs to check.** The image family is regenerable on a CPU: `scripts/verify_rebuild.py`
rebuilds a prefix from `data/class_best_2106.npy` and `data/gamma72_gens.npy`, verifies every
class exactly and checks the sizes against the shipped table. A 60-class prefix takes a minute;
`build_classes.py` then `materialise.py` then `verify_classes.py` do all 32 000 in a couple of
hours and about 1.6 GB of scratch. The fresh family is certified by its SEED rather than by its 4.6 GB of lines: every
random decision is numpy's PCG64 and the device does nothing but exact integer GEMM, so
`scripts/gpu/regenerate.py` rebuilds it and matches the SHA-256 — but that wants a GPU and a
couple of hours, and on a CPU it is not practical. Without one, `scripts/gpu/check_cert.py`
checks what needs no GPU in a second: that the shipped sizes are the certificate's, that the
disjointness invariant `|alive| + Σsizes = |pool|` holds, and that the measured `|⟨u,v⟩|`
histogram is consistent with the shell's EXACT distribution, which Venkov's theorem pins down.

Every class of the fresh family was verified in exact integer arithmetic as it was built —
each vector of norm 8, each off-diagonal `|⟨u,v⟩| ≤ 2` — and the worst value seen over all
31 000 of them was 2. That is why the lines themselves never had to be kept.

## The cap directions: what is shipped, and where it comes from

For each *k* the construction needs a 60° code of ℝᵏ split into parts. Three families were
tried and the better is taken per *k*; all are re-checked exactly by
[`scripts/verify_parts.py`](scripts/verify_parts.py) and, for *k* = 11,
[`scripts/verify11.py`](scripts/verify11.py).

| file | what it is | provenance |
|---|---|---|
| `lambda9_W.npy` | Λ₉, 272 vectors | built here: E₈ laminated at the deep hole (1,0,…,0), which has 16 lattice points at distance 1 |
| `lam10…lam16_W.npy` | Λ₁₀…Λ₁₆ | built here: BW₁₆ from RM(1,4), then descended one dimension at a time by the richest hyperplane |
| `lam17…lam23_W.npy` | Λ₁₇…Λ₂₃ | built here: Leech cross-sections by √2·R for R = E₇, E₆, D₅, D₄, A₃, A₂, A₁ |
| `rec_9_W.npy`, `rec_10_W.npy`, `rec_14_W.npy` | the record kissing configurations of ℝ⁹, ℝ¹⁰, ℝ¹⁴ (306, 510, 1932 points) | extracted from Cohn's published `dimensions1-24.txt`; the 1932-point one is the same configuration the dimension-38 package uses |
| `rec11_exact.npz` | the 604-point record configuration of ℝ¹¹, as `(A + B√2)/3` with A, B integer | same source, kept in exact ℤ[√2] form because its Gram is irrational |
| `kap12_W.npy`, `kap13_W.npy`, each with a `_c` | the Kappa sections **K₁₂** and **K₁₃** (756 and 918 vectors) in exactly *k* coordinates, with metrics diag(1,3)⁶ and diag(3,1,3,1,3,1,3,1,3,1,3,1,1) | derived by [`scripts/kappa_caps.py`](scripts/kappa_caps.py) from an order-3 automorphism of the Golay code; **larger than Λ₁₂ and Λ₁₃**, which is worth 36 and 4 triples |
| `pole_src_15.npy`, `pole_src_17.npy`, `pole_src_18.npy`, each with a `_c` | the record kissing configurations of ℝ¹⁵, ℝ¹⁷, ℝ¹⁸ (2564, 5730, 7654 points) with their inner-product coefficients, used as **pole** sources | extracted from the same `dimensions1-24.txt`, and compared with it by [`scripts/verify_record.py`](scripts/verify_record.py) whenever the data set is present |
| `parts_<k>.npz`, `rec_<k>_parts.npz`, `rec11_parts.npz` | the partitions themselves, as indices | computed here |

The record configuration wins exactly where it splits **perfectly** into zero-sum triples —
306 = 3·102, 510 = 3·170, 1932 = 3·644 — and at *k* = 11, where 84 triples plus 176 pairs give
344 units against Λ₁₁'s 292.

**The *k* = 11 certificate was corrected on 2026-08-30.** A pair part needs cos ≤ −½, which in
ℤ[√2] is `P + N/2 ≤ −Q√2`; [`scripts/verify11.py`](scripts/verify11.py) had that inequality
with its sign reversed, under an `or` whose other branch is a condition on the *conjugate*
`P − Q√2` and is satisfied by everything in that Gram. **169 of the 176 shipped pairs sat at
cos = +½ or cos = 0** and were not parts at all. The units 2T + P = 344 were nevertheless
right: 344 is the *exact* optimum over triples and pairs together and 84 the exact maximum
number of triples, both settled by CP-SAT in seconds, so dimension 83's value does not move.
`rec11_parts.npz` now holds an optimal partition whose 176 pairs are genuine — 172 at cos = −½
and 4 antipodal — and the checker demonstrates that its test can fail, on a cos = +½ and a
cos = 0 witness drawn from the same Gram, before it trusts it.

**Where it loses, it loses as a CAP set only.** For *k* = 15, 17, 18, 19, 20 the record
configuration is a layered non-lattice packing whose triple structure is too thin for the 4/3
factor, so the lattice supplies the cap directions. The measured comparison, reproducible
from `dimensions1-24.txt` with `research/dim81-95/recordW2.py`:

| k | record \|W\| | triples + pairs | units (record) | units (lattice) |
|---|---|---|---|---|
| 15 | 2564 | 508 + 511 | 1527 | **1560** |
| 17 | 5730 | 482 + 2120 | 3084 | **3564** |
| 18 | 7654 | 1475 + 1589 | 4539 | **4932** |
| 19 | 11948 | 684 + 4875 | 6243 | **7112** |
| 20 | 19448 | 760 + 8502 | 10022 | **11600** |

Dimension 12's 841-point record has **no zero-sum triples at all**.

The record column is a *packing*, so it was re-measured on 2026-08-30 with the walk of
[`scripts/oddparts.py`](scripts/oddparts.py) — the one that turned 1773 into 1782 at
*k* = 17 — in case the earlier greedy had been the limit. It had not: 520 triples against
508 at *k* = 15, 482 against 482 at *k* = 17, 1496 against 1475 at *k* = 18. At *k* = 17
the ceiling is structural — 384 of the 5730 points lie in **no** zero-sum triple, so the
packing cannot exceed 1782 there whatever the search — and the lattice keeps the column
by hundreds of units. This is a measurement of the configurations, not of the packer.

That table is about the *caps*. The **poles** are a different job: they have to be a 60°
code lying 30° from the caps and nothing else, and no partition into triples is asked of
them. So the very configurations this table rejects as cap sets are the best available pole
sources, and they are shipped as `pole_src_<k>.npy` — see
[`scripts/axis_cohn.py`](scripts/axis_cohn.py) below.

## What "previously" means, and why the range ends at 95

Above dimension 72 the published values are direct sums: **Γ₇₂ ⊥ (best k-dimensional
lattice rescaled to minimum 8)**, of size 6 218 175 600 + τ(k). The Nebe–Sloane table
records precisely that at dimension 80 ("Γ₇₂ ⊥ E₈") and records **nothing at all** for
73–79 or 81–127.

The only other candidate is Edel–Rains–Sloane, and reverse-engineering their count settles
it: it is dominated by the full-support term A(n₀, n₀/4), and codetables.de gives
A(80,20) ≥ 2³⁰, A(84,21) ≥ 2²⁸, A(88,22) ≥ 2³⁰, A(92,23) ≥ 2³⁰ — at most
2³⁰ = 1 073 741 824 all the way to n₀ = 92, with the remaining terms adding well under a
factor two (their published dimension-80 total is 1 368 532 064). So they stay far below
6 218 175 600 across 81–95. **A(96,24) ≥ 2³³ is where that stops** — codetables carries
[96,33,24], and the ERS chain (96,24,6,1) totals 12 886 999 232, a factor 2.07 on Γ₇₂ — which
is why this package ends at dimension 95. Dimensions 93, 94, 95 were re-evaluated the same
way and ERS gives 2.68e9, 2.68e9, 4.83e9 there, so the step falls between 95 and 96 and the
claims below are unaffected. See
[`../../closed/dim96-ers-takeover/`](../../closed/dim96-ers-takeover/).

## The idea

In ℝ⁷² ⊕ ℝᵏ, with every point a unit vector and every pairwise inner product ≤ 1/2, at cap
level *t*:

| family | points | count |
|---|---|---|
| equator | (v/\|v\|, 0), v minimal in Γ₇₂ and not ±u for a class line u | 6 218 175 600 − 2Σ\|C_i\| |
| caps | (s√t·û, √(1−t)·z), s ∈ {±1}, u a line of class *i*, z ∈ Z_i | 2Σ\|C_i\|·\|Z_i\| |
| poles | (0, w) | ≤ τ(k) |

giving `6218175600 + 2Σ|C_i|(|Z_i| − 1) + #poles`. The four binding constraints, with γ the
class threshold:

| pair | condition |
|---|---|
| same z, two lines of a class | t·γ + (1−t) ≤ 1/2 → t ≥ 1/(2(1−γ)) |
| same line, two directions | ⟨z,z′⟩ ≤ (1/2 − t)/(1 − t) |
| two classes, two directions | ⟨z,z′⟩ ≤ 1/2 |
| cap · pole | √(1−t)·⟨z,w⟩ ≤ 1/2 |

### The point: Γ₇₂ admits ZERO-SUM TRIPLES

Γ₇₂'s class condition is |⟨u,u′⟩| ≤ 2 of 8, i.e. **γ = 1/4** — exactly the value that
permits t = 2/3. And at t = 2/3 the second row becomes ⟨z,z′⟩ ≤ **−1/2**, so a class line
may be used at a **zero-sum triple** of directions rather than an antipodal pair. A triple
part is worth 4 points per line against a pair's 2.

**Three is the maximum possible**: n points pairwise at cosine ≤ −1/2 satisfy
n ≤ 1 + 1/(1/2) = 3. So the gain is `2·Σ|C_i|·(|Z_i| − 1)` with |Z_i| ≤ 3, and the best
possible partition of the direction set W into parts is into ⌈|W|/3⌉ zero-sum triples.

`scripts/triples.py` computes and verifies perfect zero-sum triple partitions of A₁, A₂,
A₃, D₄, D₅, E₆, E₇, E₈, reaching the optimum in every case that allows it:

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| factor \|W\| − M with **triples** | 1 | **4** | **8** | **16** | **26** | **48** | **84** | **160** |
| factor with antipodal pairs | 1 | 3 | 6 | 12 | 20 | 36 | 63 | 120 |

P₄₈ cannot do this: its γ is 1/3, which forces t = 3/4, where two directions on a line need
⟨z,z′⟩ ≤ −1. See [`../dim49-63-p48-caps/`](../dim49-63-p48-caps/).

### A caution about the direction sets for k = 19–23

Those rows take W to be the **A_j cross-section of the Leech**, which coincides with the
laminated Λ_k only for j ≤ 3 and is genuinely smaller after that — 8676 against Λ₁₉'s
10668 at k = 19, for instance. That is fine and deliberate: the direction set only has to be
a 60° code in ℝᵏ, and the gain is \|W\| − M, not \|W\|. Triples on the smaller set beat
pairs on the larger one — 5708 against 5334 at k = 19 — because a triple is worth 4 points
per class line and a pair only 2.

`final81-95.py` **asserts** this rather than assuming it: for every k where it chooses
triples, it computes what antipodal pairs on the full Λ_k would have given and checks the
triples are at least as good. At k = 23 the pair scheme is not even feasible — it would need
46 575 disjoint classes against the 32 000 available.

### Explicit classes: 2118 lines where Caro–Wei guarantees 221

Γ₇₂'s Gram matrix and six automorphisms are published in the Catalogue of Lattices (Nebe,
entry `Gamma72`). Verified here in exact integer arithmetic: 72×72, symmetric, diagonal all
8, and **determinant exactly 1** by fraction-free Bareiss elimination — hence even
unimodular.

Its minimal lines were generated by automorphism-orbit closure plus the production rule
*u, v minimal with ⟨u,v⟩ = 4 ⟹ u − v minimal*, giving a pool of 3.39 million lines, and
greedy plus local search on that pool found a class of **2106 lines** where Caro–Wei could
guarantee only 221. The greedy family built later on the 80 000 000-line pool holds a larger
one still, of **2118 lines** — a factor **9.58** over Caro–Wei, in line with the Leech's 11.3
and P₄₈'s 9.93. 2106 is the base class of the *image* family below, and the number that
appears in the file names.

Because a random automorphism image of a class is again a class and two images barely
overlap, **120 pairwise disjoint classes** (2096–2106 lines each, 252 358 in total) come
essentially free; those are shipped and verified. Scaling the same idea to **32 000
pairwise disjoint classes** (47 150 230 lines) is what dimensions 82–95 need.

### The calibration — the single most convincing check

`scripts/calibrate.py` runs this exact model over the **Leech** lattice, where every answer
is already in Cohn's table, with the real 248-line class:

| dim | 25 | 26 | 27 | 28 | 29 | 30 | 31 |
|---|---|---|---|---|---|---|---|
| this model | 197058 | 198550 | **200540** | 204520 | 209496 | 220440 | 238350 |
| Cohn's table | 197056 | 198550 | 200044 | 204520 | 209496 | 220440 | 238350 |
| | +2 | MATCH | +496 | MATCH | MATCH | MATCH | MATCH |

Five exact reproductions, and the two dimensions where it differs are **precisely the two
records this project already holds by independent, fully verified routes**
([`../dim25-cap-level/`](../dim25-cap-level/) and
[`../dim27-triple-partition/`](../dim27-triple-partition/)). The script also builds the
equator and caps in real coordinates at every k and confirms the largest pairwise inner
product is exactly 1/2.

*(This replaces `capcheck2.py` of the working tree, which crashed at k = 6 — E₆ and E₇ are
returned as subsets of the E₈ roots, so they live in ℝ⁸ rather than ℝ⁶ and ℝ⁷ — and which
used small greedily-found classes rather than the real 248-line one, so its totals fell
below the table and it never actually calibrated anything.)*

### Poles

At t = 2/3 a pole meets a cap at √(1/3)·⟨z,w⟩, which is 0.5774 > 1/2 if the pole sits at a
cap direction, so poles must be at least 30° from every cap direction. What is claimed is the
**axis layer**, and the legality argument is a shell argument: for a of norm M and w of norm m
in the same lattice, |a−w|² and |a+w|² are lattice norms, so |⟨a,w⟩| ≤ M/2, and the pole
condition 4⟨a,w⟩² ≤ 3Mm then holds automatically as soon as **M ≤ 3m**. So *every* shell of
the direction lattice with m < M ≤ 3m is legal at once, and both +a and −a are poles.

[`scripts/poles18.py`](scripts/poles18.py) used to take candidates of the form a = w₁ + w₂
with w₁ ⊥ w₂, which is a subset of one shell — and A₂ has no two orthogonal roots at all, so
dimension 74 claimed **no poles** while A₂'s norm-6 shell is legal and gives all six. Sweeping
both legal shells and searching each properly gives **374 poles against 248**: k = 2 goes
0 → 6 and k = 8 goes 98 → 174. The enumeration is checked against the theta series of each
root lattice and every set is re-verified from scratch, so this is exact rather than a
heuristic.

### The cap set is a choice, and a group can make it (2026-08-30)

The construction asks of the direction set *W* only that it be a 60° code splitting into
parts internally at cosine ≤ −½ — not that it be any particular code. A part of three is a
zero-sum triple, worth **4 per class line** against a pair's 2, and a class here is about
2100 lines: **one missing triple is about 8400 spheres**, where a missing pole is one.

Suppose an isometry *S* of the axis space satisfies

    I + S + S² = 0.

Then *S* has eigenvalues ω and ω̄ only, so it has no fixed vector and acts **freely**; and for
every *v*, the orbit {*v*, *Sv*, *S²v*} sums to (I + S + S²)*v* = 0, so it is a zero-sum
triple outright. If *W* is *S*-invariant, its partition into orbits is a **perfect** triple
partition — found by a group, not packed. Those eigenvalues come in conjugate pairs, so this
can only happen in **even** dimension; *k* = 9, 11, 13, 15, 17, 19, 21 and 23 are ruled out
by that alone, and are closed instead by
[`scripts/oddparts.py`](scripts/oddparts.py) below.

**K₁₂ and K₁₃.** Let π be an automorphism of the Golay code of cycle shape 1⁶3⁶ — an element
of order 3 in M₂₄ — acting on the Leech by permuting coordinates. Six fixed points and six
3-cycles give twelve orbits, so its fixed space is 12-dimensional, and so is the orthogonal
complement; on the complement π has no eigenvalue 1, so I + π + π² = 0 there. The 756 Leech
minimal vectors lying in it have determinant 3⁶ at minimal norm 4 — the Coxeter–Todd lattice
K₁₂ — and their 252 zero-sum triples are the orbits of π. Adjoining the one further
direction that carries the most minimal vectors gives a 13-dimensional section of 918;
13 is odd, so those 306 triples are packed rather than given.
[`scripts/kappa_caps.py`](scripts/kappa_caps.py) derives both, and asserts the determinant,
because a count does not identify a lattice.

**k = 16, 18, 20 and 22.** [`scripts/omega_parts.py`](scripts/omega_parts.py) searches for
such an *S* directly, inside the axis space, for cap sets never suspected of having one. It
finds one at all four, and the partitions go from 1434 triples + 7 pairs to 1440, from
2457 + 11 to 2466, from 5782 + 23 to 5800, and from 16 594 + 49 to 16 632. At *k* = 20 the
isometry condition alone leaves the first level of the search 664 ways wide; *S* carries
zero-sum triples to zero-sum triples, so it preserves the number of them through a vector,
and filtering candidate images by that degree narrows the level to eight. At *k* = 22 that
degree is constant at 1492 and the filter buys nothing — the isometry condition settles it
there on its own, in twelve complete assignments. Consistency on a basis is not
set-preservation: the first complete assignment the search reaches is set-preserving at
*k* = 18 and is **not** at *k* = 12 or 16, so every complete assignment is tested by applying
*S* to the whole of *W*.

**The odd *k*, where no such *S* exists.** [`scripts/oddparts.py`](scripts/oddparts.py) closes
them by search, but by a search with the right move. The zero-sum triple hypergraph is
**linear** — two vectors lie in at most one triple, since the third is forced to be −(a+b) —
and that pins the move set down completely. If a triangle {*u*,*v*,*w*} has *u* uncovered and
*v*, *w* in a common part, that part's third vector is −(*v*+*w*) = *u*, which is uncovered;
so it never happens, evicting one part always leaves the uncovered count **unchanged**, and

> every neutral or improving move is a **pair of uncovered vectors at 120°**.

With |U| in the tens that is a few hundred pairs, enumerated in one pass. A sampler that draws
a random triangle through a random uncovered vector is looking for a needle — at *k* = 23 each
vector lies in 1232 triples and all but a handful lead nowhere. All four partitions come out
perfect: **1782, 3556, 9240 and 31 050 triples**, worth 286 996 spheres across dimensions 89,
91, 93 and 95.

**One parameter decides which move the walk is using.** The pair scan is quadratic in |U|, so
the walk enumerates below a threshold and samples above it, and the first threshold written
was 150. At 150, *k* = 23 stalls at 31 047 or 31 048 — through cold restarts, through iterated
local search with kicks of up to forty triples, through a tabu on recently placed triples, and
through seeding the walk with Λ₂₂'s own perfect partition, which sits inside Λ₂₃ as the level-0
set of the functional orthogonal to its span. Four strategies, one answer, and it reads as a
property of the lattice. At 800 it closes in 55 seconds. Sampling once a packing is nearly full
is exactly where sampling is worst.

**The symmetry that is not available.** *W* = −*W* and a triple's signs are fixed up to a
global one, so its triples come in antipodal pairs and the problem looks like it halves; an
*S*-invariant partition is antipodally symmetric for the same reason, the orbit of −*v* being
the negative of the orbit of *v*. Summing [*a*] + [*b*] + [*c*] = 0 in Λ/2Λ over such a
partition gives Σ over the |*W*|/2 minimal **lines** of [*v*] = 0, and that fails for Λ_k at
*k* = 11, 13, 15, 17, 19 and 23. On the vectors the same sum is doubled and the condition is
vacuous. A line-level walk therefore stalls one triangle short at *k* = 17, at 890 of 891, for
arithmetic reasons rather than for want of looking — the `sym` column of `oddparts.py` reports
the test, and the four partitions it ships are not antipodally symmetric and cannot be.

**A layer belongs to a cap set.** An axis layer is a rotation *N*/*D* together with a set of
*indices into W*, and both halves are meaningless against a different *W*.
[`scripts/capsha.py`](scripts/capsha.py) stamps each layer with the digest of the set it was
built from and readers check it; the layers of the Kappa sections are named apart as
`poles_kap_<k>.npz`; and the Λ₁₂ and Λ₁₃ layers were retired rather than left on disk, where
a 648-pole layer read against a 756-vector set is a well-typed answer to a question nobody
asked.

### The axis layer is a rotated copy of the cap directions (2026-08-27)

The shell rule can only offer vectors **of the lattice itself**, and that is what kept these
layers small: at k = 23 it found 248 poles against a ceiling of τ(23) = 93 150, four
thousandths of a per cent, and this README used to say that the rest was "worth nothing
either way".

A pole has to be a 60° code and to clear the 30° caps around the cap directions. Take it to
be a **rotated copy of the cap direction set itself** and the first condition is free, because
a rotation preserves inner products; only the second is left. This is dimension 38's question
([`../dim38-leech-large-codimension/`](../dim38-leech-large-codimension/)) and here it is much
easier, because in these dimensions the caps are small: a single 30° cap covers
2.3×10⁻⁸ of S²², so all 93 150 of them together cover at most 2.1×10⁻³, and a
Haar-random rotation of the 93 150 directions of Λ₂₃ is expected to lose at most
about two hundred.

[`scripts/axis_rotate.py`](scripts/axis_rotate.py) takes **K** an integer skew matrix with
entries in {−1,0,1}, puts *S* = *C*⁻¹**K** for the metric *C* = diag(*c*), and forms the
Cayley transform *Q* = (*I*−*S*)(*I*+*S*)⁻¹. Then *Q*ᵀ*CQ* = *C* exactly, so *Q* is a rational
isometry of the form; its denominator *D* comes out with one to eighteen digits — the largest at
k = 14, which borrows dimension 38’s rotation, and the smallest where sampling a small
integer skew matrix happens to find one — and the test

> 4⟨a,z⟩² ≤ 3m²  ⇔  4X² ≤ 3m²D²  ⇔  2|X| ≤ isqrt(3m²D²),  X = (Nw)ᵀCz

stays in `int64` throughout — the last step because 3m²D² is never a perfect square.

**That rotation has to be an isometry of ℝ^k, and for seven of these k it was not.** The cap
direction sets are written in the coordinates of the lattice they are cut from, and several of
them span fewer: Λ₂₁'s cross-section is 27 720 vectors written in Λ₂₄'s 24 coordinates and
spanning 21 of them. A Cayley transform of those *coordinates* is an isometry of ℝ²⁴, not of
ℝ²¹, and it carries the poles out of the axis space — which is what the layers written on the
morning of 2026-08-27 did at k = 13, 15, 17, 18, 21, 22 and 23, with every inner product in
them correct.

[`scripts/axis_block.py`](scripts/axis_block.py) derives those seven instead. For *v*, *w*
rows of *W*, the matrix *S* = *v*(*Cw*)ᵀ − *w*(*Cv*)ᵀ is integral, *CS* is skew, and *S*
annihilates everything *C*-orthogonal to *v* and *w* — so its Cayley transform fixes the
complement of the span **pointwise** and cannot move a pole out of ℝ^k. Over *r* mutually
orthogonal pairs scaled by 1/λ it is a direct sum of plane rotations through 2·arctan(μ/λ) and
can be written down without inverting anything:

> *N* = *D*·*I* + Σₜ [ −2μ(*v*(*Cv*)ᵀ + *w*(*Cw*)ᵀ) + 2λ(*v*(*Cw*)ᵀ − *w*(*Cv*)ᵀ) ],  *D* = λ² + μ²

with denominator 80 to 5120. The obvious alternative — work in a lattice basis of the span,
carrying the metric across as *G* = *PCP*ᵀ — is correct and unusable: *S* = *G*⁻¹**K** has
denominators dividing det(*G*), which is 2⁵⁹ at k = 17 and 2⁷¹ at k = 23, and det(*G*) is a
lattice invariant, so no basis avoids it. That was measured, not argued: no usable rotation
exists in that family at any of k = 17, 18, 21, 22, 23.

`verify_poles.py` now checks the space as well as the angles, on a basis of the span — which
settles it, the rotation being linear — and `axis_rotate.py` rejects a candidate that fails, so
it can no longer write such a layer. The shell layers were never affected: their poles are
vectors of the same lattice and lie in the span by construction, which was checked. The
search is a search; what makes each layer a certificate is
[`scripts/verify_poles.py`](scripts/verify_poles.py), which takes the shipped *N*, *D* and
index set and checks every pair in exact integers.

### The poles do not have to be a copy of the caps

A rotated copy of *Z* cannot exceed |*Z*|, and at k = 19, 20 and 21 the cap directions are the
minimal vectors of Λ_k while τ(k) is larger — 10 668, 17 400 and 27 720 against 11 948,
19 448 and 29 768. Those three layers were short **by choice of ansatz**, not by geometry:
the construction asks of *W* only that it be a 60° code lying 30° from *Z*, and nothing in
that asks it to be a copy of anything.

[`scripts/axis_record.py`](scripts/axis_record.py) takes the **record configuration of ℝ^k**
as the source instead. At k = 19 and 20 it is already in the frame the cap file is written
in: `data/lam19_W.npy` and `data/lam20_W.npy`, halved, are integer vectors of norm 8 in
exactly *k* coordinates of two shapes — 4·C(n,2) vectors (±2,±2,0,…,0), and 128 sign patterns
on each octad of a shortened Golay code, always with an **even** number of minus signs. That
is Λ_k, and it is also the frame Cohn–Li work in. Take the octad signs **odd**: the count is
the same, the code is still 60°, and now *s* = ((−1)^c₁,…,(−1)^cⁿ) may be adjoined for every
*c* in the dual of the octad code, since ⟨s,u⟩ = ±8 would need *s* to agree with an odd
pattern on an octad and *s* has an even number of minus signs there. Two such vectors are
compatible when their Hamming distance is at least ⌈n/4⌉, so the extra term is a maximum
independent set in a Cayley graph on the dual code — 1 280 at k = 19, 2 048 at k = 20 and 21.
The sources are then 11 948, 19 448 and 29 768, which is τ(k) exactly.

**One inequality does most of the work, and it is not the one the search was built for.** A
pole needs only 30° from the caps, not 60°. At k = 19 and 20, where the source shares the cap
file's frame, both new families clear that with *no rotation at all*: in the halved frame an
odd-sign octad vector meets a cap at |X| ≤ 6 against a limit of 6, and a ±1 vector at
|X| ≤ 8 against a limit of 10. Only the 4·C(n,2) pair vectors fail, and they fail completely
— they **are** cap directions, so all 684 of them sit at 0°. The unrotated layer is therefore
already 9 984 + 1 280 = 11 264 at k = 19 against the 10 648 the rotated copy gave, with
denominator 1. A rotation moves the pair vectors off their own caps and is worth the rest.

**At k = 21 the source has to be carried in, and whether it can be is a fact about the
space.** The cap directions there are 27 720 vectors in 24 coordinates spanning 21, so a
configuration built in ℝ²¹ needs a rational similarity *M* with *M*ᵀ*M* = μ*I* — *k*
orthogonal vectors of one norm in the span. The **discriminant** decides: the determinant of
any Gram of the span, modulo rational squares, is the product of the norms of any orthogonal
basis, so a frame of common norm μ forces disc = μ^k.

| k | 12 | 13 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| discriminant | 1 | 2 | 1 | 1 | 2 | **3** | 1 | 1 | 1 | **3** | 2 |

At k = 18 and 22 the rank is even, so μ^k is a square while the discriminant is 3: **no frame
of one norm exists at all**, and those two axis spaces cannot hold a rational copy of anything
built in the *standard* ℝ^k. That is a statement about the metric, not about the
configuration, and dimension 90 turns out not to need the standard one — see the next
section. At k = 21 the discriminant is 1, so μ must be a square — which is why a
greedy frame drawn from the norm-32 cap directions stops at 20 of 21 however long it runs, and
why one drawn from norm-4 vectors of the span completes at once. This also sidesteps the
obstruction recorded above for a lattice basis of the span: in frame coordinates the metric is
μ*I*, of determinant μ^k, not the lattice's 2⁷¹.

**Dimension 83 was a ring problem, and it said so.** Its 604 cap directions live over
ℤ[√2], and the layer was empty because "the rational rotations above are not set up for it".
Only the *configuration* is irrational; the *rotation* need not be.
[`scripts/axis_sqrt2.py`](scripts/axis_sqrt2.py) takes *Q* = *N*/*D* rational with
*N*ᵀ*CN* = *D*²*C*, so the poles stay in the same ring, and the 30° condition becomes
(*p* + *q*√2)² ≤ 972*D*² for integers *p*, *q* — decided by two signs and *a*² against 2*b*².
The first random start cleared all 604, at 0.9435 of the threshold.

**The same idea reaches four more dimensions, and there the configuration is published
rather than reconstructed.** At k = 13, 15, 17 and 18 the cap file is Λ_k written in more
coordinates than it spans, and the record configurations are not of the odd-octad shape at
all — but Cohn's `dimensions1-24.txt` gives them as exact integer coordinates together with
an integer metric diag(c), so the source is on disk and the only question is whether the axis
space can hold it. [`scripts/axis_cohn.py`](scripts/axis_cohn.py) answers it.

**The metric is what makes dimension 90 work.** A frame now has to satisfy
*M*ᵀ*M* = μ·diag(c), so the discriminant condition reads μ^k·det(diag c) ≡ disc, and the
record configurations are not written in the standard metric:

| k | metric of the record | disc of that metric | disc of the axis space | μ | frame |
|---|---|---|---|---|---|
| 13 | diag(1¹³) | 1 | 2 | 8 | 13 of 13 orthogonal caps |
| 15 | diag(1¹⁴, 2) | 2 | 1 | 8 | 14 of 15 orthogonal caps + the complement line |
| 17 | diag(1¹⁶, 2) | 2 | 2 | 64 | 17 of 17 orthogonal caps |
| 18 | diag(1¹⁶, 2, 6) | 12 | 3 | 64 | 17 of 18 orthogonal caps + the complement line |

At k = 18 the standard metric is impossible — that is the paragraph above — and diag(1¹⁶, 2, 6)
has discriminant 12 ≡ 3 modulo squares, which is exactly the discriminant of the axis space.
The obstruction was never to the configuration. It was to the metric, and this configuration
does not use that one.

**The frames are built in closed form from mutually orthogonal cap directions**, as many as
the discriminant allows — fourteen of fifteen at k = 15, seventeen of eighteen at k = 18. At
k = 17 and 18 the columns are *v*₁ ± *v*₂, …, *v*₁₅ ± *v*₁₆ of norm 64 and 2*v*₁₇ of norm 128;
where the count falls one short the last column is the orthogonal complement of the chosen
caps inside the span, which is a **line**, so its norm class is forced rather than chosen. It
comes out at 1 modulo squares at k = 15 where the frame needs 1, and at 6 at k = 18 where the
frame needs 6.

**The rotation is a product of rational plane rotations, and that is what keeps the
denominator small.** A full-rank Cayley transform on a rank-17 space has denominator about
*d*¹⁷, so a usable *D* allows *d* = 4 — 12% per entry, against the 1% these margins need. One
plane rotation is different: in the plane of two frame axes with metric weights *g_i*, *g_j*,

    R = [[a, -(g_j/g_i) b], [b, a]] / D,     a² + (g_j/g_i) b² = D²

is rational, is an isometry of diag(c), and has denominator *D* alone; the rational points of
that conic are dense in the circle. A greedy pass over planes and angles — each round scored
on a subsample, the best few counted in full — does the whole job in a handful of rotations.
At k = 15 and 17 the *first* one, with *D* = 13, is worth more than the entire shipped layer's
shortfall.

| dim | k | layer before | layer now | τ(k) | source |
|---|---|---|---|---|---|
| 83 | 11 | 0 | 604 | 604 | the 604-point record over ℤ[√2], rotated |
| 87 | 15 | 2 340 | 2 562 | 2 564 | the record configuration of ℝ¹⁵, through a frame |
| 89 | 17 | 5 338 | 5 704 | 5 730 | the record configuration of ℝ¹⁷, through a frame |
| 90 | 18 | 7 374 | 7 624 | 7 654 | the record configuration of ℝ¹⁸, through a frame |
| 91 | 19 | 10 648 | 11 934 | 11 948 | the record configuration of ℝ¹⁹, rotated |
| 92 | 20 | 17 380 | 19 415 | 19 448 | the record configuration of ℝ²⁰, rotated |
| 93 | 21 | 27 682 | 29 712 | 29 768 | the record configuration of ℝ²¹, through a frame |

The search is a search. What makes each layer a certificate is
[`scripts/verify_record.py`](scripts/verify_record.py) and, at k = 11,
[`scripts/verify11.py`](scripts/verify11.py): the source is **rebuilt** from the cap file or
the Golay code and the shipped code words, checked to be a 60° code, checked to have τ(k)
points, and every kept pole is checked against every cap in exact integers.
[`scripts/seal_record.py`](scripts/seal_record.py) stamps each file with the digest of the
source its index set was made for, so a rebuild in a different order cannot pass unnoticed.

| dim | k | shell layer | rotated layer | τ(k) | fraction of the ceiling |
|---|---|---|---|---|---|
| 81 | 9 | 102 | 306 | 306 | **at the ceiling** |
| 82 | 10 | 116 | 510 | 510 | **at the ceiling** |
| 83 | 11 | 0 | 604 | 604 | **at the ceiling** |
| 84 | 12 | 0 | 756 | 841 | 89.9% |
| 85 | 13 | 0 | 918 | 1 154 | 79.5% |
| 86 | 14 | 420 | 1 932 | 1 932 | **at the ceiling** |
| 87 | 15 | 802 | 2 562 | 2 564 | 99.9% |
| 88 | 16 | 1 546 | 4 320 | 4 320 | **at the ceiling** |
| 89 | 17 | 1 228 | 5 704 | 5 730 | 99.5% |
| 90 | 18 | 1 060 | 7 624 | 7 654 | 99.6% |
| 91 | 19 | 996 | 11 934 | 11 948 | 99.9% |
| 92 | 20 | 822 | 19 415 | 19 448 | 99.8% |
| 93 | 21 | 612 | 29 712 | 29 768 | 99.8% |
| 94 | 22 | 392 | 49 852 | 49 896 | 99.9% |
| 95 | 23 | 248 | 93 074 | 93 150 | 99.9% |

That is **220 879 more spheres** across these dimensions, and the layer reaches τ(k) exactly in
dimensions 81, 82, 83, 86 and 88. Dimension 86 is the case worth naming: its cap directions are
*the same 1932-point code* dimension 38 uses, vector for vector, so dimension 38’s rotation —
the one found by descending on SO(14) rather than by

Where a layer is short of τ(k) it is short for one of two reasons, and both are measured rather
than argued. For k = 12 and 13 the ceiling τ(k) is larger than the cap direction set itself,
and a rotated copy cannot exceed the set it is a copy of — which is why the poles at k = 15,
17, 18, 19, 20 and 21 are a rotated copy of the RECORD configuration of R^k instead, and reach
within a few dozen of τ(k).

**That leaves k = 12 and 13 as the two places the same move is not available, and the reason
is not the one above.** The bound "a rotated copy cannot exceed the set it is a copy of" is a
statement about the *ansatz*, and dropping that ansatz is precisely what gained 1562 spheres
at k = 15, 17, 18, 19, 20 and 21. What actually shuts these two is arithmetic, and both were
re-measured on 2026-08-30:

| k | what the record configuration of ℝᵏ would give | why it cannot be carried in |
|---|---|---|
| 12 | τ(12) = 841 against the shipped 756 — at most **+85** on dimension 84 | Cohn's data set writes all 841 points in **floating point**; not one is integral, so there is no rational pole source to place. Nothing in this repository decides anything in floating point, and a pole layer is no exception |
| 13 | τ(13) = 1154 against the shipped 918 — at most **+236** on dimension 85 | its 1106 integral points are written in the standard metric diag(1¹³); K₁₃'s span carries diag(3,1,3,1,3,1,3,1,3,1,3,1,1), and the two forms have **different Hasse invariants at 2 and at 3**. Rank 13 is odd, so μ fixes the discriminant and cannot touch the Hasse invariant — no rational similarity exists, for any μ |

Dimension 85 *did* carry 1106 poles until 2026-08-30, when its cap set changed from Λ₁₃ to
K₁₃. That trade cost 224 poles and bought four more zero-sum triples, which is 32 350 spheres
against 224 — a triple is worth four per class line and a pole is worth one
([`scripts/axis_cohn.py`](scripts/axis_cohn.py) records the frame that did not survive it).
So the 321 spheres in the table above are the exact size of what is left here, and both
entries are closed by arithmetic rather than by search. For k = 15, 17, 18, 19, 20, 21, 22 and 23 the rotation additionally
loses a handful of points to the caps; dimension 38 shows that descent on SO(k) closes exactly
this gap, and it is affordable at the small k and not at k = 23, where one gradient step costs
8.7×10⁹ inner products.

The same trick reaches **dimensions 75–80** as well, through the root systems. There the
metric is the Cartan matrix rather than a diagonal, so the rotation is a *G*-isometry,
built the same way from S = *G*⁻¹**K**; roots are coordinate vectors in the simple-root
basis, which spares E₆ and E₇ any embedding in ℝ⁸.

| dim | k | system | shell | rotated | τ(k) | |
|---|---|---|---|---|---|---|
| 75 | 3 | A₃ | 8 | **12** | 12 | **at the ceiling** |
| 76 | 4 | D₄ | 24 | 24 | 24 | already at the ceiling |
| 77 | 5 | D₅ | 28 | **40** | 40 | **at the ceiling** |
| 78 | 6 | E₆ | 58 | **72** | 72 | **at the ceiling** |
| 79 | 7 | E₇ | 76 | **126** | 126 | **at the ceiling** |
| 80 | 8 | E₈ | 174 | **240** | 240 | **at the ceiling** |

— **+146 spheres** over the shell rule, and the layer reaches τ(k) at k = 3, 4, 5, 6, 7 and 8.
At k = 3 the 30° caps cover four fifths of S², and SAMPLING a rotation there kept 6 against the
shell rule’s 8; descent finds all 12, which is τ(3). A file is written only where the rotation
wins, so nothing regresses. Dimension 74 (A₂, six roots, six poles) was already at τ(2).

k = 11 has no layer at all: dimension 83 uses the 604-point record realised over ℤ[√2], whose
Gram is irrational, and neither rule is set up for it.
[`scripts/verify_poles.py`](scripts/verify_poles.py) checks every layer exactly, against the
same cap directions each claim uses, and the counts are the `+n poles` column of the table.
Nothing beyond the axis layer is claimed. For k = 1 no pole direction exists at t = 2/3, so
dimension 73 uses the pair scheme at t = 3/4 instead, where √(1−t) = 1/2 makes poles automatic.

## History of the lower bounds in dimensions 73–95

| when | value | who |
|---|---|---|
| 1998 | 1 368 532 064 (dim 80) | Edel–Rains–Sloane [EdRS98] — the only published construction evaluated in this range, and a factor 4.5 *below* Γ₇₂ |
| 2012 | 6 218 175 600 | Nebe [Neb12] constructs Γ₇₂; dimensions 73–95 inherit it |
| ≤2011 | 6 218 175 840 (dim 80) | Nebe–Sloane table, "Γ₇₂ ⊥ E₈" — the only entry in the whole range |
| — | 6 218 175 600 + τ(k) | the direct-sum floor, mostly never written down |
| this work | **6 218 179 814 … 6 472 914 688** | the cap construction over Γ₇₂ with explicit classes and zero-sum triples |

**No prior claim in any of these twenty-three dimensions is known to us**, other than the
dimension-80 direct sum. The range has not been looked at because the tables stop.

The dimension-80 gain grew by a factor **12.70** during this project: 9.5× from the explicit
class (2106 lines instead of Caro–Wei's 221) and 4/3 from the triples.

## Files

```
final73-80.py                     dimensions 73-80, from the 120 shipped classes
final81-95.py                     dimensions 81-95, from the 32000-class size table
README.md                         this file

data/gamma72_gram.npy             72x72 Gram matrix, Catalogue of Lattices entry Gamma72
data/gamma72_gens.npy             6 published automorphism generators
data/class_best_2106.npy          the base class: 2106 minimal lines
data/disjoint_classes.npz         120 pairwise disjoint classes, 252358 lines  [7.5 MB]
data/class_sizes_32000.npy        sizes of the 32000-class family
data/lam_factors.npy              lamtriples.py's output: (k, |W|, |W| - M) for k = 9..23
data/kap12_W.npy, kap12_c.npy     K_12's 756 minimal vectors in twelve coordinates with the
                                  metric diag(1,3)^6, and kap12_parts.npz's 252 triples
data/kap13_W.npy, kap13_c.npy     K_13's 918 in thirteen, and kap13_parts.npz's 306 triples
data/poles_kap_12.npz, _13.npz    their axis layers, each carrying the digest of the cap set
                                  it was built from
data/lambda9_triples.npy          90 disjoint zero-sum triples over Lambda_9's 272
                                  vectors, covering 270: checked to sum to zero exactly
data/hash_r.npy                   the random vector build_classes.py hashes lines with

scripts/calibrate.py              THE KEY CHECK: the model over the Leech, against Cohn's table
scripts/triples.py                zero-sum triple partitions of A_1 ... E_8, verified
scripts/lamtriples.py             triple partitions for k = 19..23
scripts/kappa_caps.py             the cap sets of 84 and 85: the Kappa sections K_12 and
                                  K_13, larger than Lambda_12 and Lambda_13 and split
                                  perfectly into zero-sum triples by an order-3 symmetry
scripts/kappa_poles.py            their axis layers, as rotated copies with one plane
                                  rotation per Eisenstein block
scripts/omega_parts.py            perfect triple partitions from an isometry with
                                  I + S + S^2 = 0, where the packing search fell short
scripts/oddparts.py               perfect triple partitions for the ODD k, where no such
                                  isometry can exist: a plateau walk whose whole move set is
                                  the pairs of uncovered vectors at 120 degrees
scripts/capsha.py                 the identity of a cap set as a digest, so a layer cannot
                                  be read against the set it was not built from
scripts/verify_classes.py         exact integer verification of any saved class file
scripts/parse_neb72.py            catalogue HTML -> Gram + generators, with exact checks
scripts/verify_parts.py           the cap partitions: each spans k dimensions, zero-sum
scripts/verify_parts.log          that script's own output, shipped so the check is readable
                                  without rerunning it
scripts/verify11.py               the k = 11 partition and its axis layer, in Z[sqrt2]
scripts/k11_optimal.py            that partition is OPTIMAL, not merely the best found:
                                  CP-SAT settles max 2T+P = 344 and max T = 84 exactly
                                  (needs OR-Tools; not part of run_all.py)
scripts/verify_poles.py           the pole layer, against the shipped candidates
scripts/axis_record.py            the axis layers of 91-93, from the RECORD configuration
scripts/axis_sqrt2.py             the axis layer of dimension 83, over Z[sqrt2]
scripts/verify_record.py          those layers, exactly: source rebuilt, every pair checked
scripts/seal_record.py            re-derives each record layer and stamps its source digest
scripts/verify_rebuild.py         rebuilds the image family on a CPU and re-checks it
scripts/poles18.py                the retired pole search, kept for the record
scripts/make_readme_table.py      regenerates the 81-95 rows above from final81-95.py's own
                                  output, so the table cannot drift from the driver again

scripts/gamma72.py                [pipeline stage] minimal-line generation from the Gram
scripts/orbit_classes.py          [pipeline stage] 120 disjoint classes as automorphism images
scripts/build_classes.py          [pipeline stage] scale to 32000, storing (automorphism, mask)
scripts/materialise.py            [pipeline stage] expand that state into explicit line arrays

scripts/verify_all32000.py        [needs scripts/classes] global disjointness of all 32000.
                                  The original verify.log covered 31050; the family was later
                                  resumed to 32000 and that extension was never re-checked,
                                  yet dimension 95 draws on 31072 of them
scripts/verify_ext.py             [needs scripts/classes] the 950 extension classes
                                  (31050..31999) one at a time, in full: verify.log's exact
                                  sample was 250 drawn from the first 31050, and dimension
                                  95 uses 497 of the extension

scripts/gpu/run.py                builds the classes on a GPU, emitting a compact certificate
scripts/gpu/check_cert.py         re-checks that certificate without a GPU
scripts/gpu/checkpool.py          the pool a certificate was built from has no repeated line
scripts/gpu/fused.py              the fused inner-product kernel the builder runs on
scripts/gpu/regenerate.py         the full GPU rebuild, matched against its SHA-256
```

**`scripts/verify_all32000.py` and `scripts/verify_ext.py` are the two jobs `run_all.py`
reports as skipped**, in the `--full` set. They read the materialised classes from
`scripts/classes/`, which are 1.6 GB and are not shipped; both exit `SKIP --` when that
directory is missing or holds fewer than 32 class files, rather than pretending to have
checked anything. To run them, build it first with `scripts/build_classes.py` then
`scripts/materialise.py` — about an hour, and a deterministic function of three small files
that *are* shipped. `scripts/verify_rebuild.py` is not a substitute for either: it checks a
prefix of the family (60 classes by default, about a minute) and answers the narrower
question of whether the family really is regenerable to the shipped sizes.

**The stages marked `[pipeline stage]` need the output of the stage before them**, and those
intermediates (a 3.39-million-line pool, then several gigabytes of explicit classes) are not
shipped. Run bare, they exit with `FileNotFoundError` on a missing `.npy` — expected, not a
defect. `final73-80.py`, `final81-95.py`, `scripts/calibrate.py`, `scripts/triples.py` and
`scripts/verify_classes.py` all run from what is here.

### What is shipped, and what is regenerated

- **Dimensions 73–80 are fully certified by the shipped data.** `data/disjoint_classes.npz`
  holds all 120 classes explicitly; `scripts/verify_classes.py` checks in exact integer
  arithmetic that every vector has norm 8, every class has maximum |off-diagonal inner
  product| ≤ 2, and all 252 358 lines are distinct with no antipodal collision.
- **Dimension 80 is claimed 20 below what the larger family would give, on purpose.**
  The two families' sorted sizes are both shipped, so the comparison is a two-line
  calculation, and it does not go one way: at *n* = 78 and 79 the 120 classes are better
  by 8 and 16, and at *n* = 80 the 32 000-class family is better by 20. Dimensions 73–80
  take the 120 throughout, because those classes are shipped and checked in exact integer
  arithmetic one by one, whereas the 32 000 are regenerable rather than shipped and were
  verified in full on a sample. Twenty points out of 6.2 billion is not worth resting a
  dimension on a weaker certificate.
- **Dimensions 82–95 need 32 000 classes**, which is 47 150 230 lines — several gigabytes.
  Those are *not* shipped. `scripts/build_classes.py` rebuilds them deterministically from
  `data/class_best_2106.npy` and `data/gamma72_gens.npy`, storing each class compactly as
  (automorphism, survivor bit-mask); `scripts/materialise.py` expands that to explicit
  lines. Exactness is guaranteed rather than assumed: coefficients of Γ₇₂ minimal vectors
  never exceed 18 in absolute value, which is asserted for every row written, so the int8
  storage and the int64 line hash are both exact (the hash is |x·r| with |r_i| ≤ 2⁵¹, giving
  |x·r| ≤ 72·18·2⁵¹ = 2.9·10¹⁸ < 2⁶³). A hash collision could only cause a line to be
  *deleted* when it need not have been, so **the output is disjoint unconditionally**; only
  the sizes could suffer, by an expected ~1 line in 65 million. The original run verified a
  sample of 250 of the 32 000 classes in full exact arithmetic.

### The provenance chain

```
Catalogue of Lattices entry Gamma72
  -> scripts/parse_neb72.py    Gram + 6 generators; symmetric, diagonal 8, det = 1 (Bareiss)
  -> scripts/gamma72.py        3.39M minimal lines by orbit closure + the u-v rule
  -> greedy + local search     a class of 2106 lines        <- data/class_best_2106.npy
  -> scripts/orbit_classes.py  120 disjoint classes         <- data/disjoint_classes.npz
  -> scripts/build_classes.py  32000 disjoint classes       <- sizes shipped, regenerable
  -> final73-80.py, final81-95.py
```

## References

See [`../../../common/published.py`](../../../common/published.py).
