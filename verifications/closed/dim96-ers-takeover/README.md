# Dimension 96 — where Edel–Rains–Sloane takes over, and the cap construction cannot follow

**Dimension 96 is claimed here, at 12 886 999 232.** The mechanism that carries dimensions
73–95 reaches only **6 480 558 568** in dimension 96, so what establishes the claim is not
this repository's cap construction but Edel–Rains–Sloane's chain (96, 24, 6, 1) — exactly
as in dimensions 62 and 63. The floor it beats is the direct sum Γ₇₂ ⊕ Λ₂₄ = **6 218 372 160**,
a factor 2.07.

The package stays in `closed/` because that is also what it establishes: dimension 96 is
CLOSED to the cap construction, which is why the range of *layered* bounds still ends at 95.
Until 2026-08-26 this package recorded 12 886 999 232 as the FLOOR at 96 rather than as the
claim, which removed the dimension from the range by definition — a baseline set to your own
new result turns a claim into a non-result, and the improvement cannot appear because the
thing it improves on IS it. Nothing published has an entry at 96: Cohn’s table stops at 48
and 72, the Nebe–Sloane table it replaced listed 1–40, 42, 44, 48, 64, 72, 80 and 128, and
Edel–Rains–Sloane’s own worked examples are dimensions 32, 36, 40, 44, 64, 80 and 128
(arXiv:math/0207291 §3, read 2026-08-26).

**Read [`ers_sweep.py`](ers_sweep.py) first.** Answering the dimension-96 question meant
asking how big `A(n₀, ⌈n₀/4⌉)` is, and that question has an answer in every dimension. Asked
everywhere rather than at one point, it finds **dimensions 62 and 63 beaten by published
codes alone** — a result about this repository, not about dimension 96, and the more
consequential of the two. It is [below](#the-same-question-asked-everywhere-dimensions-62-and-63-are-beaten).

```
python build_triples24.py   # ~3 min: zero-sum triple partition of the Leech's 196560 vectors
python cap96.py             # the Gamma_72 cap construction at k = 24
python ers96.py             # Edel-Rains-Sloane at n = 93..96 with the 2026-08-26 tables
python table48-96.py        # tau(n) for every n from 48 to 96, with its source
python ers_sweep.py         # THE SAME QUESTION IN EVERY CLAIMED DIMENSION -- read this one
python ers_exposure.py      # the FULL chain in every claimed dimension, and each exposure
python above96.py           # why the table STOPS at 96, rather than at an arbitrary place
```

[`above96.py`](above96.py) answers the question a reader asks next. K is superadditive —
two 60° codes on orthogonal subspaces meet at inner product 0 — so K(96 + m) ≥ K(96) + K(m),
and every dimension above 96 already follows from dimension 96 and a published low-dimensional
record. A row for dimension 100 reading 12 886 999 232 + 24 would be arithmetic on this
paper's own result, not a result. What is *not* obvious is whether something overtakes that
implied floor the way Edel–Rains–Sloane overtakes Γ₇₂ **at** 96; the chain, evaluated through
n = 128 from the same tables used everywhere here, exceeds it by at most **0.016%** (at
n = 114), and the layered construction has no dimension above 96 to run in at all. So 96 is
where the table stops saying anything, which is why it stops there.

[`TABLE-48-96.md`](TABLE-48-96.md) is `table48-96.py`'s output: **the extrapolation of Cohn's
table from 48 to 96**, which has 48 and 72 and nothing else in the range. It is generated,
not typed — the claims come from `RESULTS.md`, the unclaimed dimensions from
`common/published.py`'s floors, and dimension 96 from `ers96.py` — and it asserts
monotonicity as it goes.

> **RESOLVED 2026-08-23.** Dimensions 62 and 63 are now claimed FROM this
> construction, in
> [`../../improved/dim62-63-ers-chain/`](../../improved/dim62-63-ers-chain/) — and from the
> WHOLE chain, not the level-0 term this section found them beaten by: (62,15,2) and
> (63,15,2) give **71 310 732** and **138 419 844**, against the 67 108 864 and 134 217 728
> below. That is why [`ers_exposure.py`](ers_exposure.py) exists: `ers_sweep.py` asks only
> about level 0, and a claim that clears level 0 has not thereby cleared Edel-Rains-Sloane.
> The cap construction over P₄₈ could not be patched to reach either figure: its absolute
> ceiling, every one of its λ(k) classes at the maximum size 7069, is 66 075 240 and
> 70 543 480, short of level 0 alone by 1 033 624 and 63 674 248
> ([`ceiling.py`](../../improved/dim49-63-p48-caps/scripts/ceiling.py)). `audit.py` checks the
> level-0 floor against every claim as its own step 4b; the full chain is checked by
> `ers_exposure.py`.

## The three numbers

| in dimension 96 | value | what it is |
|---|---|---|
| direct-sum floor | 6 218 372 160 | Γ₇₂ ⊥ Λ₂₄ rescaled to minimum 8 — 6 218 175 600 + τ(24) |
| **this mechanism, k = 24** | **6 480 558 568** | the cap construction over Γ₇₂; +4.22% on the floor |
| Edel–Rains–Sloane, chain (96,24,6,1) | **12 886 999 232** | **a factor 1.99 above it** |
| — of which level 0 alone | 8 589 934 592 | A(96,24) ≥ 2³³, from the [96,33,24] code |

**Re-read on the day of release, 2026-08-30.** This claim is a function of two published
tables and moves when either does, so both were read again: `[96,33]` still has lower bound
24 and `[96,34]` still only 22, so 33 is still the largest dimension at that distance; and
`[24,10]` over GF(4) still has lower bound 11, so `[24,9,12]₄` is still the best grid map
with `d ≥ 12`. Neither had moved, and 12 886 999 232 is unchanged.

Even the single level-0 term beats the cap construction by 33%. There is no version of this
that is close.

## Why 96 and not 95

The ERS level-0 term is `A(n₀, ⌈n₀/4⌉)`, the size of a binary code of length n₀ and minimum
distance ⌈n₀/4⌉, and codetables.de (Grassl, read 2026-08-26) gives

| n₀ | best known code | A(n₀, ⌈n₀/4⌉) |
|---|---|---|
| 93 | [93,31,24] | 2³¹ |
| 94 | [94,31,24] | 2³¹ |
| 95 | 6 473 065 046 | 4 294 967 296 = 2³² | [95,32,24] | ×1.51 |
| **96** | **[96,33,24]** | **2³³** |

`ers96.py` evaluates the whole chain at each of those and gets 2.68e9, 2.68e9, 4.83e9 and
1.29e10. **The takeover is a single step at 96, not a trend.**

**Which way the inequality points, and where it does not.** Every input to `ers96.py` is a
construction, so everything it computes is a *lower* bound on ERS. That settles dimension 96
— a lower bound on ERS standing above Γ₇₂ proves ERS wins — and it does **not** settle
93–95, where a lower bound sitting below Γ₇₂ proves nothing at all. Saying “ERS is far below
Γ₇₂ at 93–95, so those claims are unaffected” would be arguing in the wrong direction.
What can honestly be said is the threshold, which `ers96.py` prints:

| n | claim | ERS level 0 | A(n,23,23) needed to overturn it | grid map gives | Johnson allows |
|---|---|---|---|---|---|
| 93 | 6 294 850 360 | 2³¹ | 506 210 | 65 536 | 2.7e8 |
| 94 | 6 355 739 600 | 2³¹ | 513 638 | 65 536 | 3.0e8 |
| **95** | **6 473 065 046** | **2³²** | **265 816** | **65 536** | 4.0e8 |

So **dimension 95 is exposed**: a level-1 constant-weight code only 4.06× better than the
grid map would overturn it, and for calibration ERS's own `A(64,16,16) = 30 828` beats the
same grid map (4⁷ = 16 384, since codetables' best quaternary length-16 distance-8 code is
[16,7,8]) by 1.88×. That is close enough to record as a risk rather than dismiss.

The threshold column is a fact about ONE CHAIN, not about the dimension: these are the full
chain — (93,23,4,1) at 93 and 94, (95,23,4,1) at 95 — as `ers_exposure.ranked()` computes it.
`ers96.py` now reaches the same three numbers along its own independent path, choosing the
level-2 weight rather than fixing it. Until 2026-08-26 it fixed that weight at 5 and its
constant-weight table returned 1 wherever nothing was hardcoded — it had no A(95,5,5) and no
A(95,4,4) — so it charged 206 to the lower levels instead of 536 366 and printed
506 259 / 513 684 / 265 852. Those were not a second, equally correct chain: they were this
chain costed with an incomplete table, and understating the rival OVERSTATES the threshold,
which is the direction that flatters the claim. Quote the chain with the number. What
protects the claim meanwhile is that **the value to beat throughout 73–95 is the direct sum,
not an ERS total**. ERS did publish inside this range — dimension 80, at 1 368 532 064, one of
their seven worked examples — and the direct sum 6 218 175 840 beats it by 4.5×. (An earlier
version of this line said nobody had published ERS in 73–95, which is false: 80 is in it.)
So the exposure is to a referee
with a better constant-weight code, not to the literature as it stands. It is recorded here
rather than in `improved/` because it is not a defect in that package: it is the same
missing constant-weight data that makes the dimension-96 total a lower bound too.

That the sign code is what decides it is not a new phenomenon: `A(64,16) ≥ 2²⁸`, which is
codetables' [64,28,16] and is exactly the level-0 input of ERS's own published dimension-64
total, is likewise what puts dimensions 49–63 out of reach. Section 47 of
[`../../../KNOWLEDGE.md`](../../../KNOWLEDGE.md) predicted the wall at 96 from
`A(96,24) ≥ 2³²`; the true figure is **2³³**, so the wall is twice as high as predicted.

## The same question asked everywhere: dimensions 62 and 63 are beaten

The dimension-96 work is one instance of a question that has an answer in every dimension:
**how big is A(n₀, ⌈n₀/4⌉)?** That quantity is a complete kissing configuration on its own
— put ±1/√n₀ on n₀ coordinates and 0 elsewhere, and two sign vectors at Hamming distance d
meet at 1 − 2d/n₀ ≤ 1/2 exactly when d ≥ n₀/4 — so

    tau(n)  >=  max over n_0 <= n of  A(n_0, ceil(n_0/4))

with nothing added. [`ers_sweep.py`](ers_sweep.py) evaluates it against Grassl's table in
the 43 claimed dimensions at or above 49 — the only ones it covers — and it does not come
out clean:

The column marked *then claimed* is what the repository claimed at the time; dimensions
62 and 63 are now claimed from the chain itself, at 71 310 732 and 138 419 844.

| dim | then claimed | level 0 alone | from | verdict |
|---|---|---|---|---|
| 62 | 64 217 822 | **67 108 864** = 2²⁶ | [62,26,16] | **beaten** |
| 63 | 67 365 752 | **134 217 728** = 2²⁷ | [63,27,16] | **beaten, by a factor 1.99** |
| 60 | 57 477 123 | 33 554 432 = 2²⁵ | [60,25,15] | ×1.71 of margin |
| 61 | 59 898 694 | 33 554 432 = 2²⁵ | [60,25,15] | ×1.79 |
| 68 | 361 275 480 | 268 435 456 = 2²⁸ | [64,28,16] | ×1.35 |
| 95 | 6 473 065 046 | 4 294 967 296 = 2³² | [95,32,24] | ×1.51 |

⌈62/4⌉ = ⌈63/4⌉ = 16, and both codes are in the table, so **the sign vectors of [62,26,16]
and [63,27,16] are kissing configurations in dimensions 62 and 63 already**. By the standard
this repository applies in dimensions 32–44 — where an evaluation of Edel–Rains–Sloane with
the current tables *is* the live record, and is what shows Cohn's table stale at 32, 33, 34
and 37 — those two claims are not improvements.

The faulty step is in [`common/published.py`](../../../common/published.py): the 49–63 regime
argued that Edel–Rains–Sloane "does NOT reach 52 416 000 anywhere in this range" because
"2²⁸ = A(64,16) is not available below n₀ = 64". The premise is true and the conclusion does
not follow — 2²⁷ is available at n₀ = 63 and 2²⁶ at n₀ = 62, both far above 52 416 000. That
docstring now says so. **The numeric floors have deliberately been left alone**, so that
`audit.py` still certifies dimensions 62 and 63 against the old values and the discrepancy
stays visible as something to decide rather than something already absorbed. Withdrawing the
two claims changes the headline count and touches `CITATION.cff`, both READMEs and the paper,
and that is an editorial call, not a verification one.

Dimensions 49–61 are clear by a factor 1.7 or better, and everything from 69 up is clear by
2.3 or better, so the exposure is these two plus the four listed above.

## What the cap construction does at k = 24, and why it saturates there

Dimension 72+k reads `T + P` classes of a family of pairwise-disjoint Γ₇₂ classes and pairs
each with a part of a zero-sum triple partition of a 60° code `W ⊂ ℝᵏ`, gaining
`2 Σ |Cᵢ| (|Zᵢ| − 1)`. For k ≤ 23 the binding constraint is **W**: dimension 95 uses
30 990 triples + 82 pairs, consuming 93 134 of Λ₂₃'s 93 150 directions.

At k = 24, `W` is the Leech lattice's **196 560** minimal vectors, which admit
⌊196560/3⌋ = **65 520** zero-sum triples — more than double the 32 000 disjoint Γ₇₂
classes that exist. So from k = 24 upward the class family alone decides the bound, and
dimensions 96, 97, 98, … all inherit the same gain. **k = 24 is where the ladder stops
climbing**, and it stops a factor two short.

`build_triples24.py` produces a greedy partition of 43 790 disjoint zero-sum triples (131 370
of the 196 560 directions), which is 37% more than the 32 000 that can be used;
`cap96.py` re-verifies all of them exactly — every vector of norm² 32, every triple summing
to zero with pairwise inner product −16, i.e. cos = −1/2 — and takes the first 32 000.

The surplus is not wasted. A pole `(0,w)` meets a cap `(s√t·û, √(1−t)·z)` at
`⟨z,w⟩/√3 ≤ 1/2`, i.e. `⟨z,w⟩ ≤ √3/2 = 0.866`, and distinct non-parallel Leech minimal
vectors have |cos| ≤ 1/2, so **every one of the 196 560 − 96 000 = 100 560 unused directions
is a legal pole** — the first time in this repository that the pole layer is free rather than
searched for. It is worth +100 560, which does not change the verdict.

| family | classes | triples used | poles | gain | τ(96) ≥ |
|---|---|---|---|---|---|
| `images(Aut Γ₇₂)` | 32 000 | 32 000 | 100 560 | 188 701 480 | 6 406 877 080 |
| **`greedy(spread pool)`** | 32 000 | 32 000 | 100 560 | **262 382 968** | **6 480 558 568** |

Both families now have 32 000 classes; the greedy one has 2049 lines a class against the
image family's 1473, so it wins outright. At k = 23 the greedy family had been 72 classes
short of what dimension 95 needs, which is why that dimension used to fall back on the
images; at k = 24 no such threshold applies. Both families are the ones shipped and
verified in [`../../improved/dim73-95-gamma72-caps/data/`](../../improved/dim73-95-gamma72-caps/data/).

## The ERS inputs, and which of them are ours

`ers96.py` imports the binary code table A(n,d) and the chain machinery from
[`../dim32-44-ers-audit/PIPELINE.py`](../dim32-44-ers-audit/), whose formula reproduces ERS's
own published totals at n = 64 and n = 80 to the digit. Two inputs come from outside it:

1. **the level-0 sign codes** for n₀ = 93…96, read off codetables.de — the same source and
   the same quantity ERS used, as the [64,28,16] check above shows;
2. **the level-ν support codes** `A(96,w,w)`, which no table covers at this length. Lower
   bounds here use the **q-ary grid map**: split n = w·q with q a prime power, index a w×q
   grid, and let each support take one cell per row. Two supports meet in `w − d_H` of their
   q-ary choice words, so a q-ary `[w,k,≥w/2]` code gives `A(n,w,w) ≥ qᵏ`.

   | level | grid | q-ary code | A(96,w,w) ≥ |
   |---|---|---|---|
   | w = 24 | 24 × 4 | [24,9,12]₄ (codetables.de) | 4⁹ = 262 144 |
   | w = 6 | 6 × 16 | RS[6,4,3]₁₆ | 16⁴ = 65 536 |

   The Reed–Solomon one is **built and its minimum distance verified inside `ers96.py`**, over
   GF(16) with x⁴+x+1, all 65 536 codewords enumerated. The quaternary one is cited.

So 12 886 999 232 is a **lower bound on what ERS gives**, not the optimum: no one has
searched for the best constant-weight codes at length 96, and a better `A(96,24,24)` raises
it directly (ERS's own `A(64,16,16) = 30 828` beats this grid map by 1.9× at n = 64, so the
true level-1 term is probably larger still). The verdict does not depend on any of that —
level 0 alone already settles it.

## What would have to change for dimension 96 to come back into range

The cap construction needs a factor **1.99**, and the class family is the only free
parameter: `gain = 4 Σ|Cᵢ|` over as many disjoint classes as there are triples, and at k = 24
there are 65 520 triples against 32 000 classes. Filling that out — 65 520 classes of ~2 000
lines, about 131 million of Γ₇₂'s 3 109 087 800 minimal lines — would give roughly
+536 000 000 and reach ≈ 6.75e9. **Still a factor 1.9 short.** The gap is not a compute gap;
Γ₇₂ has 6 218 175 600 minimal vectors and the cap layer can only decorate them, whereas ERS
in dimension 96 is not decorating anything — it is a 2³³-word code in its own right.

## What this handed back to dimension 95 — and it has now been collected

At k = 23 the scheme needs 31 072 parts. The `greedy(spread pool)` family had **31 000**
classes — seventy-two short — so dimension 95 fell back on the `images` family. This section
used to price those 72 classes at **+69 289 824**, on the one assumption that a continued run
yields classes no smaller than the smallest it already had (1938).

**The run was extended to 32 000 classes on 2026-08-26 and the shortfall is gone.**

| | classes | mean lines/class | τ(95) ≥ |
|---|---|---|---|
| `images(Aut Γ₇₂)` | 32 000 | 1473 | 6 402 983 248 |
| `greedy(spread pool)` | 32 000 | 2049 | **6 472 914 688** |

The realised gain is **+69 838 644**, so the forecast was 0.8% low; the continued classes ran
1932–2118 with mean 2049, dipping just below the 1938 the estimate assumed. `cap96.py` now
reports this comparison rather than the forecast.

**The lesson is that it was never about the 72 classes.** The two families differ by 28% in
mean class size, so falling one class short of 31 072 does not cost one class — it costs the
whole dimension, which drops to a family that is worse everywhere. Any dimension sitting just
above a family's length is exposed the same way; the length, not the shortfall, is what to
check. And the pool had NOT run dry when the run stopped: 14 429 398 of the 80 000 000 lines
were still alive at 32 000 classes.

## Files

```
ers_sweep.py         the level-0 floor in every claimed dimension; finds 62 and 63 beaten
build_triples24.py   greedy zero-sum triple partition of the Leech's minimal vectors
cap96.py             the cap construction at k = 24, with the triples re-verified exactly,
                     and the two families compared at dimension 95
ers96.py             Edel-Rains-Sloane at n = 93..96, with RS[6,4,3]_16 built and checked
table48-96.py        tau(n) for 48 <= n <= 96, assembled from RESULTS.md + published.py
TABLE-48-96.md       its output
data/triples24.npy   43790 verified zero-sum triples, as indices into the 196560 vectors
```

## References

See [`../../../common/published.py`](../../../common/published.py). The code tables are
M. Grassl, *Bounds on the minimum distance of linear codes*, <https://codetables.de>,
read 2026-08-26.
