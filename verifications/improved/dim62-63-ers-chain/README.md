# Dimensions 62 and 63 — the Edel–Rains–Sloane chain

| dim | previously | this work | gain | factor | chain | verified |
|---|---|---|---|---|---|---|
| 62 | 52 417 932 | **71 310 732** | +18 892 800 | 1.36 | (62, 15, 2) | levels 1–2 built here, level-0 code cited |
| 63 | 52 418 564 | **138 419 844** | +86 001 280 | 2.64 | (63, 15, 2) | levels 1–2 built here, level-0 code cited |

```
python scripts/verify.py           # < 1 s   level 0: the sign code, and the 60-degree condition
python scripts/verify_chain.py     # ~20 s   the whole chain, and tau(8) and tau(16) reproduced
```

| dim | level 0 | level 1 | level 2 | total |
|---|---|---|---|---|
| 62 | `A(62,16) ≥ 2²⁶` = 67 108 864 | 4096 × 1024 = 4 194 304 | 1891 × 4 = 7 564 | **71 310 732** |
| 63 | `A(63,16) ≥ 2²⁷` = 134 217 728 | 4096 × 1024 = 4 194 304 | 1953 × 4 = 7 812 | **138 419 844** |

**The levels above 0 are built, not cited.** Level 1 needs 4096 supports of weight 15 meeting
pairwise in at most 7, and 1024 sign patterns on each. The supports come from a **cyclic
`[15,6,8]₄`** — generator polynomial with defining set `{0,1,2,3,4,5,8,10,12}`, a union of
cyclotomic cosets mod 15 under `x ↦ 4x` — whose minimum distance 8 is established by
enumerating all 4⁶ = 4096 codewords, pushed through the **grid map** (index the cells of a
15 × 4 grid, let each support take one cell per row; two supports then meet in exactly
15 − d_H ≤ 7). All 8 386 560 pairwise overlaps are checked directly. The signs come from a
`[15,10,4]`, built as the extended Hamming `[16,11,4]` shortened once, with its minimum
distance enumerated over all 1024 codewords. Level 2 is every pair of coordinates with all
four sign patterns. `15 · 4 = 60 ≤ 62`, and `4 · 16 = 64 > 63`, so weight 15 is the largest
the chain allows in either dimension.

**Why it was missed.** `ers_sweep.py` asked only about level 0 — `A(n₀, ⌈n₀/4⌉)`, which is a
complete kissing configuration on its own and which is what these two dimensions used to
claim. It is not the whole construction: the levels above 0 are what carry Edel–Rains–Sloane's
own dimension-64 total from 2²⁸ = 268 435 456 to 331 737 984. The full sweep is now
[`../../closed/dim96-ers-takeover/ers_exposure.py`](../../closed/dim96-ers-takeover/ers_exposure.py),
which evaluates the whole chain in every claimed dimension; these two were the only ones it
found below their own rival, and the fix is to claim the rival.

## Level 0, on its own

Let `C` be a binary linear code of length `n₀`, dimension `k`, minimum distance `d`. Send each
codeword to a sign vector

    x(c) = ( (−1)^c₁, …, (−1)^c_{n₀} ) / √n₀ ,

a unit vector of ℝ^{n₀}. Two codewords at Hamming distance `t` give

    ⟨x(c), x(c′)⟩ = (n₀ − 2t) / n₀ ,

so **every** pair is at 60° or more exactly when `t ≥ n₀/4` throughout, that is when

    d ≥ ⌈n₀/4⌉ .

The `2^k` sign vectors are then a kissing configuration of ℝ^{n₀}, and of ℝⁿ for every `n ≥ n₀`,
so

    τ(n)  ≥  max over n₀ ≤ n of  A(n₀, ⌈n₀/4⌉) .

This is level 0 of Edel–Rains–Sloane with a single support: no chain, no support codes, no glue.
`⌈62/4⌉ = ⌈63/4⌉ = 16`, and codes of length 62 and 63 with `d = 16` reach dimension 26 and 27.
The rest of the chain is in `scripts/verify_chain.py` and in the table at the top; this section
is the part that carries 94% of the total.

## Why these two dimensions are here and not in the P₄₈ package

[`../dim49-63-p48-caps/`](../dim49-63-p48-caps/) builds layers over P₄₈ and carries dimensions
49 to 61 comfortably. In 62 and 63 it does not merely fall short — it **cannot** reach these
values, with any class family whatsoever. Its bound is

    τ(48+k) ≥ 52 416 000 + 2·S(λ(k)) + τ(k),

in which `52 416 000` is τ(P₄₈), `τ(k)` is the pole layer and `λ(k)` is fixed by the kissing
record in dimension `k`. The one free quantity is `S`, the total size of the `λ(k)` disjoint
classes, and `S ≤ λ(k)·7069` because 7069 is the largest class that exists. That gives a hard
ceiling, computed by [`../dim49-63-p48-caps/scripts/ceiling.py`](../dim49-63-p48-caps/scripts/ceiling.py):

| dim | k | λ(k) | cap construction reaches | its absolute ceiling | level 0 gives | short by |
|---|---|---|---|---|---|---|
| 62 | 14 | 966 | 64 217 822 | 66 075 240 | 67 108 864 | 1 033 624 |
| 63 | 15 | 1282 | 67 365 752 | 70 543 480 | 134 217 728 | 63 674 248 |

Raising `λ(k)` would itself require a new kissing-number record in dimension 14 or 15. So the
cap construction — which is what carries twenty-eight other dimensions in this repository — is
simply the wrong construction here, and the one-line one wins.

Dimensions 60 and 61 are the closest calls that go the other way: level 0 gives 33 554 432 there
(`[60,25,15]`, `⌈60/4⌉ = 15`), against 57 477 123 and 59 898 694 from the caps.

## What is verified, and what is cited

`scripts/verify.py` checks, in exact rational arithmetic, that the 60° condition holds for
every attainable Hamming distance at `n₀ = 62, 63` — with a negative control at `d = 15`, one
below the threshold, which must and does fail. It then runs the whole mechanism **end to end**
on three codes exhibited in coordinates, enumerating every codeword and computing the minimum
distance directly:

| code | `[n₀,k,d]` | needs `d ≥` | gives |
|---|---|---|---|
| extended Golay | `[24,12,8]` | 6 | τ(24) ≥ 4096 |
| Reed–Muller RM(1,5) | `[32,6,16]` | 8 | τ(32) ≥ 64 |
| Reed–Muller RM(2,5) | `[32,16,8]` | 8 | τ(32) ≥ 65 536 |

RM(2,5) is the interesting one: `d = ⌈32/4⌉ = 8` with equality, so the boundary of the condition
is exercised rather than cleared comfortably. Every *pair* is covered as well: for a linear code
the differences of codewords are again codewords, so the largest inner product is exactly
`(n₀ − 2d)/n₀` and the minimum distance settles all pairs at once. That identity is checked
against brute force over all 4096² pairs of the Golay code, where the Gram still fits in memory.

**Taken from the literature: only the two code parameters.** `[62,26,16]` and `[63,27,16]`,
from Grassl's table of best known linear codes, [codetables.de](https://codetables.de), read
2026-08-23, **re-read 2026-08-25** and again on the day of release, 2026-08-30, when
`[62,26]` and `[63,27]` still read 16 and `[62,27]` and `[63,28]` still read 15
(lower bound 16, upper bound 18 in both cases; the
`BKLC` pages state `Lb(63,27) = 16 is found by shortening of: Lb(64,28) = 16` and, for
n = 62, `Shortening of [4] at { 63 .. 64 }` from the same [64,28,16]).  `[62,27]` and
`[63,28]` still read 15, so the two values have not doubled. Both are shortenings of a single
`[64,28,16]`, built there as

> cyclic `[73,36,16]` with generator polynomial
> `x³⁷+x³⁶+x³⁴+x³³+x³²+x²⁷+x²⁵+x²⁴+x²²+x²¹+x¹⁹+x¹⁸+x¹⁵+x¹¹+x¹⁰+x⁸+x⁷+x⁵+x³+1`,
> punctured to `[72,36,15]`, Construction B2 to `[63,28,15]`, extended to `[64,28,16]`.

Shortening once gives `[63,27,16]` and twice `[62,26,16]`; shortening an `[n,k,d]` code gives
`[n−1,k−1,≥d]`, so `d = 16` survives both times. **Exhibiting that `[64,28,16]` in coordinates
is the one thing this package does not do**, and the chain is recorded above so that it can be.
That is exactly the standing of `A(38,10)` and `A(39,10)` in
[`../dim39-ers-constant-weight/`](../dim39-ers-constant-weight/), which are cited and not
rebuilt.

**Closing that gap is a bounded piece of work, if anyone wants it.** Note that Brouwer's
construction does not have to be reproduced faithfully — *any* `[64,28,16]` will do, because the
verification is self-certifying: if the minimum distance of the code you build comes out 16, it
is a `[64,28,16]` whether or not you followed the same chain. And checking that is cheap.
Represent each of the 28 generator rows as a `uint64`, split the information vector into two
halves of 14 bits, precompute the 2¹⁴ combinations of each half, and XOR the two tables against
each other: 2²⁸ ≈ 268 million codewords, one `uint64` each, popcount and take the minimum — well
under a minute in NumPy. The obstacle is purely in obtaining the generator matrix; the step this
package currently cites is `Construction B2`, whose precise definition we did not have to hand.
Shortening the result at one and two positions then gives `[63,27,16]` and `[62,26,16]`
directly, and shortening is a two-line operation.

The **whole chain** is checked end to end in coordinates at the two dimensions where it happens
to reproduce the exact kissing number — `n = 8`, chain (8, 2), giving 240 points, which is the
E₈ root system; and `n = 16`, chain (16, 4, 1), giving 4320, which is the Barnes–Wall value and
the record. Every pair of both configurations is checked directly (28 680 and 9 329 040 of
them). Landing on 240 and 4320 exactly, from a construction that knows nothing of E₈ or of
BW₁₆, is what licenses running the same chain at 62 and 63.

Both dimensions therefore carry the status `external` in
[`../../../RESULTS.md`](../../../RESULTS.md): they improve the moment the code table does, and
need no new mathematics from this repository to do so. If `[62,27,16]` or `[63,28,16]` is ever
found, these numbers double.
