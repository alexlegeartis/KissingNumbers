# SUPERSEDED — dimension 38 by Edel–Rains–Sloane with the 2026 constant-weight codes

**τ(38) ≥ 570 236**, against Cohn's table's 566 652 — a +3 584 improvement, held for about
two days.

**Superseded by this project's own
[`../../improved/dim38-leech-large-codimension/`](../../improved/dim38-leech-large-codimension/),
which gives 591 612** — a completely different construction, better by +21 376.

## What it was

Dimension 38 uses the Edel–Rains–Sloane three-level chain (n, 8, 2):

> τ_n ≥ A(n, ⌈n/4⌉) + 128·A(n,8,8) + 4·C(n,2)

W. Echols, *New lower bounds for constant-weight codes via seeded bit-swap tabu search*,
arXiv:2608.13906 (2026), improved A(38,8,8) from 2997 to **3025**. Substituting gives

```
180224 + 128*3025 + 4*703 = 180224 + 387200 + 2812 = 570236
```

a gain of exactly 128 × (3025 − 2997) = 3584.

**Why nobody had done this.** Echols derives his own kissing numbers from a fixed level-1
length n₀ = 32, which in dimension 38 gives 131072 + 128·3025 + 2812 = 521 084 — *below* the
record — so dimensions 38 and 39 were passed over. They are exactly the two dimensions where
the ERS chain uses n₀ = n rather than n₀ = 32, with the much larger level-1 codes
A(38,10) ≥ 180 224 and A(39,10) ≥ 327 680.

The same observation gives dimension 39, where it **still stands**:
[`../../improved/dim39-ers-constant-weight/`](../../improved/dim39-ers-constant-weight/).

## Where the verification lives

There is no separate script here. `improved/dim39-ers-constant-weight/scripts/verify.py`
checks both dimensions in one run — they share a single argument — and its output for
dimension 38 is this claim:

```
cd ../../improved/dim39-ers-constant-weight && python scripts/verify.py
```

It verifies exhaustively that `data/a38.8.8.3025H` contains 3025 distinct binary words of
length ≤ 38, every one of weight 8, with maximum pairwise support intersection 4, hence
minimum Hamming distance exactly 8 (all 4 573 800 pairs compared); that the level-2 sign
code has 128 words at minimum distance 2; the chain condition; and all ten geometric
conditions as exact rational inequalities.

## What is taken on trust, in both dimensions

A(38,10) ≥ 180 224 and A(39,10) ≥ 327 680, due to Zinov'ev and Litsyn, *On shortening of
codes*, Problems Inform. Transmission 20:1 (1984) 1–7, tabulated by Litsyn–Rains–Sloane. No
explicit codewords for these appear to be published anywhere.

**This is not a new assumption**: the existing records 566 652 and 755 988 in Cohn's table
rest on exactly the same two numbers with the same provenance. For reference, the best known
*linear* codes give only [38,17,10] and [39,18,10], i.e. 2¹⁷ and 2¹⁸; substituting those
yields τ(38) ≥ 521 084 and τ(39) ≥ 690 580, which are unconditional but below the current
records.

## History of the lower bound in dimension 38

| when | value | who |
|---|---|---|
| 1982 | 224 608 | Conway–Sloane, cross-section of P₄₈ |
| 1998 | 566 652 | Edel–Rains–Sloane [EdRS98] with the constant-weight codes of the day — **still the entry in Cohn's table** |
| 2026 | (skipped) | Echols [Ech26] improves A(38,8,8) but uses a fixed n₀ = 32, which misses this dimension |
| Aug 2026 | 570 236 | **this package** — the same chain with n₀ = n |
| Aug 2026 | **591 612** | this project's Leech cap construction at codimension 14 — see [`../../improved/dim38-leech-large-codimension/`](../../improved/dim38-leech-large-codimension/) |

## References

See [`../../../common/published.py`](../../../common/published.py).
