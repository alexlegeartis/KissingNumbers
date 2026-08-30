# CLOSED — dimension 37: the cap construction falls 4% short, measured end to end

**Nothing is claimed here.** Dimension 37 is the natural next target of the Leech cap
construction after dimension 38, and it was run end to end. It loses, by about 4%, and the
measurement says why — so nobody repeats it.

```
N = 196560 + 4*67122 + 2*4568 + 1152 = 475336        against the record 496232  (-20896)
```

## Why dimension 37 looked promising

Dimension 38 works because at codimension k = 14 the binding constraint flips: there are up
to 644 zero-sum triples of cap directions, far more than the maximum classes could fill, so
the question becomes *how much of the Leech shell can be partitioned into classes at all* —
and the answer turned out to be all of it. See
[`../../improved/dim38-leech-large-codimension/`](../../improved/dim38-leech-large-codimension/).

Dimension 37 is k = 13, one step down, with the 1154-point dimension-13 code as the
direction set. Same regime, so the same thing should work.

## What was measured

**Triples.** The 1154-point dimension-13 code has 3800 zero-sum triples; every direction lies
in at least one and the hypergraph is connected, so nothing *structural* blocks a perfect
packing. But the **LP relaxation of the set packing gives T ≤ 363.59**, hence T ≤ 363 — not
the naive ⌊1154/3⌋ = 384. Iterated local search reaches **T = 357**; HiGHS reaches 351 with a
3.6% gap after 25 minutes. Leftover directions pair up antipodally (the code is antipodally
symmetric, 577 antipodal pairs): **P = 41**.

**Classes.** Best disjoint-class cover found — arena pool of 4096 sign changes × 121 Golay
automorphisms, with strong per-arena search — is **400 classes, 71 894 lines, mean 179.7**.
The profile decays: the first 40 classes average 238, the first 120 average 227, the first
357 average 188. A simultaneous "maximum k-colourable subgraph" local search (insertions plus
(1,1) swaps that free a line and re-place it) adds only about **1.5%**.

**Result.** 475 336 against the record 496 232, a shortfall of 20 896.

> To win, the top 357 classes would have to average **209 lines instead of 188** — an 11%
> improvement in a packing problem where the maximum single class is 248 and the greedy
> profile is already within 4% of the best simultaneous optimisation found. **Not reachable
> by search.** It would need a better *class*, not a better cover — see
> [`../class-problem-upper-bound/`](../class-problem-upper-bound/).

## The record in dimension 37 is 496 232, not Cohn's 494 312

Worth flagging on its own. Brouwer's constant-weight table carries A(37,8,8) ≥ 2832 from
Echols (arXiv:2608.13906, 2026), which through Edel–Rains–Sloane gives **496 232**. Cohn's
table still lists 494 312.

**Cohn's table is stale in four dimensions for this reason** — 32, 33, 34 and 37 — all from
Echols' 2026 constant-weight improvements, which Brouwer's page already carries:

| dim | Cohn's table | Brouwer / correct | who |
|---|---|---|---|
| 32 | 345 408 | **346 432** | Echols 2026 |
| 33 | 360 640 | **362 048** | Echols 2026 |
| 34 | 380 868 | **381 124** | Echols 2026 |
| 37 | 494 312 | **496 232** | Echols 2026 |

None of these is a claim of this project; they are corrections to the table, and they are
worth passing on. See [`../dim32-44-ers-audit/`](../dim32-44-ers-audit/) for the full audit.

## Where the leftover material lives

The dimension-37 working directory contains the arena machinery, the cover searches, the
axis rotations and the triple packings — several hundred megabytes of intermediate `.npy` and
`.pkl` files. None of it supports a claim, so none of it is shipped. What is worth keeping is
the *structural fact* it produced, which **is** shipped and **is** used: a maximum class lies
inside X_R for a 4-dimensional totally isotropic R ⊂ Λ/2Λ, just 6120 of the 98 280 lines, and
searching inside such an arena returns 230–244 line classes where the identical search on the
full graph returns 152. That is what made dimension 38's 644-class partition findable.

## Files

```
README.md   this file
```

No scripts are shipped: the result is negative and its supporting data is large. The arena
idea it produced is implemented in
[`../../improved/dim38-leech-large-codimension/`](../../improved/dim38-leech-large-codimension/)
and described in [`../class-problem-upper-bound/`](../class-problem-upper-bound/).

## References

See [`../../../common/published.py`](../../../common/published.py).
