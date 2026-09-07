# SUPERSEDED — dimension 25 by the lens-head configuration

**τ(25) ≥ 197 058**, against the published 197 056 — a +2 improvement, held from 2026-08 until
2026-09-07.

**Superseded by this project's own
[`../../improved/dim25-lens-heads/`](../../improved/dim25-lens-heads/), which gives
197 569** — the same template (Leech equator, mirrored caps, two poles) with the caps
rebuilt: 1006 heads of squared length 3 in the lens of a minimal vector, each removing
exactly its owner, plus one non-lattice equator point. Better by +511. The level argument
below is still true; it is just no longer what carries the dimension-25 claim.

## What it was

**τ(25) ≥ 197058**, against the published 197056. An improvement of **+2**.

```
python verify.py            # ~30 s, needs numpy, no floating point in any decision
```

## The idea

Every record configuration in dimensions 25–31 has one shape. Write ℝ^(24+k) = ℝ²⁴ ⊕ ℝᵏ and
take

| family | points | count |
|---|---|---|
| equator | (y, 0), y a Leech minimal vector not used as a head | 196560 − 2c |
| caps | (√t·û, √(1−t)·ŵ), ŵ a unit vector in ℝᵏ | as many as the heads allow |
| poles | (0, â) | ≤ τ(k) |

The whole family is forced by the 60° condition, *except for one number*: the cap level
`t = α²`. The literature always takes **t = 2/3**, and for k ≥ 2 that is right — it is the
unique level at which a head class can be shared by a **zero-sum triple** of cap directions,
three unit vectors pairwise at cos = −1/2, which is worth more than an antipodal pair.

**In dimension 25 the second block is ℝ¹.** There are only two directions, +1 and −1. No
triple exists, so nothing forces t down to 2/3; a head only ever has to be shared by an
antipodal pair, and that is admissible for every t ≤ 3/4. Raising t to 3/4 costs nothing
and buys the two poles:

| | t = 2/3 (published) | t = 3/4 (here) |
|---|---|---|
| head condition | cos ≤ 1/4 | cos ≤ 1/3 — *the same condition on the Leech*, class still 496 |
| shared-head cap·cap | 1/3 | 1/2, still admissible |
| **cap·pole = √(1−t)** | **0.5774 > 1/2, inadmissible** | **1/2, admissible** |
| poles usable | no | **yes** |
| total | 196064 + 992 = 197056 | 196064 + 992 + 2 = **197058** |

The value t = 3/4 is *forced*, not chosen: `t ≥ 2/3` keeps the 496-vector head class,
`t ≤ 3/4` keeps antipodal head sharing, and `t ≥ 3/4` is exactly what admits the poles. The
three constraints meet in a single point.

That the Leech lattice makes "cos ≤ 1/4" and "cos ≤ 1/3" the same condition is not a
coincidence being exploited loosely — Leech inner products at norm 32 are
{0, ±8, ±16, ±32}, i.e. cos ∈ {0, ±1/4, ±1/2, ±1}, so there is simply nothing strictly
between 1/4 and 1/3.

### Why it does not extend to dimensions 26–31

There k ≥ 2, triples exist and are worth more than pairs (a head class serves 3 directions
instead of 2), so t = 2/3 is optimal; and at t = 2/3 the pole layer already attains its
ceiling τ(k) = 6, 12, 24, 40, 72, 126. Dimension 25 is the unique place where the
pair/triple trade-off is vacuous and the level can be raised for free. This was checked
exhaustively — see [`../dim73-95-gamma72-caps/`](../dim73-95-gamma72-caps/), whose model
reproduces dimensions 26, 28, 29, 30, 31 *exactly*.

## History of the lower bound in dimension 25

| when | class size s | τ(25) = 196560 + s | who |
|---|---|---|---|
| 1971 | — | 196560 | Leech–Sloane [LS71], inherited from Λ₂₄ by monotonicity |
| 2011 | 480 vectors (240 lines) | 197040 | Cohn–Jiao–Kumar–Torquato |
| 2016 | 488 vectors (244 lines) | 197048 | Kallal–Kan–Wang [KKW16], arXiv:1608.07270 |
| Nov 2025 | **496 vectors (248 lines)** | **197056** | Ma et al. [Ma25], arXiv:2511.13391 ("PackingStar"), game-theoretic RL — the current table entry |
| this work | 496, unchanged | **197058** | the same configuration at t = 3/4, which admits the two poles |

Only the last row and the current table entry 197056 are quoted directly; the middle two
are the class sizes recorded in the literature, put through the master formula
τ(25) = 196560 + s, which every configuration in this family satisfies.

The improvement is **not** a new spherical code. It reuses PackingStar's 496-vector class
exactly as published and only raises the cap level, so it is +2 on whatever the best class
happens to be. If someone finds a class of 500 vectors, dimension 25 becomes
196560 − 500 + 1000 + 2 = 197062 by the same argument.

The class size 496 is a *search* record, not a proven optimum. The best upper bound known
on it is 425 lines = 850 vectors, from the Delsarte LP — see
[`../../closed/class-problem-upper-bound/`](../../closed/class-problem-upper-bound/).

## Files

```
verify.py    the whole verification, self-contained apart from numpy
README.md    this file
```

`verify.py` does four things:

1. **builds the Leech lattice** from the Golay code (norm 32 scaling: 1104 vectors of shape
   (±4², 0²²), 97152 of shape (±2⁸, 0¹⁶) on octads, 98304 of shape (±3, ±1²³)) and checks
   there are exactly 196560 of them with the right inner-product spectrum;
2. **verifies the 496-vector head class** — every pair at cos ∈ {0, ±1/4}, none at ±1/2 —
   which is the condition the caps need;
3. **checks all six kinds of pair** (equator·equator, equator·cap, cap·cap same direction,
   cap·cap opposite directions, cap·pole, pole·pole) as exact rational inequalities at
   t = 3/4, so the configuration is a 60° code;
4. **spot-checks in real coordinates**: it materialises 3994 of the points as actual
   25-dimensional unit vectors and confirms numerically that no pair exceeds 1/2.

Steps 1–3 decide everything and use no floating point. Step 4 is a redundant sanity check
on the algebra of steps 1–3.

## References

`[LS71]`, `[KKW16]`, `[Ma25]` are in [`../../../common/published.py`](../../../common/published.py).
