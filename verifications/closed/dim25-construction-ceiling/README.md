# CLOSED — the ceiling of the dimension-25 lens-head construction

**Nothing is claimed here.** This package records how far the dimension-25 construction of
[`../../improved/dim25-lens-heads/`](../../improved/dim25-lens-heads/) can possibly go, so
that the next attempt on dimension 25 does not spend another session inside a mechanism whose
remaining room is about forty points.

> **The construction stops at 196 560 + 971 + α + 3**, where α is an upper bound on the
> independence number of one explicit 1188-vertex graph. The pool computation recorded
> **α ≤ 89**, giving **197 623**; the shipped configuration is **197 579**, **44** below it.
> α is a CP-SAT bound on an instance that does not close, so it moves from run to run — see
> [the warning below](#the-89-is-a-receipt-not-a-constant). `verify.py` derives the ceiling
> from the bound *it* proves, not from 89.

```
python verify.py            # ~16 min, needs numpy and ortools; 15 of those are the CP-SAT solve
python verify.py 60         # the same with a 60-second budget: a weaker bound, honestly reported
```

## What is closed, and what is not

**Closed:** *this* construction. Its shape is fixed —

```
(z, 0)        z in the Leech minimal vectors, minus the point each head removes
(x, +-1)      x a head, |x|^2 = 3, in the lens of a norm-6 lean
(0, +-2)      the two poles
(p, 0)        one non-lattice equator point
```

— and every term in its count is now bounded. **Not closed:** the kissing number of dimension
25. Nothing here says 197 623 is an upper bound for K(25); it is an upper bound for what
*equator + lens heads + poles* can produce. A head of some other shape, a lean of some other
norm, or a different equator, is untouched by all of this.

## The four terms

| term | value | what kind of statement |
|---|---|---|
| the equator | **196 560** | the Leech minimal vectors, fixed by the construction |
| the block family | **971** | a **computation**: a search over the whole norm-6 shell of 16 773 120 leans |
| the supplement | **α ≤ 89** | a **computation**: a CP-SAT bound on the independence number of one explicit 1188-vertex graph |
| poles and the extra point | **3** | two poles, plus one non-lattice equator point (see `dim25-lens-heads`) |

The sum is a count and not a sum of four sets, because each head *deletes* an equator point:
the total is

```
(196560 - H) + 2H + 2 + |E|  =  196560 + H + 3,      H = 971 + (supplement),  |E| = 1.
```

### What is a theorem

- **No removal is ever shared.** A head *x* that removes an equator point *p* has
  ⟨x, p⟩ ≥ (18 + √30)/7 = 3.353889…, from one spherical 7-design; two *compatible* heads
  sharing a removal would both have ⟨x, p⟩ ≤ 10/3 = 3.333333…. The intervals do not meet, so a
  shared removal is impossible. This is what makes `196560 + H + 3` an identity rather than an
  estimate, and it is why H is the **only** quantity in the construction that can still move.
  `verify.py` does not take it on trust: it measures ⟨x, u⟩ over all 971 core heads and all
  1188 pool heads, finds every one of them at 3.366025… > 10/3, and separately checks that
  each head's removal set is a single point — its own owner and nothing else.
- **A whole block is exactly 552.** A block — the top-level minimal vectors of one lean,
  pushed towards it by t₋ = (3 − √3)/6 — projects to a code at cosine 1/5 in ℝ²³, and
  A(23, 1/5) = 552. No lean can carry more, and the largest block on disk is 552 exactly.

### What is a computation

- **971.** A full-shell search over the norm-6 leans returns 971 and nothing larger. That
  search is not re-run here; what `verify.py` verifies is that the 971-head core on disk is a
  *legal* configuration — 971 heads of squared norm 3, pairwise at ⟨x, x′⟩ ≤ 1 in exact
  arithmetic over all 470 935 pairs (106 099 of them tight at exactly 1), over 971 distinct
  owners, each removing exactly its own owner. It is carried by **four** leans, in blocks
  **552 + 285 + 133 + 1**; the three large ones sit at Gram entries −4, −3, 0 (norm-4 units,
  diagonal 6).
- **1188.** The pool of heads still admissible against that core is a **completed
  enumeration**, over all **77 350** candidate owners — not a truncation and not a sample.
  `verify.py` re-derives that screen from the core and checks that all 1188 pool owners lie
  inside it. That also identifies *which* 971-head core the pool belongs to: a second,
  unrelated 971-head core in the working files screens to a different number and is not
  compatible with this pool at all.
- **α.** The supplement is a maximum independent set in one graph: the 1188 pool heads, with
  an edge whenever ⟨x, x′⟩ > 1 or the two heads want the same owner. 108 519 edges, density
  **0.1538**. CP-SAT does not close this instance, so what comes back is an upper bound.

### The 89 is a receipt, not a constant

The instance is left open by CP-SAT at 900 s, so `BestObjectiveBound` depends on the run.
Three solves of the identical model on one machine, while other work was on it:

| budget | best independent set found | bound proved | ceiling that bound gives |
|---|---|---|---|
| 10 s | 40 | 90 | 197 624 |
| 900 s | 43 | **88** | **197 622** — the run captured in `verify.log` |
| 900 s | 43 | 91 | 197 625 |

The pool computation that this package records used the same model at 900 s and proved **89**.
So the honest statement is **α ∈ [45, 91]** from everything measured — 45 because the shipped
configuration attains it — and the ceiling is 197 623 ± 2 depending on whose solve you quote.
Nothing here rounds that away: `verify.py` prints the bound its own solve proved, builds the
ceiling from that, and prints the recorded 89 beside it for comparison. If OR-Tools is missing,
or the solve returns no bound at all, the script fails instead of printing 89 as if proved.
This repository has been here before — see the CP-SAT bound that "ranged over 192, 193, 195 and
238 across runs" in [`../dim17-layered-family/`](../dim17-layered-family/), which was
eventually replaced by an exact six-number LP dual. Closing this graph exactly, or bounding it
by a method that does not depend on a wall clock, is the obvious next piece of work.

## Where the shipped configuration sits

The shipped 1016 heads are exactly this core plus **45** heads of this pool, so the shipped
supplement is an independent set in the very graph that is bounded above — 45 against α.
`verify.py` checks that decomposition against
[`../../improved/dim25-lens-heads/data/`](../../improved/dim25-lens-heads/data/) rather than
taking it from prose, and it re-derives 197 579 = 195 544 + 2 × 1016 + 2 + 1 from those files.

Two accountings of the same 1016 heads are in circulation and both are right: **971 core + 45
pool** (this package) and **972 class heads + 44 others** (the shipped package), because one of
the 45 supplementary heads happens to have exact class form over a further lean.

**So: at most about 44 more points are available inside this mechanism, and probably far
fewer** — α is an unclosed bound, not the independence number, and the best independent set any
solve has found in this graph is 45, which is what is already shipped. Anything beyond that
needs a head outside the block-plus-supplement family: a different height, a different lean
norm, or a head that is not in the lens of a minimal vector at all.

## Files

```
README.md               this file
verify.py               THE ONE TO RUN: every statement above, from the artefacts
verify.log              the captured output of a real run (the 900 s solve that proved 88)
data/rec971_X.npy       the 971-head core, 971 x 24 float64, norm-4 units (|x|^2 = 3)
data/rec971_own.npy     its 971 owner indices into common/leech.py's 196560 minimal vectors
data/p137_poolX.npy     the complete pool of admissible heads, 1188 x 24 float64
data/p137_poolU.npy     its 1188 owner indices, all distinct and disjoint from the core's
```

`verify.py` rebuilds the Leech lattice from the Golay code through
[`common/leech.py`](../../../common/leech.py) — the same route
[`dim22-23-maximal-cross-sections`](../dim22-23-maximal-cross-sections/) and
[`class-problem-upper-bound`](../class-problem-upper-bound/) use — so the owner indices are
meaningful in any checkout without the artefacts having to carry the lattice with them. The
twelve-line exact "P + Q√3 ≤ 0" decision is duplicated from
[`../../improved/dim25-lens-heads/verify.py`](../../improved/dim25-lens-heads/verify.py)
rather than imported, so that either package can be checked on its own.

Nothing in the package pins an expected count in an assertion. The assertions are relations —
`len(owners) == len(heads)`, `equator == 196560 − H`, `total == equator + 2H + 2 + |E|`,
`shipped ≤ ceiling` — because an earlier verifier in this repository asserted `n == 1006` and
thereby rejected a better configuration than the one it was written for.

## References

See [`../../../common/published.py`](../../../common/published.py).
