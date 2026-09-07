# Dimension 25: heads in the lens of a minimal vector

**τ(25) ≥ 197 569**, against the published 197 056 (Cohn's table; Ma et al. 2025) and this
repository's earlier 197 058 (now in [`../../superseded/dim25-cap-level/`](../../superseded/dim25-cap-level/)).
An improvement of **+513** over the table.

```
python verify.py            # ~4 min, exact arithmetic in every decision; exits non-zero on failure
python fullcheck.py         # ~2 min, an independent floating-point net over all 197 569 points
```

Both read only `data/` and rebuild the 196 560 Leech minimal vectors from the Golay code
(`lib/golay.py`, `lib/leech.py`), so nothing depends on the code that produced the
configuration.

## The configuration

Work in norm-4 units: Leech minimal vectors have squared length 4, and two points are
compatible iff their inner product is at most 2. Put the Leech shell on the equator of ℝ²⁵
and remove a set `D` of 1006 of its points, the *owners*. Above each owner `u` sits a
*head* `x` with `|x|² = 3` at height +1, and its mirror image at height −1, so that

| family | points | count |
|---|---|---|
| equator | `(z, 0)`, `z` a minimal vector not in `D` | 195 554 |
| extra equator point | `p = −(2/√6) v`, **not** a lattice vector | 1 |
| caps | `(x, ±1)`, `x` a head, `|x|² = 3` | 2 × 1006 |
| poles | `(0, ±2)` | 2 |
| **total** | | **197 569** |

A head `x` conflicts with an equator point `z` iff `⟨x, z⟩ > 2`, and the analysis in the
accompanying note proves that a head at `|x|² = 3` removing exactly one owner `p` has
`⟨x, p⟩ ≥ (18+√30)/7`, while two compatible heads sharing a removal both have
`⟨x, p⟩ ≤ 10/3 < (18+√30)/7`. Hence no removal is ever shared and the count is exactly
`196 560 + H + 2 + |E|` with `H` the number of heads and `|E|` the extra points.

**Class heads.** 971 of the heads are `x = u + t₋ v` with `u` an owner, `v` a Leech vector
of norm 6 (a *lean*) with `⟨u, v⟩ = −3`, and `t₋ = (3−√3)/6`, so `x = t₊ u + t₋ w` with
`w = u + v` minimal and `⟨u, w⟩ = 1`. Heads over one lean are pairwise compatible
(`⟨x, x'⟩ = ⟨u, u'⟩ − 1 ≤ 1`), so the conflict graph is multipartite with parts of at most
552. The record uses four leans with blocks 552 / 285 / 133 / 1 and Gram matrix

```
[[ 6 -4 -3  0]
 [-4  6  0 -1]
 [-3  0  6  2]
 [ 0 -1  2  6]]
```

**Interior heads.** The remaining 35 heads are exact rationals found by a linear programme
over the admissible region of an owner; they are invisible to any enumeration of the class
family.

**The extra equator point.** The first block is a whole `Σ(v) = {z : ⟨z, v⟩ = −3}` of 552
owners, and `p = −(2/√6) v` has `⟨p, z⟩ > 2` exactly for `z ∈ Σ(v)` (because `⟨v, z⟩` is an
integer in `[−3, 3]` and `2 < √6 < 3`), so it can be adjoined to the equator once that block
is gone. Against a class head over the same lean, `⟨p, x⟩ = √2`.

## What the verifier checks

`verify.py` rebuilds the minimal vectors; checks that the 1006 owners are distinct minimal
vectors and `|x|² = 3`; classifies a head as a class head when `(x − u)/t₋` is an integral
vector of norm 6 with `⟨u, v⟩ = −3`, and requires it to be **stored** at its exact value;
then verifies every pair. For the 971 class heads every inner product with an integer vector
is `A + t₋ B` with `A, B` integers, and `A + t₋ B ≤ 2` is `6A + 3B − 12 ≤ B√3`, decided by
comparing squares — all 971 × 196 560 head–equator tests and all head–head tests run with no
floating point. The 35 rational heads use `Fraction` arithmetic after a floating-point screen
six orders of magnitude wider than the arithmetic error. The extra point is checked the same
way. `fullcheck.py` shares no code with it: it assembles all 197 569 points of ℝ²⁵ and checks
every block of pairs, then 2 × 10⁷ random pairs from the whole.

## Provenance

This is joint work in progress with H. Cohn and B. Lindow (dimensions 25–31). The template
with poles is from this repository's `dim25-cap-level` (superseded); the idea of lifting the
block `Σ(v)` above a norm-6 direction is Cohn's; the no-sharing theorem, the depth bound, the
lean, interior and plateau analyses, and the extra equator point are new here. The full
account is the note *A kissing configuration of 197 569 points in dimension 25* (internal to
the collaboration), whose `factcheck.py` and `formulas.py` re-derive every figure and formula
in it. What the same analysis forbids — a ceiling of 203 552 for any configuration with a
Leech equator, and the measured reach of the Leech shell above 60° — is `KNOWLEDGE.md` §125.
