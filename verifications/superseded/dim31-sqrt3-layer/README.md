# Dimension 31: a deletion-free layer at height √3 in the deep holes of E7

> **SUPERSEDED on 2026-09-17** by `../../improved/dim31-frame-layer/`, which reaches
> **τ(31) ≥ 238 662** with the norm-8 frame layer at height √2. Nothing below is wrong; the
> height-√3 layer is simply the smaller one (+4 over the published 238 350, against +312), and
> it was the record while
> the frame layer's 42 disjoint type-B classes could not be packed densely enough.

**τ(31) ≥ 238 354**, against the published 238 350 (Cohn's table; Ma et al. 2025,
arXiv:2511.13391, the PackingStar configurations). An improvement of **+4**, and — unlike
every other layer of the 25–31 programme — one that does not touch the class problem at all:
the published classes are used exactly as they are.

```
python verify31.py          # ~2 min, exact arithmetic in every decision; exits non-zero on failure
```

It reads only `data/` and rebuilds the 196 560 Leech minimal vectors from the Golay code
(`lib/golay.py`, `lib/leech.py`, the coordinates of Cohn's file), so nothing depends on the
code that produced the configuration. Every pair of points is decided exactly: integer inner
products for the ℝ²⁴ parts, sympy in ℚ(√2, √3) for the ℝ⁷ parts.

## The configuration

Norm-4 units as in the dimension-30 package: a point of ℝ³¹ is `(x, y)` with `y ∈ ℝ⁷`, two
points are compatible iff their inner product is at most 2.

| | points |
|---|---|
| equator `(u, 0)`, `u` a minimal vector that is not an owner | 175 728 |
| caps `(√(2/3) u, (2/√3) z)`, `u` an owner, `z` a direction of its triangle | 3 × 20 832 |
| axis `(0, 2a)`, a rotated copy of E7, minus the four points blocked by the layer | 126 − 4 |
| **layer** `(±u₀/2, √3 w)`, `u₀` a minimal non-owner, `w` one of four deep holes | **8** |
| **total** | **238 354** |

The 126 cap directions are the E7 roots in 42 zero-sum triangles, each carrying a 496-owner
class (the 248-line class of Ma et al. and its images); the axis copy of E7 touches the caps
at 30°. That is the published template, verified here from scratch.

**The layer.** A point at height `h` (that is `|y| = h`) is compatible with an axis point iff
`⟨y/h, a⟩ ≤ 1/h`. At `h = √3` the threshold is `1/√3`, and the covering radius of the E7 root
system on `S⁶` is *exactly* `arccos(1/√3) = 54.7356°`, attained at the 56 minimal directions
of E7* — the deep holes. Take such a hole `w` and any minimal vector `u₀` that is not an owner,
and the two points `(±u₀/2, √3 w)`:

| partner | inner product | bound |
|---|---|---|
| equator `(v, 0)` | `⟨u₀, v⟩/2 ≤ 2`, equality only at `v = u₀` | touching, no deletion (Cauchy–Schwarz) |
| cap of an owner `u` at direction `z` | `⟨u₀, u⟩/(2√6) + 2⟨z, w⟩ ≤ 2/√6 + 2/√3 = 1.9712` | free: `u₀` is not an owner, so `⟨u₀, u⟩ ≤ 2`, and a deep hole has `⟨z, w⟩ ≤ 1/√3`; exact margin `√6 + 2√3 < 6` |
| axis `(0, 2a)` | `2√3 ⟨w, a⟩ ≤ 2` iff `⟨w, a⟩ ≤ 1/√3` | the four axis points above that cosine are deleted |
| the other layer points | `±1 + 3⟨w, w′⟩ ≤ 2` | the two hole lines used are at cosine ±1/3 |

So each deep-hole line costs the axis points within `arccos(1/√3)` of it and gives four points;
two lines (`±w₁, ±w₂`, at cosine ±1/3) block only four axis points between them, for a net
**+4**. That this is the layer's ceiling is established four ways in the working notes: for two
distinct deep-hole lines `w_i ∓ w_j` is a root of the cap E7, so a full 126-point axis hosts at
most one free line (the shared-root lemma); a CP-SAT over every deletion subset is optimal at 4;
3.8 million axis rotations by annealing never beat 4; and letting the directions drift inside
their admissible blobs costs at least one axis line each.

What `verify31.py` checks, in order: the 126 cap directions form an E7 root system in 42
zero-sum unit triangles and the four layer directions are deep holes (cosine `0` or `±1/√3` to
every root); the 42 classes are disjoint sets of Leech minimal vectors, each with inner products at most 1
(at least arccos(1/4) = 75.52°, not just 60°, because a class shares all three of its cap
directions), and `u₀` is
a minimal vector outside them; cap–cap through the 42 × 42 class maxima and the integer
thresholds `24 − 16 cos`; exactly four axis points are blocked by the layer, and the surviving
122 are admissible pairwise and against every cap; equator–cap; and the layer against the
equator (`max ⟨u, u₀⟩ = 32`, touching), the caps (`max ⟨u_owner, u₀⟩ = 16`, worst inner product
`2/√6 + 2/√3`), the surviving axis, and itself. `verify31.log` is the output.

## Files

```
verify31.py          the verification (exact); prints ALL CHECKS PASS  K(31) >= 238354
verify31.log         its output
data/classes.npy     42 x 496 x 24 int8: the owner classes of the published construction, Cohn coordinates
data/u0.npy          the minimal non-owner carrying the layer (norm 32, Cohn coordinates)
data/directions.json unit vectors of R^7 as sympy strings: the 42 triangles (class i on triangle i),
                     the 126 axis directions, and the 4 layer directions (two lines of deep holes)
lib/golay.py         the Golay code in the coordinates of Cohn's file (as in superseded/dim30-sqrt2-layer)
lib/leech.py         the 196560 minimal vectors from it
```

## Where it comes from

Found on 2026-09-13 by porting the dimension-30 height-√2 layer: the port died on placement
(27 disjoint class images would be needed where dimension 30 needs 16), but the height at which
the axis threshold meets the covering radius of the root system exists in dimension 31 alone —
E6 (dimension 30) has its holes at cosine 0.6124, D5 (29) at 0.6316, the 24-cell (28) at 0.7071,
all above the cap threshold `(2 − √(2/3))/2 = 0.5918` that a class-independent layer needs.
Working notes: `research/collab2531/BRAINSTORM_31.md`, FINDINGS section 47, KNOWLEDGE.md
section 129. The only large lever left in this dimension is the class: +168 per extra line.

## History of the bound in dimension 31

| year | bound | who |
|---|---|---|
| 2025 | 238 350 | Ma et al., arXiv:2511.13391 (PackingStar), the value in Cohn's table; the cap construction over the Leech lattice with its 248-line class reproduces it exactly (`common/capalgebra.py`) |
| 2026-09-13 | **238 354** | this package |
