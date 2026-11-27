# Dimension 27: obstruction to real chord lifts after the fixed Leech deletion

This is a scoped theorem about the published **210,108-point configuration of
R. Takhanov** ([source and exact certificate][takhanov]), not a new lower
bound or a maximality theorem. Its coordinates have squared norm 128 and dot
product limit 64. Keep the retained Leech shell, the axial auxiliary frames,
and the bodies of the existing paired midpoint lifts fixed. One may
recolour each existing midpoint pair with *any* signed antipodal pair
from the six lines through the roots `(±4, ±4, 0)` and their coordinate
permutations, subject only to the configuration staying valid.

**Result.** No additional point can have a body anywhere on the chord
between two already-deleted Leech vectors `r,s` with `r·s=64` and an
auxiliary tail in `R³`. The new tail is not restricted to these six lines.
The argument covers all **2,953,318** unused chords, including **2,756,992**
whose endpoints differ in both the octad and its complement. The original
selected chords are also saturated. A construction that changes the
retained Leech stock, the midpoint bodies, the axial frames, or allows
off-line continuous rearrangements of the old tails is outside the claim.

## Why the chord is constrained

The selected contact edges form a degree-two graph on all 13,248 deleted
Leech vectors. An incumbent paired midpoint lift with body
`c=(r+p)/2` requires `b·c≤64` for any new body `b`, irrespective of
the new tail: adding the two dot inequalities with opposite old tails
cancels the new tail. If `b=λr+(1−λ)s`, this gives
`λ≤(64−s·p)/(128−s·p)`. The analogous constraints at the other end give
`1−λ≤(64−r·q)/(128−r·q)`.

Let `u=max{s·p : rp is an old edge}` and
`v=max{r·q : sq is an old edge}`. A necessary condition is
`(64−u)(64−v)≥4096`. Nonincident old midpoint edges have body products
at most 64 with both endpoints by averaging, so the two local lists capture
this body-only obstruction completely. Exact incidence counts are:

| {u,v}, sorted | Unused candidate edges |
|:---|---:|
| 0,0 | 186,726 |
| 0,32 | 41,472 |
| 0,64 | 37,024 |
| 32,32 | 2,206,720 |
| 32,64 | 468,992 |
| 64,64 | 12,384 |

Only (0,0) survives, and it forces `λ=1/2`; the new tail has squared norm
32. For all mixed chords, at least one incumbent midpoint body has dot
product at least `480/7` with `b`, for every λ. The least obstructed
profile is (0,32): the lines `32+64λ` and `96−48λ` meet at `480/7`.
If a chord was already selected, its own midpoint body has product 96
with every point of the chord, so no new point can be placed there.

## Why the surviving midpoint cannot be coloured

On each positive and negative coordinate-axis auxiliary tail `±8e_i`
the same 48 body rows occur. Let `m_i` be the maximum product of the new
body with those rows. For its tail `t` we have
`8|t_i| ≤ 64−m_i`. These three constraints alone rule out squared tail
norm 32 for 184,320 of the 186,726 remaining bodies.

For every one of the remaining 2,406, there are two incumbent midpoint
bodies that both have product 64 with the new body and lie on opposite
sides in one connected component of the old product-64 graph. In any
valid six-line recolouring, a product-64 edge forces orthogonal tail
lines. Each of the six lines has a **unique** orthogonal partner within
the six-line set. The two saturated neighbours therefore have partner
lines, which span a coordinate plane. The new real tail must be
orthogonal to both, hence axial. But all three auxiliary constraints in
this remaining bank have `64−m_i≤32`, so each tail coordinate has
magnitude at most 4, ruling out an axial tail of squared norm 32.

This separates a genuinely geometric obstruction from a failed local
colour search: even coordinated discrete recolouring cannot solve it.
The next constructive path must change a supporting frame or allow
bodies beyond individual contact chords, not simply try additional tail
angles on the same chords.

## Replay

Download the original [`d27_certificate.npz`][npz] with SHA-256
`8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94`.
Then run `OPENBLAS_NUM_THREADS=2 python verify.py d27_certificate.npz`
(Python 3 and NumPy). The verifier reconstructs the Golay/Leech shell,
checks the old endpoint and opposite-tail incidences, scans the deleted
contact graph, checks the three axial envelopes, and verifies the parity
obstruction on the remaining 2,406 bodies. The binary32 matrix
multiplications are exact integer operations in this range (absolute
partial sums below 2²⁴); no floating-point tolerance is used.

[takhanov]: https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version1
[npz]: https://raw.githubusercontent.com/k-nic/KissingNumbers/77f0860bf93a745140e6e25481d3f60773b11441/version1/certificates/data/d27_certificate.npz
