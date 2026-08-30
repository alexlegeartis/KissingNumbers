# CLOSED — dimensions 22 and 23: Λ₂₁, Λ₂₂, Λ₂₃ are maximal spherical codes

**Nothing is claimed here.** Dimensions 22 and 23 were attacked from four directions and all
four are closed, with exact numbers.

```
python verify.py            # ~2 min, needs numpy; every claim below, in exact arithmetic
```

Cohn's table gives 49 896 for dimension 22 and 93 150 for dimension 23, both from the Leech
cross-sections Λ₂₂ and Λ₂₃ (Leech–Sloane 1971). They have not moved in 55 years, and this
package explains why the obvious attacks fail.

## (a) No free point exists

The trick that gains +2 in dimension 25 — find a gap and insert a point — does not repeat
here. Computing m(C) = min over unit *z* of max over *x* ∈ C of |⟨z, x̂⟩| (a point can be
added iff m(C) ≤ 1/2), by maximising |z| over the polytope {z : |⟨z,x̂⟩| ≤ 1/2}:

| code | N | m(C) | exact | margin over 1/2 |
|---|---|---|---|---|
| Λ₂₁ | 27 720 | 0.525225731 | **√(8/29)** | +5.0% |
| Λ₂₂ | 49 896 | 0.522232968 | **√(3/11)** | +4.4% |
| Λ₂₃ | 93 150 | 0.516397779 | **√(4/15)** | +3.3% |

Each optimal direction turns out to be the *next* lamination's middle-layer direction
projected down: Λ₂₁'s deepest hole is Λ₂₂'s middle layer, Λ₂₂'s is Λ₂₃'s, and Λ₂₃'s is the
Leech's level ±1. Every value is 16/√(32·|v′|²) with |v′|² = 29, 88/3, 30, checked in exact
`Fraction` arithmetic.

> **Method warning, worth keeping.** A naive iterated LP over all constraints is unreliable
> for this. On a control case (Λ₂₂ minus its two outer layers, where the axis is *provably*
> free with max|z| = √(8/3)) it returned an answer 70% low. It has to be a cutting-plane
> ascent with an explicit feasibility assertion on every returned point.

**A corollary that explains Λ₂₂'s shape.** For a dimension-22 code split over an axis with
equator E ⊂ ℝ²¹ and a layer at cosine-height h, the layer is EMPTY unless
1/(2√(1−h²)) ≥ m(E). With E = Λ₂₁ that forces h² ≥ 3/32 — *precisely* Λ₂₂'s middle-layer
height. So the layer's apparent LP slack (10 752 used against a two-sided Delsarte ceiling of
19 853) is unreachable: the lower, roomier heights have no admissible directions at all.

## (b) Level sets of the Leech: closed exactly at 93 150

In exact rational arithmetic (norm² = 32), the level pairs (0,0), (8,8), (16,16), (8,16) have
no forbidden inner product, so levels 8 and 16 together (47 104 + 4 600 = 51 704) form a
valid 23-dimensional code — but a small one. The interesting question is whether level 0
(= Λ₂₃) can be *augmented* by part of levels 8 and 16.

It cannot. Level 0 against the positive half is bipartite with **exactly regular degrees**:

```
level 0 -> level 8  = 1024,  back 2025     (93150*1024 = 47104*2025)
level 0 -> level 16 =   44,  back  891     (93150*  44 =  4600* 891)
```

Hall's condition holds *strictly*, because 1024/2025 + 44/891 = 1124/2025 < 1, so a matching
saturates the positive half, and by König the maximum independent set is **93 150 exactly**.

> No union of level sets, and no augmentation of one, can beat Λ₂₃. This killed the most
> promising untried idea in the project.

## (c) Classes over a 60° pair (dimension 22): closed exactly at 49 896

The pair splits the Leech's 196 560 minimal vectors into 6 in-plane roots plus 19 classes —
one of size 49 896, six of 20 736, six of 2 816, six of 891. All 19 are internally clean and
105 of the 171 class pairs are compatible, but the maximum compatible union is **49 896**:
the (0,0) class alone, which is Λ₂₂.

## (d) The Cohn–Li mechanism is dead here too

Dimensions 22 and 23 leave gaps of 6 733 and 27 371 against a 4 096-word ambient dual code,
so the extra term cannot possibly close them. See
[`../dim17-23-cohn-li-mechanism/`](../dim17-23-cohn-li-mechanism/).

## History of the lower bounds in dimensions 22 and 23

| when | dim 22 | dim 23 | who |
|---|---|---|---|
| 1971 | **49 896** | **93 150** | Leech–Sloane [LS71], the Leech cross-sections Λ₂₂ and Λ₂₃ |
| — | unchanged for 55 years | unchanged | |
| this work | — | — | **no improvement, and four separate routes shown to be closed** |

## Files

```
README.md          this file
verify.py          every claim above, in exact integer / rational arithmetic
base.py            Leech minimal vectors and the cross-sections
freepoint.py       the cutting-plane ascent for m(C)
lpvertex.py        its inner loop: maximise |z| over the polytope, from random starts
levelsets.py       the Hall/Koenig argument on the level graph
classes22.py       the 19 classes over a 60-degree pair
lpbounds.py        the two-point Delsarte bounds
lam21.py           the Lambda_21 case
lpz_Lambda_*.npy   the optimal directions found, one per code
```

`base.py` imports the Leech lattice from [`common/leech.py`](../../../common/leech.py), which
builds it from the Golay code and caches it; nothing is taken on trust.

## References

See [`../../../common/published.py`](../../../common/published.py).
