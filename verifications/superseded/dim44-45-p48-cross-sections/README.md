# WITHDRAWN — dimensions 44 and 45

| dim | Gram | this project certified | Sun–Wang | verdict |
|---|---|---|---|---|
| 45 | 3A₃ = 3D₃ | 6 687 320 | **6 702 080** (3D₃), and an antipode packing at **7 379 838** | theirs is better |
| 44 | 3D₄ | 4 147 967 | **4 153 600** (3D₄) | theirs is better, by 5 633 |

```
python brackets.py          # ~2 min: recomputes the brackets and checks Sun-Wang lands inside
```

Both of this project's values are valid lower bounds on τ(45) and τ(44), and both improve on
Cohn's table (3 047 160 and 2 948 552). Neither is the best known. **X. Sun and C. Wang,
*Sphere packings and kissing numbers in dimensions 39, 43 and 45 from the antipode
construction*, [arXiv:2607.20359v3](https://arxiv.org/abs/2607.20359) (18 August 2026),
Table 3**, compute the kissing numbers of the twelve Conway–Sloane (1982) cross-sections of
P₄₈ directly, by sweeping all 52 416 000 minimal vectors, and beat both.

**Dimensions 44 and 45 are therefore withdrawn from this project's claims.**

## Why this directory still exists

Because the comparison is a sharp, independent check on machinery that four surviving claims
depend on — dimensions 46, 47, 70 and 71 all come from the same k-point moment LP.

| dim | complement | Sun–Wang | this project's LP bracket | inside? |
|---|---|---|---|---|
| 47 | 3A₁ | 23 766 960 | exact, 23 766 960 | **identical** |
| 46 | 3A₂ | 12 309 600 | exact, 12 309 600 | **identical** |
| 45 | 3A₃ | 6 687 648 | [6 687 320, 6 702 736] | ✔ |
| 45 | 3D₃ | 6 702 080 | [6 687 320, 6 702 736] | ✔ |
| 44 | 3D₄ | 4 153 600 | [**4 147 967**, 4 165 632] | ✔ — a bracket 0.4% wide |

Two independent methods — linear programming on design moments here, direct enumeration of
minimal vectors there — agree exactly where the LP is exact and bracket each other where it
is not. That is about as good a cross-validation as this kind of machinery admits.

## Two things the comparison settles

**The count depends on the embedding, not just the abstract Gram matrix.** Sun–Wang's 3A₃
and 3D₃ are the same abstract lattice, yet give 6 687 648 and 6 702 080. The k-point LP sees
only the Gram matrix, so it can *never* do better than bracket them. The width of the
bracket [6 687 320, 6 702 736] is not slack in the method — it is the real spread of the
geometry, and it contains both of their values almost exactly at its two ends.

**Ozeki's Theorem 6.3 falls out of the LP.** With the mixed moments (see the warning in
[`common/PROOF-kpoint.md`](../../../common/PROOF-kpoint.md) §2), the three-point system for
a 60° triple in an extremal 48-dimensional lattice has affine solution space of dimension
exactly **1**:

```
n[0,0,0] + 164 * n[3,3,3] = 6732912       (per triple, exact)
```

Non-negativity forces 184 ≤ n[3,3,3] ≤ 278, which is exactly Ozeki's bound on
τ = a(T₄,3)/a(T₃,1) — obtained here by linear programming with an exact rational dual
certificate rather than by modular forms. The two ends of that interval are the two ends of
the dimension-45 bracket. `brackets.py` verifies all of this.


## Getting the Gram matrix right

Sun–Wang's dimension-44 complement is **3D₄**, the D₄ root lattice scaled to minimum 6, whose
Gram has off-diagonals in {0, −3}. It is *not* 3A₄, whose Gram has all off-diagonals 3. The
two are different lattices and give different cross-sections:

| Gram | LP bracket | contains Sun–Wang's 4 153 600? |
|---|---|---|
| 3A₄ | [3 694 470, 3 833 121] | no — and it should not, it is a different configuration |
| **3D₄** | [**4 147 967**, **4 165 632**] | **yes** |

This mistake was made once here, and it announces itself loudly: the "independent check"
appears to fail, with a real configuration sitting above what is supposed to be an upper
bound. If that ever happens, suspect the Gram matrix before suspecting the LP.

## History of the lower bounds in dimensions 44 and 45

| when | dim 44 | dim 45 | who |
|---|---|---|---|
| 1982 | | | Conway–Sloane cross-sections of P₄₈, *SPLAG* — the configurations, though not their kissing numbers, were known |
| 1998 | 2 948 552 | 3 047 160 | Edel–Rains–Sloane [EdRS98] — **still the entries in Cohn's table** |
| Aug 2026 (this project, later withdrawn) | 4 147 967 | 6 687 320 | k-point moment LP on P₄₈ |
| Aug 2026 | **4 153 600** | **7 379 838** | Sun–Wang [SW26] — 4 153 600 and 6 702 080 by direct enumeration of the cross-sections, 7 379 838 by the antipode construction |

The antipode value 7 379 838 lies *outside* this project's bracket, correctly — that
configuration is not a cross-section of P₄₈ at all, so the bracket says nothing about it.

## A note on the earlier, larger claim

Before the mixed-moment bug was found, this project reported **6 462 480** for dimension 45
and listed it in its results table. That number came from an LP missing 105 of its 161
moment equations; the symptom was that two equivalent Gram matrices, (3,−3,0) and 3·A₃, gave
6 462 480 and 5 194 969. With the correct moment family both give 6 687 320. The 6 462 480
figure appears in early drafts of this project and **should be disregarded** — it happens to
be a valid bound, but it was produced by code that was wrong.

## Files

```
brackets.py   recomputes the LP brackets, checks Sun-Wang lands inside each,
              and rederives Ozeki's 184 <= n[3,3,3] <= 278
README.md     this file
```

`brackets.py` uses [`common/kpoint_lp.py`](../../../common/kpoint_lp.py) and
[`common/published.py`](../../../common/published.py).
