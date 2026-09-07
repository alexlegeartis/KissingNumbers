# RECOVERED — dimensions 46 and 47, cross-sections of the extremal 48-dimensional lattice

**Nothing here is claimed.** Both values were in the literature before this project reached
them, and neither is counted among the 48 in [`../../../RESULTS.md`](../../../RESULTS.md)
or presented as a result in the paper.

| dim | Cohn's table | value | who has it |
|---|---|---|---|
| 47 | 9 741 412 | 23 766 960 | Boyvalenkov–Cherkashin, *Results in Mathematics* **80** (2025), Paper No. 3, equation (3) — [arXiv:2312.05121](https://arxiv.org/abs/2312.05121), December 2023. The same paragraph observes that those A₀ vectors form a 47-dimensional kissing configuration |
| 46 | 5 318 060 | 12 309 600 | Sun–Wang, [arXiv:2607.20359v3](https://arxiv.org/abs/2607.20359) (18 Aug 2026), Table 3, by direct enumeration of all 52 416 000 minimal vectors. The ingredients were in Ozeki's 2016 degree-3 Siegel theta tables a decade earlier, though not stated as a kissing number |

```
python derive.py            # a few seconds, needs numpy + scipy
```

**Why the package is kept.** Both numbers were obtained here independently, by a method
that uses no coordinates, no enumeration, and no knowledge of either paper — only the
design property of the minimal shell. Reproducing a published count from a route that
cannot have been fitted to it is the sharpest check the cross-section machinery has, and it
is the only place where that machinery can be tested against a published count in the
dimensions it was written for. Everything else the method produces (dimensions 70 and 71)
is in a range no table covers.

What remains true, and is a fact about the table rather than about this project, is that
**Cohn's table carries neither value**, and dimension 47's has been in the literature since
December 2023.

## The idea

Let *L* be an extremal even unimodular lattice of dimension 48 — *P₄₈p*, *P₄₈q*
(Conway–Sloane, *SPLAG* ch. 7), *P₄₈m* or *P₄₈n*, the fourth of them Nebe 2014. Three
classical facts:

- **(C1)** *L* has minimum μ = 6 and exactly *N* = 52 416 000 minimal vectors, forced by
  the extremal theta series (recomputed in `derive.py` from E₄ and Δ).
- **(C2)** [Venkov] the minimal vectors, rescaled to the unit sphere, form a spherical
  **11-design**.
- **(C3)** distinct minimal vectors have ⟨v,w⟩ ∈ {0, ±1, ±2, ±3}, and ±6 exactly when
  w = ±v — because |v−w|² = 12 − 2⟨v,w⟩ is an even integer ≥ 6.

**Lemma (a cross-section is a kissing configuration).** Fix minimal vectors u₁…u_k and let
S = {v minimal : ⟨u_i,v⟩ = 0 ∀i}. Then S lies in a (48−k)-dimensional subspace, all its
members have norm 6, and for distinct v,w ∈ S, v−w is a non-zero lattice vector so
|v−w|² = 12 − 2⟨v,w⟩ ≥ 6, i.e. cos ≤ 1/2. Hence **τ(48−k) ≥ |S|**.

### k = 1: dimension 47 is a closed-form calculation

Fact (C2) fixes every even moment Σ_v ⟨û,v̂⟩^{2p} = N(2p−1)!!/(n(n+2)…(n+2p−2)) for
2p ≤ 11. Fact (C3) says there are four unknowns. Six equations, four unknowns — the system
is **over-determined**, solves in non-negative integers, and the surplus equations hold:

| ⟨u,v⟩ | 1 | 2 | 3 | **0** |
|---|---|---|---|---|
| count (each sign) | 12 608 784 | 1 678 887 | 36 848 | **23 766 960** |

### k = 2: dimension 46 needs genuine two-point information

Everyone applies Venkov's theorem to powers of a *single* inner product, which gives only
the one-point distribution. But an 11-design integrates **every** polynomial of degree ≤ 11
correctly, including products across two fixed vectors:

> Σ_v x₁^{p₁} x₂^{p₂} = N · (sphere moment), x_i = ⟨u_i,v⟩/6, p₁+p₂ ≤ 11 and even,

the sphere moment evaluated exactly from the normalised Gram of (u₁,u₂) by Wick's theorem.
Writing n[c₁,c₂] for the joint distribution, these moments — **including the mixed ones**,
see the warning below — plus the exact one-point marginals, antipodal symmetry, lattice
translations and lattice-range constraints give a linear program whose minimum over
n[0,0] bounds the 46-dimensional cross-section.

Two ingredients take this past the classical k = 1 case:

**Translation constraints.** If ⟨u_i,v⟩ = 3 then |u_i − v|² = 6, so u_i − v is *again a
minimal vector*. Since v ↦ u_i − v is an involution, n[c] = n[c′] with c′_i = 3 and
c′_j = G_ij − c_j. This is lattice structure that the moments cannot see, and it is what
makes dimension 46 come out **exact**: on the Leech it turns the k = 3 sixty-degree-triple
bound from 26946 into the true 27720, and here it raises dimension 46 from 12 086 040 to
12 309 600.

**Lattice range.** For every lattice vector w = Σa_i u_i in the span, |v−w|² is a lattice
norm, so |⟨v,w⟩| ≤ |w|²/2 unless v = ±w. Strictly stronger than positive semidefiniteness:
for a 60° pair, u₁−u₂ is itself minimal, forcing |c₁−c₂| ≤ 3 and killing ten cells that PSD
permits. (Those ten are exactly the zeros of Ozeki's Lemma 4.2, independently rederived.)

### Realizability

A bound for a prescribed Gram is vacuous unless two minimal vectors with that Gram exist.
Here Step 1 proves it: a 60° pair needs n[3] > 0, and n[3] = 36 848. A sweep over all Gram
matrices throws up larger apparent bounds whose witnessing cell has LP minimum 0; those are
**not** used, and `derive.py` claims nothing it has not witnessed.

### Why "exact" here, and what it does not mean

For k = 1 the linear system has a unique solution, so that cross-section contains
**exactly** 23 766 960 minimal vectors. For k = 2 with a 60° pair the LP minimum and
maximum **coincide** at 12 309 600, so that cross-section has exactly that size — 31 of the
system's 47 cells are pinned to single values, and the 16 exceptions all have both
coordinates in {±1,±2}. `derive.py` computes both directions and checks they agree.

What is exact is the **size of the configuration**, not τ of the dimension. The consequence
for the kissing number is an inequality, τ(48−k) ≥ |S|. τ(46) and τ(47) remain unknown;
Cohn's table gives upper bounds 441 900 184 and 621 658 419, so there is a great deal of
room above these constructions.

### The bound cannot be corrupted by floating point

The floating-point LP is used only to *find* a dual vector y. y is rounded to rationals,
repaired into exact dual feasibility by lowering the coordinate of the all-ones moment row,
and the weak-duality bound bᵀy − Σ u_j max(0, −s_j) is evaluated in exact rationals and
rounded to an integer in the safe direction. Any exactly dual-feasible y gives a valid
bound, so an imprecise solver can only produce a *weaker* result, never an invalid one.

### ⚠ A bug worth knowing about

Earlier revisions of the LP imposed only those moments in which **every** exponent is even —
56 equations for k = 3 instead of the correct 161. The symptom was visible in the output and
went unnoticed for a day: the Gram matrices (3,−3,0) and 3·A₃ are equivalent
(UᵀGU = 3A₃ for an integer U of determinant −1) yet gave 6 462 480 and 5 194 969. **Any
future version of this machinery should be checked on a pair of equivalent Gram matrices.**
Dimensions 46 and 47 were unaffected — k = 1 has no mixed moments, and dimension 46 was
already exact — but dimensions 44, 45 and 70 all moved. The fix is in
[`../../../common/kpoint_lp.py`](../../../common/kpoint_lp.py).

## History of the lower bounds in dimensions 46 and 47

**Dimension 47**

| when | value | who |
|---|---|---|
| 1998 | 9 741 412 | Edel–Rains–Sloane [EdRS98] evaluated with Brouwer's constant-weight code table — **still the entry in Cohn's table** |
| Dec 2023 | **23 766 960** | Boyvalenkov–Cherkashin, arXiv:2312.05121, equation (3): the one-point distribution of P₄₈'s minimal shell, with the observation *in the same paragraph* that its A₀ = 23 766 960 vectors "define a 47-dimensional kissing configuration" |
| Aug 2026 | 23 766 960 | Sun–Wang, arXiv:2607.20359v3, Table 3 (complement 3A₁), by sweeping all 52 416 000 minimal vectors directly |
| this work | 23 766 960 | independent rederivation from Venkov's theorem alone |

Three independent routes, one number. **The table has simply not absorbed it.**

**Dimension 46**

| when | value | who |
|---|---|---|
| 1998 | 5 318 060 | Edel–Rains–Sloane [EdRS98] with Brouwer's codes — **still the entry in Cohn's table** |
| 2016 | (implicit) | Ozeki [Oze16], Tsukuba J. Math. 40, Table 3 tabulates the degree-3 Siegel Fourier coefficients from which 12 309 600 follows as a quotient — but does not state it as a kissing number |
| Aug 2026 | **12 309 600** | Sun–Wang, arXiv:2607.20359v3, Table 3 (complement 3A₂), by direct computation |
| this work | **12 309 600** | independently, by the two-point moment LP |

Sun–Wang's paper is dated 18 August 2026 and this work reached the same number
independently a few days later, with Ozeki's 2016 table containing the ingredients a decade
earlier. Priority is theirs and the value is **not claimed here**; what this package adds is
that the same number falls out of a method that uses no coordinates at all.

**Independent confirmation, checked in `derive.py`.** Ozeki's Table 3 gives

```
a((3,3,3,3,0,0)) = 23775066324172800000     a((3,3,3)) = 1931424768000
23775066324172800000 / 1931424768000 = 12309600      exactly, no remainder
```

That is a modular-forms computation with no linear programming in it at all. Two further
entries of the same table corroborate the Gram sweep: a pair at inner product 2 gives
11 306 800 and at inner product 1 gives 10 811 360, both below 12 309 600, so the 60° pair
really is the best rank-2 cross-section; and the **orthogonal** pair's row divided by
a((3,3,0)) is 458386320/43 = 10 660 146.98…, *not an integer* — that is an exact average
over ordered orthogonal pairs, so the count genuinely varies from pair to pair and the LP
interval [10 659 120, 10 697 520] cannot be collapsed. That is geometry, not a defect of
the method.

## What went wrong next door, and why 44 and 45 are not here

The same machinery at k = 3 and k = 4 gives dimensions 45 and 44. Those claims were
**withdrawn**: Sun–Wang compute the twelve Conway–Sloane cross-sections of P₄₈ directly and
beat them. Every one of their values lands *inside* the LP brackets computed here, which is
a sharp independent check on the machinery — see
[`../../superseded/dim44-45-p48-cross-sections/`](../../superseded/dim44-45-p48-cross-sections/).

## Files

```
derive.py    the whole derivation and every check
README.md    this file
```

The mathematics is written out as numbered lemmas, with proofs, in
[`common/PROOF-kpoint.md`](../../../common/PROOF-kpoint.md) — one document, because the
same eight lemmas carry dimensions 70 and 71 as well.

`derive.py` imports three shared modules:

| module | what it supplies |
|---|---|
| [`common/theta.py`](../../../common/theta.py) | extremal theta series; the Venkov one-point distribution, with its surplus check |
| [`common/kpoint_lp.py`](../../../common/kpoint_lp.py) | the k-point moment LP with exact rational dual certificates |
| [`common/published.py`](../../../common/published.py) | Cohn's table, so the "previously" column is never typed by hand |

`derive.py` **validates before it claims**: it first runs the identical code on E₈ → E₇
(126), Leech → Λ₂₃ (93 150), and at the two-point level on the Leech's orthogonal (43 164),
cos 1/4 (44 550) and 60° (49 896 = Λ₂₂) pairs, all of which are known independently, and all
of which come out exact.

## References

See [`../../../common/published.py`](../../../common/published.py) for the full list.
