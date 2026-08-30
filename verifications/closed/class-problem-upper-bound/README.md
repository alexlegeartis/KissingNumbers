# CLOSED (to two- and three-point methods) — the class problem

**Nothing is claimed here.** This records how far the central open computational question of
the whole cap family has been pushed, so the next attempt starts where this one stopped.

> **The class problem.** How many lines of Leech minimal vectors can be pairwise at
> |cos| ≤ 1/4? Equivalently: how large can a *class* be?
>
> Best known: **248 lines** (496 vectors). Best upper bound: **425.45 lines**.

This number is the whole game for dimensions 25–31: raising the maximum class size by one
line gives +2, +8, +16, +32, +52, +96, +168 points in dimensions 25…31 respectively.
(Per extra *vector* the figures are half of those, +1, +4, +8, +16, +26, +48, +84; a line is
an antipodal pair, so mind which unit is meant.) See
[`../../improved/dim25-cap-level/`](../../improved/dim25-cap-level/) and
[`../../improved/dim27-triple-partition/`](../../improved/dim27-triple-partition/), both of
which are +2 and +496 on *whatever* the best class turns out to be.

**Two questions were attacked — beat 248, or beat 425. Neither moved.** But both are now
sharply delimited, and the negative results are strong enough to redirect effort.

```
python certificate.py       # ~1 min: the 425 bound, twice, in exact arithmetic
```

## The literature

- **248 lines / 496 vectors is the published world record**: Ma et al., *PackingStar*,
  [arXiv:2511.13391](https://arxiv.org/abs/2511.13391) (Nov 2025). Predecessors: 240
  (Cohn–Jiao–Kumar–Torquato), 244 (Kallal–Kan–Wang, arXiv:1608.07270, quoted there as "488
  vectors"). Their 496 is **maximal in the Leech shell** — no minimal vector can be added —
  re-verified here. Their structural description: 28 × X₈ + X₂₄, i.e. 28 eight-dimensional
  "symmetric frame" blocks of 16 vectors each plus one 24-dimensional Conway–Curtis cross
  of 48.
- **425 is the published upper bound**, in the form |S| ≤ ⌊9360/11⌋ = 850 *vectors*
  (Cohn–Jiao–Kumar–Torquato, *Rigidity of spherical codes*, Geom. Topol. 15 (2011), §7).
  Nothing better exists in print. (SPLAG Table 9.2 is an LP table whose entries stop at
  n = 10 and is irrelevant here.)
- Exact degree-4 certificate for 425:
  `f(t) = (4992/11) t^2 (t^2 - 1/16) = P_0 + (5083/77) P_2^24 + (27600/77) P_4^24`,
  with f ≤ 0 on [−1/4, 1/4], giving N ≤ f(1)/f₀ = **4680/11 = 425.4545**. It uses no Leech
  assumption: **425 bounds any line system in ℝ²⁴ with |cos| ≤ 1/4.**

## The 425 bound, verified exactly, by two independent routes

`certificate.py` checks both in `Fraction` arithmetic, with no floating point deciding
anything.

**Route 1 — the Hoffman ratio bound.** The 98280 lines carry a 4-class association scheme
(valencies 1/4600/47104/46575). The conflict relation is |⟨u,u′⟩| = 16, of valency 4600, and
its eigenvalues are the eigenvalues of the 4×4 integer **intersection matrix**

```
B_1 = [[   0, 4600,    0,    0],
       [   1,  892, 2816,  891],
       [   0,  275, 2300, 2025],
       [   0,   88, 2048, 2464]]
```

whose characteristic polynomial factors exactly as (x−4600)(x−1000)(x−76)(x+20) over ℤ —
checked by Faddeev–LeVerrier in exact rationals. Hoffman then gives

```
alpha <= N * (-lambda_min) / (k - lambda_min) = 98280 * 20 / 4620 = 4680/11 = 425.4545...
```

**Route 2 — the degree-4 polynomial certificate, which never mentions the Leech lattice.**
f(t) = (4992/11)·t²·(t² − 1/16) expands in the Gegenbauer basis for S²³ as

```
f = 1*P_0 + (5083/77)*P_2 + (27600/77)*P_4
```

— all coefficients non-negative, odd ones vanishing, f_0 = 1 — and f ≤ 0 wherever
t² ≤ 1/16, vanishing exactly at t = 0 and t = ±1/4. Delsarte gives N ≤ f(1)/f_0 = 4680/11.
The recurrence, the expansion, the non-negativity, the sign of f and the final quotient are
all verified exactly. **This bounds any line system in ℝ²⁴ with |cos| ≤ 1/4**, Leech or not.

> ⚠ **Do not quote `scheme_lp.py`'s output.** It builds the same scheme correctly — its
> valencies and intersection numbers are right, and it verifies them constant over random
> base pairs — but its eigenvalue extraction is numerically unstable at this size and
> returns garbage: eigenvalues near 8801.87, −0.0468 and 19489.5 where the true ones are
> 4600, 1000, 76, −20, and a bound of 6.22 where the true one is 425.45. It is kept because
> the scheme construction in it is sound and reusable; `certificate.py` is what to run.

## Every two-point method gives exactly 4680/11

On the same 4-class scheme (eigenvalues 4600/1000/76/−20, multiplicities
1/299/17250/80730):

| method | value |
|---|---|
| Delsarte LP = Schrijver θ′ | 425.4545 |
| Lovász θ | 425.4545 |
| Hoffman ratio bound | 425.4545 |

All three coincide, so **no two-point relaxation can do better**. The LP optimum has inner
distribution a₂ = 267.64 (lines at cos 1/4) and a₃ = 156.82 (orthogonal), whereas the real
248-line class has a₂ = 151.2, a₃ = 95.8 — so the LP optimum is nowhere near a real class.

## The three-point SDP does not improve it either

`bv3pt.py`. The Bachoc–Vallentin three-point bound in the *prescribed inner-product* form
(the variant Ma et al. use in their supplement D.2): antipodal codes in S²³ whose
distinct-point inner products lie in {−1, −1/4, 0, 1/4}. Primal feasibility SDP in the triple
distribution F_{u,v,t}, bisected on M = 2N. After reduction by S₃ (on the three edge labels)
and the antipodal sign action there are only **49 realisable label triples in 8 orbits**, of
which 5 are free.

At truncation d_max = 4, 6, 8, 10, 12, with and without the subconstituent LP constraints and
the doubly-stabilised ("4-point") constraints:

```
M <= 850.9091 = 9360/11    ->    N <= 425.4545     -- EXACTLY the two-point LP value.
```

An ablation shows the Y_k positivity constraints are **never active**: the bound is set by
the two-point LP alone. The three-point relaxation simply has too few degrees of freedom,
because the inner-product set has only three values.

**The implementation is validated, not broken.** The same script on the inner-product set
{−1, −1/5, 1/5} (equiangular lines at arccos 1/5) returns **276** in dimensions 24, 25 and
30, where the two-point LP gives 576 — i.e. it reproduces the hard Lemmens–Seidel / Barg–Yu
value N_{1/5}(d) = 276. At n = 16 with {0, ±1/4} it returns 144 lines (Nordstrom–Robinson,
sharp).

> **Conclusion: closing the 248–425 gap needs a four-point / Lasserre-3 relaxation, or a
> lattice-structural argument. A three-point bound cannot do it.**

## Why the 248-line class resists search from below

It is **160-rigid**: no line outside a maximum class conflicts with fewer than 4 of its
members, so any local move that removes *k* lines can add at most … nothing useful. Local
search therefore stalls immediately. The one structural handle found: every known 496-vector
class consists of 124 orthogonal **pairs** whose differences span a 4-dimensional totally
isotropic subspace R of Λ/2Λ, so the whole class lies inside

```
X_R = { line l : b(l, R) = 0 },     b(x,y) = x . y  mod 2
```

which is just **6120** of the 98 280 lines — a 16-fold reduction. Local search restricted to
such an *arena* returns classes of 230–244 lines; the identical search on the full graph
returns 152. That is what made the 644-class partition of
[`../../improved/dim38-leech-large-codimension/`](../../improved/dim38-leech-large-codimension/)
findable — but it has not yet produced a class of 249.

## Subconstituent data (new, and negative)

Fix one line *u*. The other lines split into:

- **N₃(u)**: the 46 575 lines orthogonal to *u* (the lines of the shorter Leech O₂₃).
  Relations by |ip| ∈ {16, 8, 0} form a **3-class association scheme**, valencies
  1/2464/22528/21582; its eigenvalues involve √23. Delsarte LP (forbidding |ip| = 16):
  **352.11**.
- **N₂(u)**: the 47 104 lines at |ip| = 8, signed so that u·v = +8. Relations by the *exact*
  inner product {16, 8, 0, −8, −16} form a **5-class scheme**, valencies
  1/2025/15400/22275/7128/275, multiplicities 1/23/275/2277/12650/31878 (all integers).
  Delsarte LP (forbidding ±16): **267.64**.

So N ≤ 1 + 267.64 + 352.11 = 620.7: "fix one line and split" is far worse than 425.
**Splitting loses too much.**

## History

| when | lines | vectors | who |
|---|---|---|---|
| 2011 | 240 | 480 | Cohn–Jiao–Kumar–Torquato |
| 2016 | 244 | 488 | Kallal–Kan–Wang [KKW16] |
| Nov 2025 | **248** | **496** | Ma et al. [Ma25], "PackingStar", game-theoretic RL |
| 2011 (upper) | **425.45** | 850.91 | CJKT, *Rigidity of spherical codes* — **still the best known** |
| this work | no change | no change | but the 425 bound is shown to survive every two-point method, the three-point SDP, and the subconstituent split |

Apparently the Delsarte LP value 4680/11 had not been recorded for this specific problem
before; if so, that is a small contribution of this package.

## Files

```
README.md         this file
certificate.py    THE ONE TO RUN: the 425 bound, exactly, by both routes
scheme_lp.py      the 4-class scheme built from scratch; its LP/theta/ratio output is
                  numerically unstable and should not be quoted -- see the warning above
bv3pt.py          the Bachoc-Vallentin three-point SDP, with its validation cases
octlp.py          the subconstituent schemes and their Delsarte bounds
data/octads.npy   the 759 Golay octads, used by octlp.py
data/meet4.npy    the octad meet-4 relation, precomputed, also for octlp.py
```

`certificate.py` needs only numpy; it builds the Leech lines from the Golay code via
[`common/leech.py`](../../../common/leech.py). `bv3pt.py` needs `cvxpy` and `octlp.py`
needs OR-Tools; neither is required for the bound itself.

## References

See [`../../../common/published.py`](../../../common/published.py).
