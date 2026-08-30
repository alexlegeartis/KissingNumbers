# Dimensions 70 and 71: cross-sections of Nebe's extremal 72-dimensional lattice

| dim | previously published | this package | factor |
|---|---|---|---|
| 71 | 331 737 984 | **2 603 658 750** | **7.85** |
| 70 | 331 737 984 | **1 249 778 250** | **3.77** |

```
python derive.py            # a few seconds, needs numpy + scipy
```

These are the largest improvements in the repository, and the reason is not that the
mathematics is deep — it is Venkov's theorem applied to a lattice everyone knows — but that
**the published tables have no entry for any dimension between 65 and 71**, so the value
being improved on is a 1998 bound inherited from dimension 64 by monotonicity.

## What "previously" means here

This is the single easiest thing to get wrong above dimension 48, and it was got wrong in
an earlier draft of this project.

- Cohn's *Table of kissing number bounds* covers dimensions **1–48, plus 72. Nothing else.**
- The Nebe–Sloane "Table of the highest kissing numbers presently known" that it replaced
  (last modified Feb 2011, retrievable from the Wayback Machine — the live page now just
  points at Cohn's) lists dimensions **1–40, 42, 44, 48, 64, 72, 80, 128** and nothing in
  between. Its dimension-64 entry is 138 458 880 (lattice) / **331 737 984** (non-lattice,
  Edel–Rains–Sloane 1998).
- OEIS does not fill the gap either: A001116 stops at n = 9, A002336 at 25, A028923 at 17.

So for every dimension from 65 to 71 the best published value is **331 737 984**, inherited
from dimension 64. The ceiling is dimension 72's 6 218 175 600, from Γ₇₂.

### Could Edel–Rains–Sloane beat that if someone evaluated it at n = 70, 71?

No, and this is checked two ways. Their construction is generic in *n* but was evaluated by
its authors only at 32, 36, 40, 44, 64, 80 and 128. It reverse-engineers exactly — both of
their large published values reproduce to the digit,

```
2^28 + 30828*2048 + 10416*16 + 128 =  331 737 984      (their dimension 64)
2^30 + 143780*2048 + 20540*16 + 160 = 1 368 532 064    (their dimension 80)
```

— and evaluating it at n = 70 and 71, *generously handing it the dimension-80
constant-weight code size*, gives at most **563 125 646**. Dimension 70 is a factor 2.2
above that and dimension 71 a factor 4.6. Independently, the construction is monotone in
*n*, so it is at most its dimension-80 value 1 368 532 064 for every n ≤ 80; dimension 71's
2 603 658 750 clears even that, beating it in *every* dimension up to 80.

## The idea

Γ₇₂ (Nebe 2012) is the extremal even unimodular lattice of dimension 72: minimum 8,
6 218 175 600 minimal vectors. Because 72 ≡ 0 (mod 24), **Venkov's theorem applies in its
strongest form** — the minimal shell, rescaled to the unit sphere, is a spherical
**11-design**.

### Dimension 71: five unknowns, six equations

For a fixed minimal vector *u*, the inner product ⟨u,v⟩ is an integer *c* with |c| ≤ 4,
because |u−v|² = 16 − 2c must be a lattice norm and hence ≥ 8; the exception is v = ±u,
where c = ±8. That is five unknowns. The design property supplies six exact moment
equations (degrees 0, 2, 4, 6, 8, 10). The system is **over-determined** and solves
consistently in non-negative integers:

| c | 0 | ±1 | ±2 | ±3 | ±4 | ±8 |
|---|---|---|---|---|---|---|
| n[c] | **2 603 658 750** | 1 512 243 200 | 280 928 256 | 13 959 168 | 127 800 | 1 |

The surplus equation is a genuine check, and it holds. Since Γ₇₂ ∩ u⊥ is a 71-dimensional
lattice of minimum 8 whose minimal vectors are exactly the n[0] minimal vectors orthogonal
to *u*, and any two of them are at angle ≥ 60°,

> **τ(71) ≥ 2 603 658 750.**

Nothing about Γ₇₂ beyond extremality is used — no coordinates, no automorphism group.

### Dimension 70: a 60-degree pair, and the two-point LP

Take two minimal vectors at 60°, i.e. ⟨u₁,u₂⟩ = 4. **Such a pair exists**, and the table
above is what proves it: n[4] = 127 800 > 0. The two-point LP of
[`common/PROOF-kpoint.md`](../../../common/PROOF-kpoint.md) — design moments *including the
mixed ones*, exact one-point marginals, antipodal symmetry, lattice translations and
lattice range — then gives, with an exact rational dual certificate,

> **τ(70) ≥ 1 249 778 250.**

`derive.py` sweeps all five possible inner products and confirms the 60° pair is the best:

| ⟨u₁,u₂⟩ | 0 | 1 | 2 | 3 | **4** |
|---|---|---|---|---|---|
| certified | 1 073 551 500 | 1 089 616 500 | 1 117 183 200 | 1 167 565 500 | **1 249 778 250** |

Unlike dimension 46, the LP minimum and maximum do **not** meet here — the bracket is
[1 249 778 250, 1 250 506 441] — so only the lower end is claimed and the size of that
cross-section is not pinned down.

### ⚠ This number changed late, because of a real bug

An earlier version of this package reported **1 170 823 530**. The k-point LP was imposing
only those moments in which *every* exponent is even, rather than every polynomial of degree
≤ 11 as Venkov's theorem allows — 56 equations at k = 3 instead of 161. Adding the missing
mixed moments raises the minimum to 1 249 778 250. Dimension 71 is unaffected: k = 1 has no
mixed moments. See the warning box in
[`common/PROOF-kpoint.md`](../../../common/PROOF-kpoint.md) §2 for the cheap test that
catches this class of error.

## History of the lower bounds in dimensions 70 and 71

| when | value | who | how |
|---|---|---|---|
| 1971– | (various) | Leech–Sloane and successors | laminated lattices, all far below |
| 1998 | **331 737 984** | Edel–Rains–Sloane [EdRS98] | their dimension-64 value, inherited upward by monotonicity — **still the best published value in dimensions 65–71** |
| 2012 | (ceiling) | Nebe [Neb12] | Γ₇₂ gives τ(72) ≥ 6 218 175 600, so 65–71 are squeezed between 331 737 984 and that |
| this work | **1 249 778 250** (dim 70), **2 603 658 750** (dim 71) | | cross-sections of Γ₇₂, from Venkov's theorem |

**No prior claim in either dimension is known to us.** Unlike dimensions 46 and 47, where
the same construction had already been carried out by others on P₄₈, we have found no
published cross-section computation for Γ₇₂. The one-point distribution above may well be
folklore — it is a five-minute calculation once you know Venkov's theorem — but the
consequence for τ(71) does not appear to have been recorded. **This is the claim in the
repository most in need of a literature check by someone who knows the area**, and if it
turns out to be known, the same is likely true of dimension 70.

Dimensions 65–69 were also attempted. They do not reach 331 737 984: see
[`../../closed/dim69-gamma72-three-point/`](../../closed/dim69-gamma72-three-point/).

## Files

```
derive.py    the whole derivation and every check
README.md    this file
```

`derive.py` imports three shared modules:

| module | what it supplies |
|---|---|
| [`common/theta.py`](../../../common/theta.py) | extremal theta series; the Venkov one-point distribution with its surplus check |
| [`common/kpoint_lp.py`](../../../common/kpoint_lp.py) | the k-point moment LP with exact rational dual certificates |
| [`common/published.py`](../../../common/published.py) | the published floors above dimension 48, with their provenance |

The mathematics is in [`common/PROOF-kpoint.md`](../../../common/PROOF-kpoint.md).

**The dimension-70 certificate is re-derived independently**, by
[`../dim68-69-gamma72-cross-sections/scripts/recheck_cert.py`](../dim68-69-gamma72-cross-sections/scripts/recheck_cert.py), which takes the integer
dual this programme settles on, maps it back to the *original* rows and re-evaluates
weak duality in `Fraction` arithmetic — none of the integralising, gcd-dividing or
common-denominator packing that earns the fast path its speed. Dimension 71 needs no
such pass: its claim is the exact size of the cross-section, a count rather than an LP
bound.

**`derive.py` validates before it claims.** It first runs the identical code on E₈ → E₇
(126), Leech → Λ₂₃ (93 150), P₄₈ → dim 47 (23 766 960), and at the two-point level on the
Leech's orthogonal (43 164) and 60° (49 896 = Λ₂₂) pairs and E₈'s orthogonal (60 = D₆) and
60° (72 = E₆) pairs. All seven are known independently and all seven come out exact.

## References

See [`../../../common/published.py`](../../../common/published.py).
