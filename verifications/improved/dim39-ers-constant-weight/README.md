# Dimension 39: Edel-Rains-Sloane with the 2026 constant-weight codes

$$τ_{39} \ge 756116$$

improving the published record 755988 by **+128**.

The same script also derives $τ_{38} \ge 570236$, a +3584 improvement on 566652 that
held the record briefly during this project.  It is **superseded** by the Leech cap
construction of [`../dim38-leech-large-codimension/`](../dim38-leech-large-codimension/),
which gives 591612 -- see
[`../../superseded/dim38-ers-constant-weight/`](../../superseded/dim38-ers-constant-weight/).
Both dimensions are checked here because they share a single argument and a single script;
only dimension 39 is claimed from this package.

## History of the lower bound in dimension 39

| when | value | who |
|---|---|---|
| 1998 | 755 988 | Edel-Rains-Sloane [EdRS98], the chain (39, 8, 2) evaluated with the constant-weight codes of the day -- **still the entry in Cohn's table and in Brouwer's** |
| Aug 2026 | (skipped) | Echols [Ech26], arXiv:2608.13906, improves A(39,8,8) from 3323 to 3324 but derives his kissing numbers from a fixed level-1 length n0 = 32, which in dimension 39 gives less than the record -- so this dimension was passed over |
| this work | **756 116** | the same chain with n0 = n = 39 and Echols' A(39,8,8) >= 3324 |

The gain is exactly 128 x (3324 - 3323) = 128.  Dimensions 38 and 39 are the only two where
the ERS chain uses n0 = n rather than n0 = 32, which is precisely why they are the two that
Echols' fixed-n0 formula misses.

**Check before citing.** Sun-Wang [SW26] (arXiv:2607.20359) improve dimension 39 by the
antipode construction in a later revision than the one consulted here.  This is the claim in
the repository most likely to have already moved.

```
python scripts/verify.py
```

runs in about ten seconds, uses only integer and exact-rational arithmetic, and needs no
third-party packages.

## Where the improvement comes from

Both records come from the Edel–Rains–Sloane construction
(*On kissing numbers in dimensions 32 to 128*, Electron. J. Combin. **5** (1998) R22,
doi:10.37236/1360, arXiv:math/0207291).

Fix a dimension $n$ and a chain of support sizes $n \ge n_0 \ge 4n_1 \ge 16 n_2 \ge \dots$.
The level-$\nu$ centres are the vectors of shape $\pm a_\nu^{\,n_\nu} 0^{\,n-n_\nu}$ with
$a_\nu = \sqrt{n_0/n_\nu}$, supports drawn from a constant-weight code
$\mathcal C(n,n_\nu,n_\nu)$ and signs from a code $\mathcal C(n_\nu,\lceil n_\nu/4\rceil)$.
Every centre has squared norm $n_0$ and every two are at squared distance $\ge n_0$, so
the whole set is a 60° spherical code and

$$\tau_n \;\ge\; \sum_\nu A(n,n_\nu,n_\nu)\,A\!\left(n_\nu,\lceil n_\nu/4\rceil\right).$$

Dimensions 38 and 39 use the three-level chain $(n,8,2)$, i.e.

$$\tau_n \;\ge\; A\!\left(n,\lceil n/4\rceil\right) \;+\; 128\,A(n,8,8) \;+\; 4\binom{n}{2}.$$

The middle term is the one that moved. W. Echols, *New lower bounds for constant-weight
codes via seeded bit-swap tabu search*, arXiv:2608.13906 (2026), improved $A(n,8,8)$ for
several $n$; Brouwer's table now carries

| $n$ | value used by the published record | current best known |
|-----|-----------------------------------|--------------------|
| 38  | 2997                              | **3025**           |
| 39  | 3323                              | **3324**           |

Echols' own paper derives its kissing numbers from a fixed $n_0 = 32$, which for $n = 38$
gives only $131072 + 128\cdot 3025 + 2812 = 521084$ — *below* the record — so dimensions 38
and 39 were passed over. They are exactly the two dimensions in which the ERS chain uses
$n_0 = n$ rather than $n_0 = 32$, with the much larger level-1 codes
$A(38,10) \ge 180224$ and $A(39,10) \ge 327680$. Substituting the new constant-weight codes
into *that* chain gives

$$\tau_{38} \ge 180224 + 128\cdot 3025 + 4\cdot 703 = 180224 + 387200 + 2812 = 570236,$$
$$\tau_{39} \ge 327680 + 128\cdot 3324 + 4\cdot 741 = 327680 + 425472 + 2964 = 756116.$$


`A(39,8,8) = 3324` was **re-fetched live from Brouwer's page on 2026-08-25 and again on
2026-08-30, the day of release**, and both entries still read 3025 and 3324, still
marked `E` for Echols. It still reads
3324, alongside 2832, 3025, 3683 and 4510 at n = 37, 38, 40, 41 - the same values the cached
`Andw.html` in `../../closed/dim32-44-ers-audit/data/` carries. The whole +128 is that one
entry, so it is the thing to re-check whenever this dimension is quoted.

The gains are exactly $128\cdot(3025-2997) = 3584$ and $128\cdot(3324-3323) = 128$.

## What the script verifies, and what it assumes

**Verified computationally, exhaustively:**

* `data/a38.8.8.3025H` contains 3025 distinct binary words of length $\le 38$, every one of
  weight 8, with maximum pairwise support intersection 4, hence minimum Hamming distance
  exactly 8. Same for `data/a39.8.8.3324H` with 3324 words of length $\le 39$.
  (All $\binom{3025}{2} = 4{,}573{,}800$ and $\binom{3324}{2} = 5{,}522{,}826$ pairs are
  compared.) These files are Yves Edel's, taken from Brouwer's table
  <https://aeb.win.tue.nl/codes/Andw.html>.
* the level-2 sign code — the even-weight binary code of length 8 — has 128 words and
  minimum distance $2 = \lceil 8/4\rceil$;
* the level-3 data — all $\binom n2$ two-element supports and all 4 sign patterns of
  length 2, minimum distance $1 = \lceil 2/4 \rceil$;
* the chain condition $n \ge n_0 \ge 4n_1 \ge 16 n_2$;
* **all ten geometric conditions**, as exact rational inequalities: for each ordered pair
  of levels, that the largest possible inner product is at most $n_0/2$. Three of them are
  tight at equality ($19$ vs $19$ in dimension 38, $39/2$ vs $39/2$ in dimension 39), which
  is what makes the chain condition necessary rather than decorative;
* the arithmetic of the three totals.

**Taken from the literature, not verified here:**

* $A(38,10) \ge 180224$ and $A(39,10) \ge 327680$, due to
  V. A. Zinov'ev and S. N. Litsyn, *On shortening of codes*, Probl. Peredachi Inf. **20**:1
  (1984) 3–11 (Engl. transl. Problems Inform. Transmission **20**:1 (1984) 1–7), tabulated
  by Litsyn–Rains–Sloane. No explicit codewords for these appear to be published anywhere.

  This is not a new assumption: the **existing** records 566652 and 755988 in H. Cohn's
  table rest on exactly the same two numbers, with the same provenance. The improvement
  claimed here is entirely in the level-2 term and is fully machine-checked.

  For reference, the best known *linear* codes give only $[38,17,10]$ and $[39,18,10]$
  (Grassl, codetables.de), i.e. $2^{17}$ and $2^{18}$; substituting those yields
  $\tau_{38} \ge 521084$ and $\tau_{39}\ge 690580$, which are unconditional but below the
  current records.

## Context: these are the two dimensions nobody has updated

As of 2026-08-19:

| dim | Cohn's table | Brouwer's table | correct value | who claimed it |
|-----|--------------|-----------------|---------------|----------------|
| 32 | 345408 | 346432 | 346432 | Echols 2026 |
| 33 | 360640 | 362048 | 362048 | Echols 2026 |
| 34 | 380868 | 381124 | 381124 | Echols 2026 |
| 37 | 494312 | 496232 | 496232 | Echols 2026 |
| **38** | 566652 | 566652 | **591612** | this project, but by a different construction — the ERS value here is 570236 |
| **39** | 755988 | 755988 | **756116** | this project |

Dimensions 32, 33, 34 and 37 are stale only in Cohn's table; Brouwer already carries the
corrected values. Dimensions 38 and 39 are stale in **both**, because both tables inherit
Echols' fixed-$n_0$ formula.

Dimensions 35, 36, 40, 41, 42, 44, 45, 46 and 47 give no gain: Brouwer's current
$A(n,8,8)$ is identical to the value already used there. Dimension 43 is held by a
different construction (Sun–Wang, arXiv:2607.20359, and the 1982 Conway–Sloane
cross-section of $P_{48}$), against which ERS is not competitive.

The full audit of this construction across dimensions 32–44 — maximised over **all**
admissible chains rather than assuming $n_0 = n$, which matters, because dimensions 41 and
42 optimise at $n_0 = 40$ — is in
[`../../closed/dim32-44-ers-audit/`](../../closed/dim32-44-ers-audit/). It also identifies
exactly where a future code improvement would pay off: every dimension except 38, 43 and 44
sits *at* the ERS value, so $+1$ in $A(n,8,8)$ is immediately $+128$ in $\tau_n$.

## Files

```
data/a38.8.8.3025H     Brouwer/Edel constant-weight code, A(38,8,8) >= 3025
data/a39.8.8.3324H     Brouwer/Edel constant-weight code, A(39,8,8) >= 3324
scripts/verify.py      the checker
```

The code files are in Brouwer's format: a `$BASE=16` header, then one codeword per line as
a hexadecimal integer whose bit $i$ is coordinate $i$.
