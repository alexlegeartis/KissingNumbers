# CLOSED — dimensions 17 to 24 are one construction, and it has one free integer

**Nothing is claimed here.** This package records what the whole range 17–24 actually *is*,
which component of it each published record moves, and — the point of a `closed/` package —
which components are at a ceiling that is proved, which at a ceiling that is only measured,
and which are genuinely open.

```bash
python layers.py      # the layer identity, as a structure    (~30 s)
python flats.py       # the flat layer, and dims 20 and 21 rebuilt exactly   (~80 s)
python equator.py     # the deep-hole rule and A(16,8,6) = 16  (~20 s)
```

Each ends in `ALL CHECKS PASSED` and exits non-zero otherwise.

---

## 1. The layer identity

Fix an octad `O` of the Golay code and split the 24 Leech coordinates as ℝ²⁴ = ℝ¹⁶ ⊕ ℝ⁸, the
ℝ⁸ being `O`'s. Write `p(x)` for the ℝ⁸ part of a Leech minimal vector (norm² = 32 throughout).
`layers.py` establishes, by exhaustive enumeration in integer arithmetic:

| |p(x)|² | positions it takes | how many | vectors per position |
|---|---|---|---|
| 0 | the origin | 1 | **4320** |
| 8 | the roots of 2E₈ | 240 | **512** |
| 16 | the norm-16 shell of 2E₈ | 2160 | **32** |
| 32 | the **doubled roots** — 240 of the 17520 vectors of that norm | 240 | **1** |

and no other value occurs, 24 included. The content is the last column: **the multiplicity is
constant on each shell.** That makes the kissing number of a coordinate section *additive over
positions*, so for **every** subspace `W ⊆ ℝ⁸`

```
tau( Leech ∩ (R^16 ⊕ W) )  =  4320 + 512·a(W) + 32·b(W) + c(W)
```

with `a`, `b`, `c` the number of roots, norm-16 vectors and doubled roots of 2E₈ inside `W`.
When `W` is spanned by roots and `L = E8 ∩ W` is the root lattice they generate, `c = a =
N₂(L)` and `b = N₄(L)`, which is the identity in its familiar form

```
tau(Lambda_{16+k})  =  4320 + 513·N_2(L_k) + 32·N_4(L_k),
L_1..L_8  =  A1, A2, A3, D4, D5, E6, E7, E8.
```

`layers.py` evaluates this for the nested chain A₁ < A₂ < A₃ < D₄ < D₅ < E₆ < E₇ < E₈ of
Bourbaki sub-diagrams, **both by the formula and by counting the section directly**, and both
agree with the tabulated kissing numbers of the laminated lattices at all eight values:

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| L | A₁ | A₂ | A₃ | D₄ | D₅ | E₆ | E₇ | E₈ |
| τ(Λ₁₆₊ₖ) | 5346 | 7398 | 10668 | 17400 | 27720 | 49896 | 93150 | 196560 |

### Two things to get right about it

**The 513 is an accident of counting.** It is 512 + 1 — a *tier* at a root and a *pole* at
twice that root. They are different positions at different radii; they are grouped only
because both sets happen to be indexed by the roots (`c = a` is checked separately). Written
as one coefficient the identity looks like a two-term formula fitted to eight data points. It
is not: it is a four-term decomposition with constant multiplicities, and the script prints
`a`, `b` and `c` separately for that reason.

**It is about Leech cross-sections, not about "all constructions".** The identity says what the
*laminated lattices* are. The published records in 17 through 21 are strictly larger, and the
extra points are a flat layer that is not in any of these sections.

---

## 2. What each record adds

| dim | τ(Λ) | record | difference | what the difference is |
|---|---|---|---|---|
| 17 | 5346 | 5730 (Cohn–Li) | 384 | flat layer — **arithmetic only here**, the decomposition is not verified in this package |
| 18 | 7398 | **8358** (this project) | 960 | flat layer, built and checked in `improved/dim18-bent-hexagon/` |
| 19 | 10668 | 11948 (Ho) | 1280 | flat layer, 4 × α = 4 × 320, certified in `closed/dim17-23-cohn-li-mechanism/` |
| 20 | 17400 | 19448 (Cohn–Li) | 2048 | flat layer, **rebuilt and proved maximum in `flats.py`** |
| 21 | 27720 | 29768 (Cohn–Li) | 2048 | flat layer, **rebuilt and proved maximum in `flats.py`** |
| 22 | 49896 | 49896 | 0 | the record *is* Λ₂₂ |
| 23 | 93150 | 93150 | 0 | the record *is* Λ₂₃ |
| 24 | 196560 | 196560 | 0 | the record *is* the Leech lattice, and it is optimal |

---

## 3. The flat layer, exactly (`flats.py`)

Cohn–Li's second family is, in ℝⁿ, `4·C(n,2)` vectors `(±2,±2,0ⁿ⁻²)` plus `128·|C|` vectors
`(±1)⁸` on the blocks of the shortened Golay code with an **odd** number of minus signs, plus
flats `√(8/n)·(±1)ⁿ`. Everything has norm² = 8, so the condition is pairwise inner product ≤ 4,
and two derivations follow from that alone:

* a flat meets a block in `√(8/n)(8−2k)`, `k` = disagreements, and the block carries *every*
  odd pattern, so `k = 0` occurs unless `⟨c, O⟩ = 0`. The blocks span the shortened Golay code
  (checked: its rank is `12 − (24−n)`), the Golay is self-dual, and shortening-then-dualising
  is puncturing — so the flat words are exactly the **(24−n)-punctured Golay code**, 4096 words;
* two flats at Hamming distance `d` meet in `(8/n)(n−2d) ≤ 4` iff `d ≥ n/4`.

So the flat layer is a maximum independent set in `Cay(F₂¹², S)`, `S` = the punctured code's
nonzero words of weight `< ⌈n/4⌉`. That graph is computed, and in two dimensions its structure
settles α outright:

| n | \|S\| | the graph | α |
|---|---|---|---|
| 19 | 21 (one of weight 3, twenty of weight 4) | not bipartite; four components of 1024 | 4 × 320 = **1280** (elsewhere) |
| 20 | 5, independent, all weight 4 | **128 disjoint 5-cubes**, 5-regular bipartite | **2048** |
| 21 | 21, all weight 5 | four components of 1024, 21-regular **bipartite** | **2048** |
| 22, 23 | 0 — the punctured code already has `d_min ≥ ⌈n/4⌉` | — | 4096 |

A regular bipartite graph has a perfect matching, so by **König** α is exactly half the
vertices: 2048 in both 20 and 21, which is Cohn–Li's value. The script then **builds both
configurations from the Golay code and checks every pair in integer arithmetic** — base against
base (chunked Gram), base against flat (`⟨u,s⟩ ≤ ⌊√(2n)⌋`), flat against flat (`d ≥ ⌈n/4⌉`) —
and reaches 17 400 + 2048 = **19 448** and 27 720 + 2048 = **29 768** exactly.

**Where this family stands relative to the lattice, dimension by dimension.** Its base equals
τ(Λₙ) *only* for n = 19, 20, 21. In 17 and 18 the base is smaller (the record there is the
other Cohn–Li family, ℝ¹⁶ ⊕ ℝᵏ over the odd Barnes–Wall lattice); in 22 and 23 it is far
smaller, and even with the whole 4096-word ambient as a flat layer it reaches 47 260 and
69 876 against the lattice's 49 896 and 93 150. That — not an empty flat layer — is why this
family gives nothing in 22 and 23.

One more caution the script prints: for n ≤ 18 the size of `C` depends on *which* coordinates
are deleted (25 or 30, 45 or 46), because the Golay octads are only a 5-design. For n ≥ 19 the
count is independent of the choice: 78, 130, 210.

---

## 4. The deep hole fixes the radii (`equator.py`)

A layer at radius `r` over the equator carries vectors `x = (v, w)` with `|w|² = r²`,
`|v|² = 8 − r²`, and needs `⟨v/|v|, ê⟩ ≤ 4/√(8(8−r²))` against every equator point. So with
`m(E) = min_z max_e ⟨z, ê⟩`, a layer at radius `r` can exist **only if**

```
r^2  >=  8 - 2/m(E)^2 .
```

An explicit direction certifies the side that says a layer *does* open, and both certificates
are exact integer conditions:

| equator | certificate | m(E) ≤ | smallest radius² | which layer |
|---|---|---|---|---|
| even BW₁₆ | `z = v/(2√2)`, `v` a ±1 vector on a six-set, `⟨v,e⟩ ≤ 8` for all 4320 | `1/√3` | **2** | the **tier** — and `|v|² = 6` *is* a tier's v-part |
| odd BW₁₆ | `z = s√2/6`, `s ∈ {±1}¹⁶`, `⟨s,e⟩ ≤ 12` for all 4320 | `3√2/8` | **8/9** | the **flats** — v-part `(2/3)(±1)¹⁶` |

**The two hole sets are one code, and the script names it rather than asserting it.** Over the
odd shell exactly 2048 sign vectors work; their minus-sign indicators are checked to be
*linear*, of length 16 and dimension 11, with weight distribution 1/140/448/870/448/140/1 —
minimum distance 4, so the code is the extended Hamming code of length 16, which is unique:
**RM(2,4)**. Over the even shell exactly 28 672 = **448 × 64** six-set directions work, and
those 448 six-sets are checked to be *exactly that same code's weight-6 words* — while **no**
all-±1 direction works over the even shell at all. (28 672 is also the tier ground set of the
dimension-18 work, which is where the two meet.) So Cohn–Li's sign flip is precisely "deepen
the equator's holes by 8.2%", and it buys exactly the one layer that opens.

> **Not proved here.** That these `m(E)` are *minimal* — hence that no layer sits closer in —
> is a minimax. It is measured by LP ascent (`research/dim20-21/maximal2.py`) at exactly
> 0.577350269 and 0.530330086, and reported as measured. The certificates above are one-sided.
> A method warning that cost a day: a smoothed log-sum-exp minimax returns 0.5807 on the even
> shell where the truth is 1/√3 = 0.57735. Use the LP ascent, and calibrate it on that shell.

### The equator is free for the tier

Without an equator a tier of full parity classes could use any six-subset of [16] and would be
`32 · A(16,8,6)`.

> **Theorem. A(16,8,6) = 16.** With `b` blocks of size 6 pairwise meeting in ≤ 2 and `r_p` the
> number through `p`: `Σ_p r_p = 6b`, and for each block `X`,
> `Σ_{p∈X} r_p − 6 = Σ_{Y≠X} |X ∩ Y| ≤ 2(b−1)`, so `Σ_{p∈X} r_p ≤ 2b + 4`. Summing over blocks,
> `Σ_p r_p² ≤ b(2b+4)`, while Cauchy–Schwarz gives `Σ_p r_p² ≥ (6b)²/16 = 9b²/4`. Hence
> `9b²/4 ≤ 2b² + 4b`, i.e. **b ≤ 16**. ∎

`equator.py` exhibits the attainment: the sixteen weight-6 words of the bent coset
`x₁x₂ + x₃x₄ + RM(1,4)`, which meet pairwise in exactly 2 — the 2-(16,6,2) biplane. So a tier
is 512 with or without an equator, and **deleting equator points to make room for more six-sets
cannot help**. That closes what had looked like the most promising remaining lever.

---

## 5. What is actually left

Every component above is pinned except one: **the tier size**.

* 512 is attained in two structurally different ways, and both are maximal
  (`research/dim18/TIER.md`, `research/dim22/tiermax22.py`: 0 of 61 440 candidates addable,
  with a drop-one calibration that recovers exactly the removed vector).
* The 30-class orbital association scheme of the 28 672-vertex conflict graph gives
  **T ≤ 154112/289 = 533.26** with an exact rational certificate, so `T ≤ 533`.
* `T ≥ 513` would improve **every one** of dimensions 17 through 23 at once, by `N₂(L_k)` =
  2, 6, 12, 24, 40, 72, 126 respectively.

**But dimension 24 says a free tier of 513 cannot exist.** The k = 8 instance of the identity
is 196 560, which is *proved optimal* (Odlyzko–Sloane, Levenshtein). If some tier in the full
E₈ arrangement could hold 513 points while staying compatible with the other 239 tiers, the
2160 half-tiers and the 240 poles, the result would be a 196 561-point kissing configuration in
ℝ²⁴. So **any 513-point tier must be incompatible with the rest of the E₈ arrangement.** That
does not close 18 through 23, where far fewer positions have to be satisfied simultaneously —
but it does say the tier cannot be enlarged locally and for free, and it is the calibration any
future tier search should be run against.

Two further cautions on the `T ≤ 533` bound, both recorded in `research/dim18/TIER.md`:

* it bounds the natural **discrete** model. The admissible set is not finite — each flat has
  slack 0.326 and carries an open neighbourhood of admissible non-uniform profiles — so only
  the continuous spherical LP bound, 672, covers everything;
* it is computed over the **odd** equator. Over the even one (which is what dimensions 22 and
  23 use) the ground set is different, and the only bound available there is the continuous
  `A(16,1/3) ≤ 751`. 2281 large-neighbourhood iterations over the full 61 440-vertex ground
  set found nothing above 512.

---

## Files

| file | what it does | time |
|---|---|---|
| `layers.py` | the layer identity as a structure: multiplicities, position sets, the chain both ways, and a 2000-restart search for a better section | ~30 s |
| `flats.py` | the flat layer as a Cayley-graph independence number; dimensions 20 and 21 rebuilt from the Golay code and checked pair by pair; and \|C\| over **every** deletion set | ~80 s |
| `equator.py` | the deep-hole radius rule with exact one-sided certificates, and A(16,8,6) = 16 with its attainment | ~20 s |

All three import only `numpy` and this repository's `common/leech.py` and `common/golay.py`.

## Working notes

`research/notes/KNOWLEDGE.md` section 93; `research/dim22/NOTES.md`; `research/dim18/TIER.md`;
`research/noequator/NOTES.md`; `research/dim20-21/flatcap.py`.
