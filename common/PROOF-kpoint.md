# Cross-sections of extremal even unimodular lattices: the k-point moment method

This is the mathematics behind two of the claims in this repository, and behind the two
values it recovers rather than claims:

```
claimed        tau(71) >=  2 603 658 750       (k = 1, dimension 72)
               tau(70) >=  1 249 778 250       (k = 2, dimension 72)

recovered      tau(47) >=       23 766 960     (k = 1, dimension 48)
               tau(46) >=       12 309 600     (k = 2, dimension 48)
```

The two dimension-48 numbers were published before this project reached them —
Boyvalenkov–Cherkashin, *Results in Mathematics* **80** (2025), Paper No. 3, and Sun–Wang,
arXiv:2607.20359 — and are derived here because agreeing with a published count, from a
method that uses no coordinates and knows nothing of either paper, is the only external
check this machinery has. Nothing about them is claimed.

Every numerical step is verifiable in exact rational arithmetic by
[`kpoint_lp.py`](kpoint_lp.py) and [`theta.py`](theta.py), and by the `derive.py` of the two
packages
[`verifications/recovered/dim46-47-p48-cross-sections/`](../verifications/recovered/dim46-47-p48-cross-sections/)
and
[`verifications/improved/dim70-71-gamma72-cross-sections/`](../verifications/improved/dim70-71-gamma72-cross-sections/).

---

## 0. Notation and the three classical inputs

Let *L* be an extremal even unimodular lattice of dimension *n* ≡ 0 (mod 24), Λ its set of
minimal vectors, μ its minimum and *N* = |Λ|.

**(C1)** The theta series of *L* is the unique modular form of weight *n*/2 for SL₂(ℤ) with
constant term 1 and vanishing coefficients of q¹ … q^(n/24). Hence (μ, N) is determined by
*n* alone:

| *n* | 8 | 24 | 48 | 72 |
|---|---|---|---|---|
| μ | 2 | 4 | 6 | 8 |
| *N* | 240 | 196 560 | 52 416 000 | 6 218 175 600 |

*Computed from E₄ and Δ in `theta.py`. For n = 48 four such lattices are known — P₄₈p, P₄₈q
(Conway–Sloane, SPLAG ch. 7), P₄₈m, P₄₈n (Nebe, Discrete Math. 331 (2014) 133–136); for
n = 72 one, Nebe's Γ₇₂ (J. Reine Angew. Math. 673 (2012) 237–247). Nothing below
distinguishes them.*

**(C2)** Λ, rescaled to the unit sphere, is a spherical **11-design**.
*Venkov, "Even unimodular extremal lattices", Trudy Mat. Inst. Steklov 165 (1984) 43–48;
reproduced as Theorem 4.1 of Boyvalenkov–Cherkashin, arXiv:2312.05121. For n = 8 it is a
7-design, which is all the E₈ validation needs.*

**(C3)** For distinct *v*, *w* ∈ Λ with *w* ≠ −*v*, the inner product *v*·*w* is an integer
with |*v*·*w*| ≤ μ/2; and *v*·*w* = ±μ exactly when *w* = ±*v*.
*Integrality because L is integral; the range because |v−w|² = 2μ − 2v·w is the norm of a
non-zero lattice vector, hence an even integer ≥ μ, and likewise for |v+w|².*

---

## 1. A cross-section is a kissing configuration

**Lemma 1.** Let *u*₁, …, *u*_k ∈ Λ be linearly independent and put
*S* = {*v* ∈ Λ : ⟨*u_i*, *v*⟩ = 0 for all *i*}. Then *S* lies in an (*n*−*k*)-dimensional
subspace, its members all have the same length, and any two distinct *v*, *w* ∈ *S* satisfy
cos(*v*,*w*) ≤ 1/2. Hence **τ(*n*−*k*) ≥ |*S*|**.

*Proof.* *S* lies in the orthogonal complement of span(*u*₁…*u*_k), of dimension *n*−*k*.
All members have norm μ. For distinct *v*, *w* ∈ *S*, *v* − *w* is a non-zero lattice
vector, so |*v*−*w*|² = 2μ − 2*v*·*w* ≥ μ, i.e. *v*·*w* ≤ μ/2, i.e. cos ≤ 1/2. A set of
equal-length vectors with pairwise angle ≥ 60° is a kissing configuration. ∎

For *k* = 1 and *n* = 48 this is the classical construction, and |*S*| = 23 766 960 is
recorded in arXiv:2312.05121, equation (3).

---

## 2. The moment identities

**Lemma 2.** Fix *u*₁…*u*_k ∈ Λ and let *n*[*c*₁…*c*_k] be the number of *v* ∈ Λ with
⟨*u_i*,*v*⟩ = *c_i* for every *i*. Then for **every** exponent vector *p* with
Σ*p_i* ≤ 11 and Σ*p_i* even,

```
    sum_cells  n[c] * prod_i (c_i / mu)^{p_i}
        = N * E[ prod_i <u_i-hat, Z>^{p_i} ] / ( n (n+2) ... (n + sum p - 2) )
```

where *Z* is a standard Gaussian in ℝⁿ and the numerator is evaluated by Wick's theorem
from the normalised Gram matrix of the *u_i*.

*Proof.* The monomial ∏⟨û_i, x⟩^{p_i} has degree Σp_i ≤ 11, so by (C2) its sum over the
rescaled Λ equals *N* times its average over the sphere. Writing *x* = *Z*/|*Z*| and using
independence of |*Z*| and *Z*/|*Z*| gives the normalisation E[|Z|^{2m}] = n(n+2)…(n+2m−2).
The Gaussian moment is a sum over perfect matchings, each pair contributing a Gram entry
(Wick's theorem). Odd total degree vanishes by antipodal symmetry. ∎

> **⚠ The mixed moments are not optional.** It is tempting to impose only those *p* in which
> every exponent is even — that is the family you get by squaring, and it is what an earlier
> version of this code did. It is wrong, and it is *silently* wrong: it gives a valid but
> weaker LP whose answer depends on the presentation of the Gram matrix. The symptom is
> decisive and cheap to test for: the Gram matrices (3,−3,0) and 3·A₃ are equivalent
> (UᵀGU = 3A₃ for an integer *U* of determinant −1) and must give the same answer. With
> even-exponent moments only, they gave 6 462 480 and 5 194 969. With the full family, both
> give 6 687 320. For *k* = 3 the correct family has **161 equations, not 56**. These are
> exactly what Ozeki obtains by expanding Σ_c l_c (au+bv+ct)^{2k}.

By (C3) the sum runs over *c* ∈ {0, ±1, …, ±μ/2, ±μ}^k, restricted to those *c* for which
the Gram matrix of (*u*₁…*u*_k, *v*) is positive semidefinite. Cells with |*c_i*| = μ mean
*v* = ±*u_i* and contribute exactly 1 each; cells where that Gram matrix is **singular**
force *v* into span(*u*₁…*u*_k) and so contribute **at most 1**.

**Lemma 2b (translation constraints).** If *c_i* = μ/2 then |*u_i* − *v*|² = μ, so
*w* := *u_i* − *v* is again a minimal vector, with ⟨*u_i*,*w*⟩ = μ/2 and
⟨*u_j*,*w*⟩ = *G_ij* − *c_j*. Since *v* ↦ *u_i* − *v* is an involution,

```
    n[c] = n[c']    where  c'_i = mu/2,  c'_j = G_ij - c_j   (j != i),
```

and symmetrically *n*[*c*] = *n*[*c*″] via *w* = *u_i* + *v* when *c_i* = −μ/2, with
*c*″_j = *G_ij* + *c_j*. If the image cell is not a legal cell the constraint reads
*n*[*c*] = 0; if it is a ±*u_j* cell it reads *n*[*c*] = 1.

*Proof.* Immediate from |*u_i* − *v*|² = 2μ − 2*c_i* = μ and bilinearity. ∎

These use the *lattice* structure, which the design property alone does not see, and they
are strong: on the Leech lattice they turn the *k* = 3 sixty-degree-triple bound from 26 946
into the exact value 27 720, and they raise dimension 46 from 12 086 040 to 12 309 600.

**Lemma 2c (lattice range).** For every lattice vector *w* = Σ*a_i u_i* in the span and
every minimal *v*, |*v*−*w*|² = μ + |*w*|² − 2⟨*v*,*w*⟩ is a lattice norm, so
|⟨*v*,*w*⟩| ≤ |*w*|²/2 unless *v* = ±*w* (possible only when |*w*|² = μ). A cell violating
this for some *w* is empty.

*Proof.* |*v*−*w*|² ≥ μ for *v* ≠ *w*, and likewise for |*v*+*w*|². ∎

Strictly stronger than positive semidefiniteness: for a 60° pair, *u*₁−*u*₂ is itself
minimal, forcing |*c*₁−*c*₂| ≤ μ/2 and killing ten cells that PSD permits. Those ten are
exactly the zeros of Ozeki's Lemma 4.2, independently derived.

Two further exact constraints are used: the **one-point marginals** (Lemma 2 with *k* = 1,
solved in closed form — see §4), and **antipodal symmetry** *n*[*c*] = *n*[−*c*], since
Λ = −Λ.

---

## 3. The bound, and why floating point cannot corrupt it

Collect the constraints of §2 as *Ax* = *b*, 0 ≤ *x* ≤ *u* (with *u_j* = 1 on the singular
cells and +∞ elsewhere), *x* indexed by cells, and let *c* be the indicator of the all-zero
cell. The true distribution is feasible, so

```
    |S| = x_true[0..0]  >=  min { c^T x : A x = b, 0 <= x <= u }.
```

**Lemma 3 (weak-duality certificate).** For any rational vector *y*, put *s* := *c* − *A*ᵀ*y*.
If *s_j* ≥ 0 for every *j* with *u_j* = +∞, then

```
    min { c^T x : A x = b, 0 <= x <= u }
        >=  b^T y  -  sum_{j : u_j finite} u_j * max(0, -s_j).
```

*Proof.* For feasible *x*, *c*ᵀ*x* = *y*ᵀ*Ax* + *s*ᵀ*x* = *b*ᵀ*y* + *s*ᵀ*x*, and
*s*ᵀ*x* ≥ −Σ_{u_j finite} *u_j* max(0, −*s_j*) because 0 ≤ *x_j* ≤ *u_j*. ∎

`kpoint_lp.py` obtains *y* from a floating-point LP, rounds it to rationals, restores the
hypothesis of Lemma 3 by decreasing the dual coordinate of the all-ones moment row (whose
coefficients are all strictly positive), and then evaluates the right-hand side **in exact
rational arithmetic**, rounding to an integer in the safe direction. An imprecise *y* can
only weaken the bound, never invalidate it. `run3` returns the certified integer as its
third component; that, never the floating-point LP value, is what is claimed.

The *k* ≥ 3 programs are badly conditioned (moment coefficients span 6⁻¹⁰) and HiGHS
converges only on the relaxed form |Ax − b| ≤ ε. Its dual is automatically feasible for the
equality problem, which is exactly what the certification needs.

---

## 4. The one-point distribution in closed form

For *k* = 1 the system of Lemma 2 has μ/2 + 1 unknowns and 6 equations, so it is
over-determined and can be solved by exact Gauss–Jordan; the surplus equations are a genuine
check. `theta.py` does this and finds:

| lattice | *n* | μ | n[0] | n[1] | n[2] | n[3] | n[4] | surplus eqs |
|---|---|---|---|---|---|---|---|---|
| E₈ | 8 | 2 | 126 | 56 | | | | 2 |
| Leech | 24 | 4 | 93 150 | 47 104 | 4 600 | | | 3 |
| P₄₈ | 48 | 6 | **23 766 960** | 12 608 784 | 1 678 887 | 36 848 | | 2 |
| Γ₇₂ | 72 | 8 | **2 603 658 750** | 1 512 243 200 | 280 928 256 | 13 959 168 | 127 800 | 1 |

The first two rows are E₇ (126) and Λ₂₃ (93 150), both known independently — that is the
validation. The third and fourth rows are dimensions 47 and 71.

---

## 5. Realizability

A bound for a prescribed Gram matrix is vacuous unless *k* minimal vectors with that Gram
actually occur in *L*. The same machinery proves this.

**Lemma 4.** If the minimum of *n*[*c*₁…*c*_{k−1}] over the (*k*−1)-point system is strictly
positive, then some *v* ∈ Λ has ⟨*u_i*,*v*⟩ = *c_i* for all *i*, so the *k*-tuple exists.

*Proof.* The true distribution is feasible, so its value at that cell is at least the
minimum, which is positive. ∎

| configuration | witnessing cell | certified minimum | exists |
|---|---|---|---|
| P₄₈, *u*₁ | one-point n[0] | 23 766 960 | yes |
| P₄₈, *u*₁·*u*₂ = 3 (60°) | one-point n[3] | 36 848 | yes |
| Γ₇₂, *u*₁ | one-point n[0] | 2 603 658 750 | yes |
| Γ₇₂, *u*₁·*u*₂ = 4 (60°) | one-point n[4] | 127 800 | yes |

A sweep over all Gram matrices produces larger apparent bounds whose witnessing cells have
LP minimum 0 — most temptingly 527 597 644 in dimension 69, from the Γ₇₂ Gram (4,−2,0).
Lemma 4 cannot license those. It is not needed for them either: realizability is a fact about
the lattice, and where the lattice is available in coordinates a k-tuple with a prescribed Gram
can be **exhibited** instead, which is what
[`verifications/improved/dim68-69-gamma72-cross-sections/`](../verifications/improved/dim68-69-gamma72-cross-sections/)
does for dimensions 68 and 69 — and it exhibits the *determinant-minimal* Gram, A₃ at det 256,
which Hermite shows is the best a triple can have. See
[`verifications/closed/dim69-gamma72-three-point/`](../verifications/closed/dim69-gamma72-three-point/)
for the withdrawal this replaced.

---

## 6. The results, and what "exact" means

| dim | *n* | *k* | Gram | certified \|S\| | previously published | factor |
|---|---|---|---|---|---|---|
| 47 | 48 | 1 | — | **23 766 960** | 9 741 412 | 2.44 |
| 46 | 48 | 2 | *u*₁·*u*₂ = 3 | **12 309 600** | 5 318 060 | 2.31 |
| 46 | 48 | 2 | orthogonal | [10 659 120, 10 697 520] | | |
| 71 | 72 | 1 | — | **2 603 658 750** | 331 737 984 | 7.85 |
| 70 | 72 | 2 | *u*₁·*u*₂ = 4 | **1 249 778 250** | 331 737 984 | 3.77 |

**A point of care about the word "exact".** Everything in the table is a **lower bound on
the kissing number**. What is exact in three of the rows is the *size of the particular
configuration*, not τ of that dimension:

* for *k* = 1 the linear system has a unique solution, so those cross-sections contain
  **exactly** 23 766 960 and 2 603 658 750 minimal vectors;
* for dimension 46 the LP minimum and maximum coincide, so that cross-section contains
  **exactly** 12 309 600 minimal vectors (31 of the system's 47 cells are pinned to single
  values; the 16 exceptions all have both coordinates in {±1,±2});
* for dimension 70 they do not coincide — the bracket is [1 249 778 250, 1 250 506 441] —
  so only the lower end is claimed;
* for the orthogonal pair in dimension 46 the size is only bracketed, and Ozeki's degree-3
  table shows it genuinely **varies from pair to pair**, so the bracket cannot be collapsed.
  That is geometry, not a defect of the method.

In every case the consequence for the kissing number is an **inequality**,
τ(*n*−*k*) ≥ (size of the configuration). τ(46), τ(47), τ(70) and τ(71) themselves remain
unknown, unlike τ(8) = 240 and τ(24) = 196 560, which are theorems.

---

## 7. Independent confirmation

**Dimension 47.** The distribution of §4 is equation (3) of Boyvalenkov–Cherkashin,
arXiv:2312.05121 (Dec 2023), which states in the same paragraph that these A₀ vectors
"define a 47-dimensional kissing configuration". Sun–Wang, arXiv:2607.20359v3 (Aug 2026),
Table 3, obtain 23 766 960 by sweeping all 52 416 000 minimal vectors directly.

**Dimension 46.** Ozeki, *Tsukuba J. Math.* 40 (2016) 139–186, Table 3, gives the degree-3
Siegel Fourier coefficients, whose values count triples of lattice vectors with a prescribed
Gram matrix:

```
a( (3,3,3,3,0,0) ) = 23775066324172800000        a( (3,3,3) ) = 1931424768000
23775066324172800000 / 1931424768000 = 12309600      exactly, no remainder
```

His degree-2 table likewise reproduces the one-point distribution of §4 exactly:
a((3,3,c)) = 52 416 000 × A_c for c = 3, 2, 1, 0. Sun–Wang, Table 3 (complement 3A₂), give
12 309 600 by direct computation.

**Ozeki's Theorem 6.3, rederived.** With the mixed moments, the three-point system for a
60° triple in an extremal 48-dimensional lattice has affine solution space of dimension
exactly **1**, with

```
n[0,0,0] + 164 * n[3,3,3] = 6732912   (per triple, exact)
```

Running the LP on the cell n[3,3,3] of the Gram 3·A₃ certifies min = **184**, and on the
all-zero cell certifies the bracket **[6 687 320, 6 702 736]**. The two are exactly the two
ends of that line: 6 732 912 − 164·278 = 6 687 320 and 6 732 912 − 164·184 = 6 702 736, with
nothing lost to the outward rounding of a maximisation certificate. So non-negativity forces 184 ≤ n[3,3,3] ≤ 278, which is Ozeki's bound on
τ = a(T₄,3)/a(T₃,1) — obtained here by linear programming with an exact rational dual
certificate rather than by modular forms. **Sun–Wang's two dimension-45 values, 6 687 648
(3A₃) and 6 702 080 (3D₃), both lie inside that bracket**, which is a sharp independent
check on the machinery.

**The validation suite.** Run on lattices whose cross-sections are known, `kpoint_lp.py`
reproduces every one exactly:

| lattice | *k* | Gram | LP | true |
|---|---|---|---|---|
| E₈ | 1 | — | 126 | 126 (E₇) |
| E₈ | 2 | orthogonal | 60 | 60 (D₆) |
| E₈ | 2 | 60° | 72 | 72 (E₆) |
| Leech | 1 | — | 93 150 | 93 150 (Λ₂₃) |
| Leech | 2 | orthogonal | 43 164 | 43 164 |
| Leech | 2 | cos 1/4 | 44 550 | 44 550 |
| Leech | 2 | 60° | 49 896 | 49 896 (Λ₂₂) |
| Leech | 3 | 60° triple | 27 720 | 27 720 (Λ₂₁) |
| Leech | 3 | orthogonal | 19 530 | 19 962 — a valid bound, 2.2% short |

The right-hand column is computed directly from explicit Leech coordinates.

---

## 8. Where this stops, and what would go further

Dimensions 44 and 45 are where the method stops being competitive, though not where it stops
being informative. At k = 3 the Gram 3A₃ certifies 6 687 320 against Sun–Wang's 6 702 080;
at k = 4 the Gram 3D₄ certifies 4 147 967 against their 4 153 600. Both of their values land
*inside* the brackets computed here — [6 687 320, 6 702 736] and [4 147 967, 4 165 632] — which
is a sharp check on the method but not a competitive bound. (Use the right Gram: 3A₄, a
different lattice, gives only [3 694 470, 3 833 121], and their value then appears to fall
outside.) See
[`verifications/superseded/dim44-45-p48-cross-sections/`](../verifications/superseded/dim44-45-p48-cross-sections/).

The natural next input is the Siegel theta series of degree ≥ 3. Ozeki computes degrees 2
and 3 explicitly and degree 4 almost explicitly, and the degree-2 series is now known to be
uniquely determined by extremality. Its Fourier coefficients count *g*-tuples of lattice
vectors with a prescribed Gram matrix — exactly the quantities bounded here — so they would
replace these bounds by exact values in dimensions 44–47 and 68–71. Ozeki's Theorem 6.2
shows every degree-4 coefficient has the form panθ(*T*) + *t*·*J*³(*T*) with one global
unknown *t*, bounded only within a factor 278/184 by his Theorem 6.3, so the degree-4 counts
are in general lattice-dependent and the linear program is doing real work at *k* = 3.
