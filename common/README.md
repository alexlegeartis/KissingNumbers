# `common/` — shared modules

Six modules, two validation suites and one proof document, used by more than one
package. Each is runnable on its
own and self-validating: running it exercises the code on cases whose answers are known
independently, and exits non-zero if any of them fails.

```bash
python theta.py         # E_7 = 126, Lambda_23 = 93150, dim 47, dim 71   (seconds)
python kpoint_lp.py     # the LP, TIGHT on the Leech through k = 4, and on P_48  (13 s)
python m24.py           # automorphisms of the Golay code                (seconds)
python published.py     # prints the published table with its provenance

python validate_lp.py   # ADVERSARIAL: 15 known truths inside their LP brackets  (~4 min)
python capalgebra.py    # the cap construction's algebra and count, exactly      (seconds)
```

The last two are the ones to run if you doubt the mathematics rather than the bookkeeping.
Between them they attack the two things every one of the 47 claims rests on.

## The modules

### [`theta.py`](theta.py) — extremal theta series and the Venkov one-point distribution

`extremal_theta(n)` computes, by exact linear algebra over ℚ from E₄ and Δ, the unique
weight-n/2 modular form for SL₂(ℤ) with constant term 1 whose first ⌊n/24⌋ coefficients
vanish. That determines (minimum, kissing number) from the dimension alone: (2, 240),
(4, 196 560), (6, 52 416 000), (8, 6 218 175 600). **The lattice never has to be named**,
which is why the dimension 46/47 statements hold for all four of P₄₈p, P₄₈q, P₄₈m, P₄₈n.

`onepoint(n, mu, N, t)` solves the Venkov moment system for the one-point distribution. It is
**over-determined** — μ/2 + 1 unknowns against t/2 + 1 equations — so the surplus equations
are a genuine consistency check, and the function returns whether they hold. Its self-test
reproduces E₇ = 126 and Λ₂₃ = 93 150, both known independently, and then gives dimensions 47
and 71.

### [`kpoint_lp.py`](kpoint_lp.py) — the k-point moment LP with an exact integer certificate

`run3(n, N, mu, IPS, M1, G, ...)` is the engine behind dimensions 46, 47, **68, 69**, 70 and 71
and behind the withdrawn 44 and 45. It builds the linear program described in
[`PROOF-kpoint.md`](PROOF-kpoint.md) — design moments (**including the mixed ones**), exact
one-point marginals, antipodal symmetry, lattice translations, lattice-range constraints,
unit upper bounds on singular cells — solves it in floating point *only to find a dual vector*,
then proves a bound from that vector in exact integer arithmetic.

**It is the third return value, never the second, that is ever claimed.** The second is the
floating-point LP value; the third is the certified integer. An imprecise solver can only make
the certificate weaker, never invalid.

Three things about the certificate stage are worth knowing before touching it.

*It is integer, not rational.* Every constraint is scaled by the lcm of its denominators and
divided by its gcd, and the dual is carried as one integer vector over a **single** common
denominator `D`, so `D·slack_j = D·c_j − Σᵢ aᵢ Aᵢⱼ` is an integer: feasibility is an integer
comparison and the bound is one division. The previous version carried the dual as unrelated
`Fraction`s from `limit_denominator`; against row denominators up to `μ^tdes` the running
denominator of one inner product passed 10²⁸, and a single feasibility check on the k = 4 system
took a quarter of an hour. Same mathematics, 33× faster.

*There is more than one dual, and the best is kept.* The certificate validates whatever dual it
is handed, so an extra candidate can only raise the bound — never invalidate one. Two sources
are tried: the equality form `A_eq = A`, which is half the rows and has no tolerance to guess,
and the ±EPS relaxed inequality form that was here before. Taking the better of the two is worth
+1 259 008 in dimension 68 by itself. `KN_CERT=rational` selects the old rational tail for
comparison, `KN_BUDGET` sets the wall-clock allowance (default 900 s), and a per-**solve**
`time_limit` is passed to HiGHS as well — without it the A₄ four-point Gram spends tens of
minutes inside one call and the caller has no way out.

*Minimising and maximising round differently.* For `sense = -1` the tight integer bound on the
cell is `floor(−B)`, not `−floor(B)`; the two differ by one whenever `B` is not an integer. That
is never invalid — it is still an upper bound — but it costs the `min = max` test on the P₄₈
two-point cross-section, which needs both ends to meet at 12 309 600 exactly.

The dimension 68–69 package used to carry a byte-identical copy of this file as
`scripts/kgram9.py`. It imports this one now, so there is a single engine, and
`../verifications/improved/dim68-69-gamma72-cross-sections/scripts/recheck_cert.py`
re-derives both of its certificates by a second, independent arithmetic path.

> ⚠ **If you modify this file, test it on a pair of equivalent Gram matrices.** The Grams
> (3,−3,0) and 3·A₃ satisfy UᵀGU = 3A₃ for an integer U of determinant −1 and must give the
> same answer. An earlier revision imposed only exponent vectors with every entry even — 56
> equations at k = 3 instead of 161 — and gave 6 462 480 and 5 194 969 for those two. The bug
> was visible in the output and went unnoticed for a day. The check is cheap and decisive.

### [`published.py`](published.py) — the published state of the art, with provenance

`COHN` is the lower-bound column of Cohn's *Table of kissing number bounds*, re-fetched and
re-read 2026-08-20, with `COHN_SOURCE` naming who holds each entry this repository touches.
`floor_for(d)` returns the best previously published value in dimension *d* **together with
an explanation of where it comes from**, which is what every "previously" figure in the
repository is generated from — none is typed by hand.

The point of this module is the range above 48, where getting it wrong is easy: **Cohn's
table covers dimensions 1–48 plus 72 and nothing else**, and the Nebe–Sloane table it
replaced lists 1–40, 42, 44, 48, 64, 72, 80, 128. The four regimes above 48 — direct sums
over P₄₈, monotone inheritance from dimension 64, direct sums over Γ₇₂, and
Edel–Rains–Sloane itself from dimension 96 — are documented in the module docstring, with
the reverse-engineering of the Edel–Rains–Sloane count in each.

**It does not rule ERS out everywhere, and the docstring now says where it fails.** At
n₀ = 96 the code [96,33,24] makes ERS the record outright, and at n₀ = 62 and 63 the codes
[62,26,16] and [63,27,16] beat this repository's own claims. See
[`../verifications/closed/dim96-ers-takeover/ers_sweep.py`](../verifications/closed/dim96-ers-takeover/ers_sweep.py),
which asks the question in every claimed dimension instead of arguing it in prose.

**Cohn's table moves.** Two entries (dimensions 12 and 19) changed while this work was in
progress. Re-fetch before claiming anything.

### [`validate_lp.py`](validate_lp.py) — adversarial validation of the k-point LP

The LP produces **lower** bounds. A lower bound that is too small is merely weak; one that is
too **large** is wrong, and it becomes too large the moment any single constraint in the
program is not actually valid. This file attacks that failure mode the only way that is
convincing: it runs the program on every configuration whose true size is known
independently and checks

```
LP minimum   <=   the true value   <=   LP maximum
```

**Fifteen data points**, from four lattices and three separate sources of truth — direct
enumeration from coordinates, the classical laminated values, and other people's published
computations (Boyvalenkov–Cherkashin, Ozeki, Sun–Wang). If any constraint were too strong,
some truth would fall below its own minimum. None does. Several are pinned exactly, and one
— the Leech's orthogonal triple — has LP maximum equal to the truth, 19962.

It also runs the **equivalent-Gram regression test** first: (3,−3,0) and 3·A₃ satisfy
UᵀGU = 3A₃ for an integer U of determinant −1, so they describe the same configuration and
must give the same answer. That is the cheap decisive test for the missing-moment class of
bug, and it is the test the README asks people to run — so it is now run automatically.

### [`capalgebra.py`](capalgebra.py) — the cap construction, from first principles

Thirty-nine of the forty-seven claims come from one construction, and all of them rest on two
things:
that a list of inequalities holds at the chosen cap level, and that a count is right. Both
are re-derived here in exact `Fraction` arithmetic, **independently of every script that uses
them** — the seven kinds of pair, each reduced to one inequality, with square roots handled
by squaring so nothing is decided in floating point.

It also derives, rather than asserts, why P₄₈ gets antipodal pairs and the Leech and Γ₇₂ get
zero-sum triples:

| lattice | γ | forces t ≥ | ⇒ ⟨z,z′⟩ ≤ | directions per line |
|---|---|---|---|---|
| Leech | 1/4 | 2/3 | −1/2 | **3** |
| P₄₈ | 1/3 | 3/4 | −1 | 2 |
| Γ₇₂ | 1/4 | 2/3 | −1/2 | **3** |

and reproduces Cohn's table in dimensions 26, 28, 29, 30, 31 exactly from the count alone.

### [`capcheck.py`](capcheck.py) — the cap construction in coordinates, at t = 3/4

The *pair* half of the cap family, built and checked in real coordinates over the Leech
lattice: equator, caps at level t = 3/4 with each class line carried on an antipodal pair of
directions, and poles. It is the analogue of
[`../verifications/improved/dim73-95-gamma72-caps/scripts/calibrate.py`](../verifications/improved/dim73-95-gamma72-caps/scripts/calibrate.py),
which does the *triple* half at t = 2/3, and it supplies `golay()`, `leech()`,
`strong_classes()` and `zlines()` to the packages in `verifications/closed/` that need Leech
coordinates.

Gamma_72's 6 218 175 600 minimal vectors cannot be written down, so the construction used for
dimensions 73–95 is exercised here on its exact analogue over the Leech, where every step of
the argument — the deletion rule, the within-group condition, the cross-group condition and
the poles — is identical and only the numerical values of mu and the class size differ.

### [`leech.py`](leech.py), [`golay.py`](golay.py), [`m24.py`](m24.py)

The extended binary Golay code from the cyclic [23,12,7] QR code, the 196 560 Leech minimal
vectors at squared norm 32 built from it (1104 of shape (±4², 0²²), 97 152 of shape
(±2⁸, 0¹⁶) on octads, 98 304 of shape (∓3, ±1²³)), and M₂₄ generators. `leech.load()` caches
to `data/leech196560.npy` on first use.

Nothing is taken on trust: the Golay weight enumerator is checked to be
1 + 759x⁸ + 2576x¹² + 759x¹⁶ + x²⁴, and the Leech vectors are checked to be 196 560 distinct
vectors of squared norm 32.

### [`data/golay_basis.txt`](data/golay_basis.txt)

Twelve generators of the Golay code, as used by the dimension-25 and dimension-38 packages
and by the calibration script. Those scripts expand it to 4096 words and verify the weight
enumerator before using it.

## The proof

[`PROOF-kpoint.md`](PROOF-kpoint.md) is the mathematics behind dimensions 68, 69, 70 and 71,
and behind the recovered 46 and 47, written out as eight numbered lemmas with proofs:
cross-sections are kissing configurations (Lemma 1), the moment identities (Lemma 2) with the
translation (2b) and lattice-range (2c) constraints, the weak-duality certificate (Lemma 3),
realizability (Lemma 4), the one-point distribution in closed form, what "exact" does and does
not mean, and the independent confirmations from Boyvalenkov–Cherkashin, Ozeki and Sun–Wang.

The cap construction's own accounting is derived in the READMEs of
[`../verifications/improved/dim49-63-p48-caps/`](../verifications/improved/dim49-63-p48-caps/)
and
[`../verifications/improved/dim73-95-gamma72-caps/`](../verifications/improved/dim73-95-gamma72-caps/),
and checked in coordinates against Cohn's table by the latter's `scripts/calibrate.py`.
