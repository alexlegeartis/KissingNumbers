# CLOSED — dimensions 32 to 44 are exactly at the Edel–Rains–Sloane value

> **Correction, 2026-09-14 — this package's verdict is stale in dimensions 40, 41, 42 and 45.**
> I. Dorofeev, X. Sun and C. Wang, *Optimal Extensions of Cross-Sections: Sphere Packings in
> Dimensions 38 to 43*, arXiv:2607.20359v4 (29 Aug 2026), beat the Edel–Rains–Sloane values
> with **lattice cross-sections** of an extremal even unimodular 48-dimensional lattice:
> **1 092 000** in dimension 40 (this package said 1 064 368), **1 324 472** in 41,
> **1 792 386** in 42, and **7 379 838** in 45. Their mechanism takes a saturated rank-$k$
> cross-section $K\subset\Lambda$ and counts $\#\{w\in\Lambda: |w|^2=6,\ w\perp K\}$.
>
> **Strengthened, 2026-09-14 — dimension 32 is now known to be a *maximal* spherical
> code.** Scaled by $\sqrt{32}$ the three layers are integer vectors of norm 32 with pairwise
> inner product $\le 16$, the same framework as the Leech. The pair and octad layers cut out a
> polytope whose maximum Euclidean norm is exactly 1, attained only at the flat vector, so every
> addable unit vector must be a flat $\pm 1/\sqrt{32}$ one; and the $[32,17,8]$ cube code has
> **covering radius 6**, so no flat vector is addable either. Nothing at all can be added.
> The two coding numbers are independent levers, and the neglected one is $A(32,8)=2^{17}$,
> which is 35-optimal (every outside word has $\ge 36$ codewords within distance 7) and linearly
> exhausted ($[32,18]$ has $d\le 7$). See KNOWLEDGE.md section 132.
>
> **Dimensions 32–37 are unaffected**, and now for a measured reason rather than an untested
> one: the kissing number of a bare cross-section tracks
> $R\cdot\mathrm{vol}(B_n(\sqrt6))/\sqrt{\det K}$ with $R=6.92$ — calibrated on their own three
> results, reproducing dimension 40 to 3 % — while their Proposition 7 floors
> $\det K\ge(6/\gamma_k)^k$. At $k=16$ that floor is 168 151, capping dimension 32 at about
> **205 000**, short of its record 346 432 by a factor **1.69**; the rank-16 section that
> actually exists, $\sqrt3(E_8\perp E_8)$ with $\det 3^{16}$, gives only about **12 800**. The
> crossover sits at $k\approx8$–$10$, exactly where their chain stops.
>
> **Dimension 32 is also now known to be rigid.** Its record is a *free product* of three
> separately-maximal objects — the cube $(\pm1)^{32}/\sqrt{32}$ on an $A(32,8)$ code, a flat
> layer on each octad of an $A(32,8,8)$ family, and all $4\binom{32}{2}$ two-sparse vectors —
> whose cross-conditions are satisfied automatically, two of the three with equality. No fourth
> layer exists. Neither coding number has a published upper bound, so both were computed here:
> $A(32,8,8)\le 3540$ (Johnson, correctly nested -- it returns $A(24,8,8)=759$ exactly; the
> Delsarte LP in $J(32,8)$ gives 3596, which is just the counting bound) and $A(32,8)\le 351\,908$. Hence
> $\tau_{32}\le 586\,176$ with $A(32,8,8)$ at its cap and $807\,012$ with both — against Cohn's
> SDP ceiling of 3 162 316. Full account: `KNOWLEDGE.md` §130–131.


**Nothing is claimed here**, but this is the most immediately *useful* package in the
`closed/` tree, because it identifies four dimensions where **Cohn's table is stale** and
says exactly where any future code improvement will pay off.

```
python PIPELINE.py          # ~1 min, pure integer arithmetic, no third-party packages
```

## What was done

The Edel–Rains–Sloane construction

> τ_n ≥ N(n) = Σ_ν A(n, n_ν, n_ν) · A(n_ν, ⌈n_ν/4⌉)

over a chain n ≥ n₀ ≥ 4n₁ ≥ 16n₂ ≥ … was **maximised over all admissible chains**, fed with
every currently published lower bound on A(n,d,w) and A(m,d) (Brouwer's `Andw.html`, fetched
2026-08-21 and citing Echols' emails through 2026-08-08; Grassl's codetables.de; the
Litsyn–Rains–Sloane table).

The formula is confirmed three ways before it is used: it reproduces the two large published
ERS totals to the digit (n = 64 and n = 80), it reproduces every entry of Cohn's table in
dimensions 32–44 where the code tables have not moved, and its geometry is re-derived from
scratch as exact inequalities.

## The result

| n | best chain | N(n) | Cohn's table | best known | verdict |
|---|---|---|---|---|---|
| 32 | (32,8,2) | **346 432** | 345 408 | 346 432 | at the value — **table stale** |
| 33 | (32,8,2) | **362 048** | 360 640 | 362 048 | at the value — **table stale** |
| 34 | (32,8,2) | **381 124** | 380 868 | 381 124 | at the value — **table stale** |
| 35, 36, 40 | (32/40,8,2) | = table | = | = | at the value |
| 37 | (32,8,2) | **496 232** | 494 312 | 496 232 | at the value — **table stale** |
| 38 | (38,8,2) | 570 236 | 566 652 | **591 612** (this project) | below |
| 39 | (39,8,2) | **756 116** | 755 988 | 756 116 (this project) | at the value |
| 41, 42 | **(40,8,2)** | = table | = | = | at the value |
| 43 | (43,8,2) | 1 745 692 | 2 060 399 | 2 060 399 (Sun–Wang) | below |
| 44 | (44,8,2) | 2 948 552 | 2 948 552 | **4 153 600** (Sun–Wang) | below |

## Three things worth keeping

**Cohn's table is stale in four of these dimensions** — 32, 33, 34 and 37 — from Echols'
2026 constant-weight improvements (arXiv:2608.13906), which Brouwer's page already carries.
In particular **dimension 37's live record is 496 232, not 494 312.** None of this is a claim
of this project; it is a correction to the table, and it is worth passing on.

**The chain choice is non-obvious.** Dimensions 41 and 42 optimise at n₀ = 40, *not* n₀ = n,
because A(41,11) = 2¹⁸ and A(42,11) = 2¹⁹ both sit **below** A(40,10) = 589 824. Anyone
evaluating this construction should enumerate the chains rather than assume n₀ = n. That
same observation, applied at n₀ = n, is what gives dimensions 38 and 39 — see
[`../../improved/dim39-ers-constant-weight/`](../../improved/dim39-ers-constant-weight/).

**Where the leverage is.** Every dimension except 38, 43 and 44 sits *exactly* at the ERS
value, so **any** code improvement is immediately a record:

- +1 in A(n,8,8) gives **+128** in τ_n;
- +1 in A(32,8) gives +1 in τ₃₂ … τ₃₇ **simultaneously**;
- +1 in A(40,10) gives +1 in τ₄₀, τ₄₁, τ₄₂ **simultaneously**.

The middle levels are hopeless: weight 9 is killed outright by the Johnson bound, and for
weights 10 and 11 nothing at all is published above n = 32 — beating the weight-8 level at
n = 44 would need A(44,10,10) ≥ 11 773, about 24 times the largest published value.
codetables.de gives d(42,20) = d(43,21) = d(44,22) = 10 exactly, so no linear code can lift
the top-level sign codes either.

**The one place ERS leaves slack, and why it does not pay.** ERS forces two level-0 supports
to overlap in ≤ n₀/2, which is what makes the sign patterns of different supports
unconstrained and gives the clean product A(n,n₀,n₀)·A(n₀,⌈n₀/4⌉). For n ≤ 44 and n₀ ≥ 30 two
n₀-subsets of [n] meet in ≥ 2n₀ − n > n₀/2, so A(n,n₀,n₀) = 1 and level 0 is a *single*
support. One could instead take several supports with larger overlap and constrain the signs
across them. Worked out in the only case where it could matter, n = 44 with n₀ = 40: a
589 824-word code projected to the 36 shared coordinates **saturates the cube** — a Hamming
ball of radius 7 in {0,1}³⁶ holds 10 739 176 points, and 589 824 · 10 739 176 / 2³⁶ = 92, so
an average point already has ~92 first-support words within distance 7. Essentially no room
is left for a second support, and it would take 3.56 supports' worth merely to match the
single-support n₀ = 44 term already in use. At n = 44 the winning chain has n₀ = n anyway, so
there is no freedom of this kind at all.

## Files

```
PIPELINE.py             the whole audit; pure integer arithmetic, self-checking
run_2026-08-21.out      its output on the date the code tables were fetched
README.md               this file
data/Andw.html          Brouwer's constant-weight table, as fetched
data/andw_parsed.json   parsed A(n,d,w) values
data/bklc_*.html        Grassl's codetables.de pages for the sign codes
data/bklc_parsed.json   parsed
data/LRS_tableand.html  Litsyn-Rains-Sloane
data/lrs_parsed.json    parsed
data/echols.html        Echols arXiv:2608.13906 in HTML, used to audit Andw.html
data/echols_table1.txt  its Table 1, extracted -- the 124 new constant-weight codes
data/binary-1.html      Brouwer's A(n,d) table for n <= 28
data/cohn_kissing.html  Cohn's table, as fetched 2026-08-20 -- the provenance for every
                        "previously published" figure in this repository
data/agrell_unr_*.html  Erik Agrell's unrestricted-code table, https://codes.se/bounds/unr.html
data/cohn_upper.json    Cohn's table as [lower, upper] pairs, 49 entries -- the same
                        1-48 and 72 covered elsewhere here.  The first of each pair
                        agrees with common/published.py's COHN in all 49; the second
                        is 1 in exactly the six dimensions where tau is known exactly
                        (1, 2, 3, 4, 8, 24) and the upper bound in the rest.  Nothing
                        in this repository reads it -- every claim here is a LOWER
                        bound -- so it is provenance for the other half of the table,
                        not an input.
scripts/parse_andw.py   Brouwer's constant-weight table  -> andw_parsed.json
scripts/parse_bklc.py   Grassl's codetables.de pages     -> bklc_parsed.json
scripts/parse_lrs.py    Litsyn-Rains-Sloane              -> lrs_parsed.json
```

`PIPELINE.py` reads only `binary-1.html`, `LRS_tableand.html`, `Andw.html` and the three
`*_parsed.json` files. The rest are **archived evidence rather than inputs**: `echols*.*`
is what `Andw.html` was audited against, `cohn_kissing.html` is the fetched table itself,
and the two `agrell_unr_*.html` are the source for the statement in the header that Agrell's
table stops before the range that matters here. None of them is loaded, and none needs to be.

The HTML pages are shipped deliberately: **the code tables move**, and a year from now the
only way to reproduce this audit will be to know exactly what they said when it was run.

## References

See [`../../../common/published.py`](../../../common/published.py).
