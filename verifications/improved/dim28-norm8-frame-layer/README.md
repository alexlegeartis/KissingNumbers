# Dimension 28: a full frame of norm-8 heads at height √2

**τ(28) ≥ 204 896**, against the published **204 520** (Cohn's table; Ma et al. 2025,
arXiv:2511.13391, the PackingStar configurations). An improvement of **+376**. Dimension 28 was the
one dimension of 25–31 that this programme had not moved.

```
python verify28.py          # ~1 min, exact arithmetic in every decision; exits non-zero on failure
```

It reads only `data/` and rebuilds the 98 280 Leech minimal lines from the Golay code
(`lib/golay.py`), checking the inner-product signature 46575 / 47104 / 4600 / 1 before it trusts
them, so nothing depends on the code that produced the configuration. `data/classes.npy` stores
the owner classes as **vectors**, not indices, so no construction order is assumed. Every pair is
decided exactly: integer inner products for the ℝ²⁴ parts, sympy in ℚ(√2, √3) for the ℝ⁴ parts.

## The configuration

Norm-4 units: a point of ℝ²⁸ is `(x, y)` with `y ∈ ℝ⁴`, and two points are compatible iff their
inner product is at most 2.

| | points |
|---|---|
| equator `(u, 0)`, `u` a minimal vector that is not an owner | 192 592 |
| caps `(√(2/3) u, (2/√3) z)`, `u` an owner, `z` a direction of its triangle | 11 904 |
| axis `(0, 2z)`, `z` one of the 16 half-vectors of the dual 24-cell | 16 |
| **layer `(v/2, √2 w)`, `v = ±f₁…±f₂₄` a Leech frame, `w` one of `±e₁..±e₄`** | **384** |
| | **204 896** |

The classical construction is the same equator and caps with an axis of 24 and no layer:
`192 592 + 11 904 + 24 = 204 520`. The layer costs 8 axis points and pays 384.

## Why it works

The head is a **norm-8** Leech vector. Every earlier layer in this programme used a norm-6 head
(`f/√3` in dimensions 28, 30, 31; `v/2` with `|v|² = 6` in dimension 29), and for a norm-6 head the
`m = ⟨f,u⟩ = 2` case never clears D4's covering cosine `1/√2`. For a norm-8 head it does: the
cap–layer constraint `√(2/3)·⟨u,v⟩/2 + (2√2/√3)·⟨z,w⟩ ≤ 2` has binding value

    (√6 + 2√3)/3 = 1.971197119…  <  2        ⟺   18 + 12√2 < 36   ⟺   √2 < 3/2

so only `m ≥ 3` is forbidden — 2 094 lines per head instead of 11 730. Also `⟨v/2, u⟩ ≤ 4/2 = 2`
exactly, so **no equator point is deleted**.

A head carries directions with `⟨w,w'⟩ ≤ 0`, i.e. at most a cross-polytope (8), and two heads share
a direction iff `⟨v,v'⟩ ≤ 0`, so a set of mutually orthogonal head lines and all their negatives
share all 8. With `H` head lines the layer holds `16H` points and the axis drops to 16, giving
`NET = 16H − 8` with `H ≤ 24`, since 48 pairwise non-positive vectors in ℝ²⁴ form a frame.

`H = 24` is reached here. It requires the 8 owner classes to avoid `{u : |⟨u,f_i⟩| ≥ 3}` for every
frame vector, which in frame coordinates says exactly: **each class consists entirely of type-B
lines, shape `(2⁸, 0¹⁶)`**. Such classes exist — the six tetrads of a *sextet* each carry 4
mutually orthogonal `(4⁴)` sign-lines, those 24 lines form a single Leech frame (two tetrads of a
sextet union to an octad, so every `(f_i − f_j)/2` is minimal), and in that frame the record
248-line class is all type-B.

## Files

| | |
|---|---|
| `verify28.py` | the verification |
| `data/classes.npy` | the 8 owner classes, 8 × 248 × 24 integer vectors at norm 32 |
| `data/heads.npy` | the 24 frame heads, norm 64 |
| `lib/golay.py` | the extended binary Golay code |
| `verify28.log` | a recorded run |

The construction scripts live in `research/collab2531/dim28/`; the working note is FINDINGS
section 58.

**Verifier strengthened 2026-09-15.** One check could not fail, two facts the case list used
were asserted rather than checked, and one check was misnamed:

* `heads are lattice vectors` read `... in KEY or True`, which is true for any input — and could
  not have passed as written either, since `KEY` holds minimal vectors and a head has twice the
  norm. It is now the Conway–Sloane membership criterion, against the same Golay code the shell
  is built from, and it rejects 3 of 7 probe vectors.
* `layer x equator, <v,u> = 4` assumed a bound on `max <u,v>` over the **whole** shell. That scan
  now runs (the answer is exactly 32, so the layer touches the equator and deletes nothing).
* `cap x axis, <z,z'> = 1/sqrt2` assumed a cosine between two direction sets that section [4]
  never compared. It does now.
* `class k is 60-degree free` named the wrong condition: what the cap–cap bound needs is
  `|<u,u'>| <= 8`, i.e. at least arccos(1/4) = 75.52°, since the owners of one class share all
  three of its cap directions. The check was right; only its name was wrong.

The pair type `equator x axis` (identically 0) was also missing from the case list. Separately,
an independent rebuild of all 204 896 points — from this package's own data and Golay code, every
pair checked except equator-against-equator, which the shell signature settles — gives a maximum
inner product of exactly 2 with no violation.

## History of the bound in dimension 28

| year | bound | who |
|---|---|---|
| 2025 | 204 520 | Ma et al., arXiv:2511.13391 (PackingStar), the value in Cohn's table; the cap construction over the Leech lattice with its 248-line class reproduces it exactly (`common/capalgebra.py`) |
| 2026-09-15 | **204 896** | this package |
