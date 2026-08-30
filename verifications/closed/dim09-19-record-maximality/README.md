# CLOSED — every published record in dimensions 9 through 19 is MAXIMAL

**Nothing is claimed here.** This records that the +2 trick of dimension 25 — find a hole
and drop a sphere into it — has **no analogue anywhere in dimensions 9 through 19**, with
exact margins showing how far each record is from admitting one.

```
python sweep.py             # needs numpy and Cohn's dimensions1-24.txt (see below)
```

`sweep.py` ends with a verdict and an exit code. It requires every record configuration it
reaches to be maximal and its m(C) to match the table below to 1e-6, and exits 1 otherwise,
so a regression in the ascent or in the parser cannot pass unnoticed. Without Cohn's file it
prints SKIP and exits 0. It is part of `run_all.py`'s fast set.

## The test

For a configuration C of unit vectors, a further point can be added iff

> m(C) := min over unit *z* of max over *x* ∈ C of |⟨z, x̂⟩| ≤ 1/2.

m(C) is computed by maximising |z| over the polytope {z : |⟨z,x̂⟩| ≤ 1/2}, by cutting-plane
ascent with an explicit feasibility assertion on every returned point. Cohn's
`dimensions1-24.txt` gives explicit coordinates for every record, so the whole low range can
be swept at once.

> **Method warning, worth keeping.** A naive iterated LP over all constraints is
> **unreliable** for this. On a control case (Λ₂₂ minus its two outer layers, where the axis
> is provably free with max|z| = √(8/3)) it returned an answer 70% low. It has to be a
> cutting-plane ascent with feasibility asserted each time.

## The result: none of them admits an extra sphere

`python sweep.py` covers dimensions **9, 10, 12, 13, 14, 15, 16** directly from Cohn's
coordinates. Dimensions 11, 17, 18 and 19 are covered by separate runs, noted below.

| dim | record | source | m(C) | exact | margin over 1/2 | swept by |
|---|---|---|---|---|---|---|
| 9 | 306 | | 0.577350269 | **1/√3** | 15.5% | `sweep.py` |
| 10 | 510 | Ganzhinov 2025 | 0.670820393 | | 34.2% | `sweep.py` |
| 11 | 604 | Bianchi et al. 2026 | 0.577350269 | **1/√3** | 15.5% | separately — see below |
| 12 | 841 | Takhanov et al. 2026 | 0.604964075 | | 21.0% | `sweep.py` |
| **13** | **1154** | Zinoviev–Ericson **1999** | 0.554700196 | **2/√13** | **10.9%** | `sweep.py` |
| **14** | **1932** | Ganzhinov 2025 | 0.534522484 | **√(2/7)** | **6.9%** | `sweep.py` |
| **15** | **2564** | Leech–Sloane **1971** | 0.534522484 | **√(2/7)** | **6.9%** | `sweep.py` |
| 16 | 4320 | Barnes–Wall **1959** | 0.577350269 | **1/√3** | 15.5% | `sweep.py` |
| 17 | 5730 | Cohn–Li 2024 | 0.522232968 | | 4.4% | separately — KNOWLEDGE §67a |
| 18 | 7654 | Cohn–Li 2024 | 0.530330086 | | 6.1% | separately — KNOWLEDGE §67a |
| 19 | 11948 | Ho 2026 | 0.521948010 | | 4.4% | separately — KNOWLEDGE §67a |

All of them have maximum inner product exactly 0.5, so they are valid 60° codes, and **none
admits a free point**.

**Dimension 11 is skipped by `sweep.py`** — it reports `parsed 502 of 604 -- SKIP`, because
the dimension-11 block of Cohn's file uses coordinate tokens the parser in `extract.py` does
not handle, and the script correctly refuses to test a partial configuration rather than
report a wrong answer on 502 of 604 points. Its value in the table above comes from a
separate run over the 604-point configuration reconstructed in the working tree. If you want
that row reproduced from Cohn's file, extending `extract.py`'s token grammar is the job.

**Cohn's file lists several configurations per dimension**, and `sweep.py` tests every one
— which is why the raw output repeats a dimension several times. A free point in *any* of
them would be a record, so none can be skipped.

**Dimensions 14 and 15 share the same minimax exactly**, which points at a common structure,
and they are the tightest in the range — their deep holes are only 6.9% too shallow. So
*they*, not 11 or 12, are the natural surgery targets if anyone wants one.

Note how old some of these are: dimension 16 is Barnes–Wall **1959**, dimension 15 is
Leech–Sloane **1971**, and dimension 13 is Zinoviev–Ericson **1999** — the last surviving
1999 record in the low range.

This is **maximality, not optimality**: a different configuration could still be larger, but
not one built by adding to these.

## A second dead end, measured: layering over a record equator

The low-dimensional records are layered, N = |E| + 2|L| + 2, with the caps at height 2
obeying ⟨x̂,ŷ⟩ ≤ 1/√3 against the equator and ⟨ŷ,ŷ′⟩ ≤ 1/3 among themselves. In both
decoded cases the equator sits well **below** the record of the dimension underneath:

```
dim 13: equator 816   against dimension 12's record 841   (25 short)
dim 15: equator 1484  against dimension 14's record 1932  (448 short)
```

so the obvious move is to substitute the record. The cap region is non-empty exactly when
m(E) ≤ 1/√3 = 0.5774, and both records qualify (0.5547 and 0.5345) — **but the layer
collapses**:

| target | equator used | cap layer found | cap layer needed | result |
|---|---|---|---|---|
| dim 14 | dimension-13 record, 1154 | **23** | 389 | 1202 against 1932 |
| dim 15 | dimension-14 record, 1932 | **28** | 316 | 1990 against 2564 |

An order of magnitude short, and the reason is exactly the maximality measured above: a
record configuration is maximal with a small margin, so its 1/√3-free region is a thin
sliver that holds only a couple of dozen directions at pairwise 1/3. **The published
constructions use a deliberately sub-maximal equator — 816 rather than 841 — precisely to
keep the cap region big, and that trade is already optimised.**

> **The rule this yields.** In a layered kissing construction the equator should be chosen to
> maximise |E| + 2|L(E)|, and |L(E)| falls off a cliff as m(E) approaches 1/√3. Substituting
> a record equator is always a mistake.

## History note: two records landed mid-project

| dim | was | became | when |
|---|---|---|---|
| 11 | 592 (Ganzhinov) → 593 (AlphaEvolve) | **604** | Bianchi et al. 2026, arXiv:2606.10402 |
| 12 | 840 | **841** | Takhanov–Assylbekov–Yun 2026, arXiv:2606.18984 |

Both were pursued in this project before the papers were found. Dimension 11's progression
is 582 (Best 1980 / Zinoviev–Ericson 1999) → 592 → 593 → 604; dimension 12's two-block
ansatz is capped at **840**, so the 841st sphere genuinely needs the different structure
Takhanov et al. found. **Cohn's table moves — re-fetch before claiming anything.**

## Files

```
README.md         this file
sweep.py          the cutting-plane ascent for m(C), over dimensions 9-16
extract.py        parses Cohn's coordinate file into configurations
maximal841.py     the same test on the dimension-12 record, where it was first run
data/c841.npy     Takhanov-Assylbekov-Yun's 841 points in R^12, unit vectors
```

`extract.py` caches each configuration it parses as `c<dim>.npy` beside itself; those are
regenerable and are gitignored. `data/c841.npy` is not a cache -- it is the dimension-12
record itself, and `maximal841.py` runs from it directly.

`sweep.py` needs H. Cohn's coordinate file `dimensions1-24.txt` from
<https://hdl.handle.net/1721.1/153312> — 24 MB, not shipped, because it is not ours to
redistribute. It is found automatically if it sits anywhere at or above this repository;
otherwise pass its path or set `D124`:

```bash
python sweep.py                                  # finds it if it is nearby
python sweep.py /path/to/dimensions1-24.txt      # or say where it is
python sweep.py 13 14 15                         # only these dimensions
```

Without it the script exits with instructions rather than a traceback.

## References

See [`../../../common/published.py`](../../../common/published.py).
