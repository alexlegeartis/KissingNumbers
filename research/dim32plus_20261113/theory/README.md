# Gamma72 head compression: 13 November 2026

## Certified result

The research certificate proves **tau(69) >= 743,522,304**, an increase of
**54,069,382 (7.84%)** over the preceding repository certificate, 689,452,922.
The final maintained verification package is assembled separately by the main
agent; this directory preserves the discovery route and the first exact proof.

```
python research/dim32plus_20261113/theory/verify_head69.py
```

The verifier uses Python's standard library only. It reconstructs 1,681 possible
D4 profiles, their sixteen averaged/translation groups, ten exact moment rows,
the embedded normal Gram matrix, a rational population dual, and all 55 geometric
pair types. Its exact dual objective is 743,522,304, with no rounding loss.

The construction starts with Gamma72's norm-eight shell. With four actual
minimum-vector normals U and Gram G, write c=Uv and

    b(v) = v - U^T G^(-1)c,  q(c,d) = c^T G^(-1)d.

For each of ten selected profiles c, replace v by

    (alpha_c b(v), h_c) in R^68 x R,
    alpha_c^2 = (8-h_c^2)/(8-q(c,c)).

Thus each image has squared norm eight. Positive alpha values and the old
inner-product bound four give the sufficient cross-profile inequality

    alpha_c alpha_d (4-q(c,d)) + h_c h_d <= 4.

All heights and alpha-squares in `head69_certificate.json` are rational. The
checker squares this inequality only after checking the necessary signs.
The same inequalities exclude image collisions. Averaging over orientations of
the embedded D4 normal lattice gives at least one orientation attaining the
certified population; no invariance of the ambient Gamma72 shell is assumed.

The published inputs remain Gamma72's existence, integrality, minimum eight,
6,218,175,600 minimal vectors, its spherical 11-design property, and the catalogue
Gram identity. This is an exact existence construction and an algebraic recipe,
not an explicit enumeration of hundreds of millions of coordinates. It proves
a lower bound for the unrestricted kissing number; it makes no lattice-packing
claim.

## Discovery and bounded follow-ups

| Search | Discrete vertices | Best numerical population lower bound | Compared repository bound | Outcome |
|---|---:|---:|---:|---|
| D4, one new head coordinate, 69D | 16,849 | 743,522,304 | 689,452,922 | Exact proof above |
| A3, one head coordinate, 70D | 4,861 | 1,331,478,528 | 1,360,986,112 | No improvement |
| A2, one head coordinate, 71D | 1,135 | 2,603,642,568 | 2,603,658,750 | Truncated central-section recovery |
| D5, one head coordinate, 68D | 29,269 | 358,397,339 approximately | 471,884,288 | No improvement |
| D4, two head coordinates, 70D | 20,329 | 1,248,262,038 | 1,360,986,112 | No improvement |
| D4, two heads, refined radii and zero profile excluded | 26,376 | 1,268,662,272 | 1,360,986,112 | No improvement |

The search scores use an interior feasible moment distribution to rank candidates;
the last column of the numerical comparisons uses a separate minimization for
each selected objective. The scores themselves are not lower-bound certificates.

The final 70D attempt jointly optimized all head coordinates with SLSQP while
adding each of seventeen unused norm-quarter profiles and optionally removing
one or two low-weight norm-one profiles. All 51 starts retained violations; the
best cross-profile margin was approximately -0.06614 in norm-eight units.
These failures describe bounded searches and do not prove an obstruction or
optimality. Details are in `continuous70_result.json`.

The C++ discovery engines use a fixed random seed but a wall-clock stopping
condition. Exact replay of the certified result uses the saved certificate,
not a rerun of the heuristic search. Rebuild the engines with a C++17 compiler
before running `run_head_search.py`, `run_dk_head_search.py`, or
`run_twohead_search.py`. Generated vertex lists and binaries are reproducible
intermediates and are not tracked.

## Source audit

The current repository snapshot is substantially stronger than the older
comparison values in recent papers. In particular, Lin's 32--39D results and
Kamenetsky's 49--61D constructions must not be used as the immediate baseline
for a new repository gain. `source_audit.json` records the reviewed primary
sources and the relevant current values.

- [Lin, arXiv:2611.03682](https://arxiv.org/html/2611.03682v1): constant-weight
  codes and the 38D Leech equator; the best stated 35D support size is 2,202.
- [Kamenetsky, arXiv:2611.01906](https://arxiv.org/html/2611.01906v1): the
  7,193-line P48 class and optimized image coverage; the repository already
  has stronger classes, refills and 61D code constructions.
- [PackingStar](https://github.com/CDM1619/PackingStar): the November 11 release
  improves dimensions 26--31 and adds no 32+ record in its current README.
- [Nebe, Gamma72](https://arxiv.org/abs/1008.2862) and the
  [catalogue entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html)
  supply the published parent lattice.

The live Cohn table returned HTTP 502 through the retrieval tool. Search exposed
an older cached table, which was not treated as current. No independent 69D
construction exceeding the repository's preceding value was found in the
targeted primary-source review. This bounded review is not a proof of worldwide
priority.
