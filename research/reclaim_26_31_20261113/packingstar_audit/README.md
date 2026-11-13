# PackingStar audit and a 30D Barnes–Wall module replacement

13 November 2026. The new explicit exchange has **248,020 points in 30D**,
**28 more** than PackingStar's 247,992-point November release. It removes
932 old points and inserts 960 new points, with zero deletions outside one
14-dimensional module. The search script checked every changed pair directly
with integer arithmetic. Unchanged pairs inherit the pinned source until its
separate complete replay is performed; this directory does not itself claim
to supply that full-source replay.

The later Leech–Sloane section improves the campaign result to **248,084**.
Both complete final proofs are supplied in the portable
[Barnes–Wall](../../../verifications/improved/dim30-bw14-module/) and
[Leech–Sloane](../../../verifications/improved/dim30-leech-sloane-section/)
packages. The experiments below document the earlier discovery stage.

## Source and the actual search restriction

The reviewed PackingStar commit is
[`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6).
The README labels the current update 11 November 2026. Its 26–31D counts are
201,180; 209,862; 222,104; 231,774; 247,992; and 251,254.

The premise that PackingStar is restricted to rational coordinates is false.
In the latest paper,
[`arXiv:2511.13391v5`](https://arxiv.org/html/2511.13391v5), the separate 13D
experiment explicitly restricts **Gram entries** to rational cosine sets.
The general algorithm restricts products against its basis to a finite set
`C1` (Methods, equation 18); the remaining product constraint `C2` may be
continuous (equation 19). It can also impose a structural constraint `C*`,
such as membership in a Leech shell. Earlier paper Table 3 already includes
irrational pairwise cosines in its 25–31D constructions. Moreover, the current
26D coordinate archive contains 16 points over `Q(sqrt(3))`, as recorded by
the [preceding exact source audit](../../reclaim_26_31_20261112/source_and_layers/).

Changing coordinates by an irrational rotation leaves the Gram matrix
unchanged and creates no new kissing geometry. The useful opening is to
change the admitted substructures or Gram entries. This successful 30D
replacement itself uses integer coordinates after the usual norm-32 scaling.

The exact current 30D source is
`PackingStar/25D-31D/new_bounds_20261111/D30_247992_coordinates.npz`, SHA256
`5e46243e5f3cf4790d56abc4008c9cebeb45bdde2fd3080e44a0546d59cccc5a`.

## Construction

Work with vectors of squared norm 32 and pairwise product at most 16. Set

```text
I = {0,2,4,5,6,10,11,23}
```

using zero-based indices among the 24 Leech coordinates. The six tail
coordinates are 24 through 29. The incumbent module supported on these
fourteen coordinates has

```text
240 E8 equator points + 60 D6 pure-tail points + 932 mixed heads = 1232.
```

The fourteen-coordinate section of the standard Barnes–Wall configuration
in sixteen dimensions contains 1,260 points. Generate all 364 coordinate
roots `(±4,±4,0,...)`, then the seven affine hyperplanes of `F2^4` that avoid
coordinate labels 14 and 15. Each hyperplane supplies 128 even sign patterns
with eight entries `±2`. Thus

```text
364 + 7 × 128 = 1260 = 240 + 60 + 960.
```

The script independently checks all pairs in this small module. Its 960
mixed heads have core and tail squared norms both equal to 16. Rotate their
eight-coordinate core by the exact matrix `R_core/4`, where the rows of
`R_core` are rows `[0,2,1,5,7,6,3,4]` of `J−4I8`. Thus this is an E8 root
reflection combined with a coordinate permutation. The tail map has numerator

```text
 0 -4  0  0  0  0
 0  0  0  0 -4  0
 0  0  0  0  0 -4
 0  0  0 -4  0  0
 4  0  0  0  0  0
 0  0  4  0  0  0
```

and denominator four. The machine-readable module includes the original
heads and the complete transformation, satisfying `R R^T = 16 I14` and
`base_heads R^T = 4 heads` exactly.

The resulting 960 heads are compatible with every source point outside the
932 replaced heads. Consequently

```text
247992 − 932 + 960 = 248020.
```

`bw14_module.npz` stores the module and exact transformation.
`exchange30_gain28.npz` stores all 932 removed and all 960 inserted vectors.
`bw14_search_report.json` records the deterministic search, which first found
zero external blockers at iteration 11,025. Its original run continued to
200,000 iterations; the current reproducer stops immediately at zero because
the full-module objective cannot improve further.

## Follow-up computations and their limits

The obvious fourteen extensions of this orientation to a BW15 module in
31D were tested: seven choices for the omitted tail axis, with both signs.
Every admitted two-copy union had gain zero, and the external deletion costs
were 156–178. This is recorded in `bw14_extensions31_report.json`; it does not
exclude other orientations or combinations of modules.

All 759 Golay octads were scanned in the improved 248,020-point source.
The successful normal octad has 152 potentially obstructing external points.
The next two octads have 8,256 each, with 288 incumbent mixed heads and hence
a nominal full replacement gain of 672. A bounded 5,000-step search on
`{3,7,8,9,17,18,20,21}` attained external cost 4,448 and no gain in the tested
unions. Thus directly repeating the module on another octad is not supported
by this pilot. No global exclusion is claimed.

## A precise ceiling for complete fibers in the current tail alphabet

The successful module reaches a restricted ceiling. This explains which
assumption must change before another full-fiber search can improve it.

The complete integer tail alphabet of squared norm 16 admitted here contains
12 axes `(±4,0^5)` and 240 signed tetrads `(±2)^4` on four of six coordinates.
Among the tetrads, join two vertices when their product exceeds eight. This
240-vertex graph is 16-regular. `verify_tail60_ceiling.py` checks entirely
with integers that its adjacency matrix is annihilated by

```math
\prod_{r\in\{-4,-2,0,2,4,8,12,16\}}(A-rI)=0.
```

Because A is symmetric, its least eigenvalue is at least −4. If an independent
set has m vertices, decomposing its indicator along the all-ones vector gives

```math
0\geq \frac{16m^2}{240}-4\left(m-\frac{m^2}{240}\right),
\qquad m\leq48.
```

The 12 axes can all be added, so this tail alphabet admits at most 60 vectors
with pairwise products at most eight. Three tetrad supports complementary to
the pairs of a matching, together with all their sign patterns and the axes,
attain 60. The graph certificate is a finite exact computation, with explicit
integer-overflow bounds; numerical eigenvalues are not used as proof.

At one fixed tail, any norm-16 core bodies must have mutual products at most
zero. There can be at most 16 in eight dimensions, and equality forces a
cross-polytope. An elementary proof normalizes their Gram matrix G: all
off-diagonal entries lie in [-1,0], hence
`sum_ij (G_ij^2+G_ij) <= 2m`, while the same sum is at least `m^2/8` by the
frame-matrix trace inequality and nonnegativity of the squared vector sum.
For m=16, every off-diagonal entry is therefore zero or minus one.

Two complete 16-point fibers cannot have tails of product 12. For any unit
vector in the first cross-polytope, its squared products with the eight
orthonormal axes of the second sum to one. One absolute product is at least
`1/sqrt(8)`, and antipodality makes the positive sign available. The original
norm-16 body product is at least `16/sqrt(8)>4`, violating the required
`core_product <= 16−12=4`. Thus the active tails of complete fibers must form
the at-most-60-point tail code above, and their heads total at most
`60×16=960`. The new module attains this bound.

This ceiling does **not** cover partial fibers, different tail vectors,
different radii, or changes to the surrounding configuration. In particular,
it is not a universal 14D or 30D kissing-number upper bound.

A complementary complete scan priced all `2160×252=544,320` integer
norm-16/core, norm-16/tail heads against the new 248,020-point configuration.
There were no free insertions. Exactly 49,772 new heads each have four old
blockers; all four are among the new 960 module heads. Every one of these
four-blocker sets is distinct. The scan leaves their joint selection open.
The report is `module_head_completion30_report.json`; the catalogue generator
and all omitted-layer guards are in `module_head_completion30.py`.

### Executed individual-fiber searches

`partial_fiber_search30.py` solved 18 bounded local MILPs. Six models used
four incumbent tails (64 potentially deleted heads, 1,062–1,098 candidates);
twelve used six tails (96 heads, 1,948–2,188 candidates). All shared deletions
and candidate/candidate conflicts were represented. Exact fixed-tail
orthoplex-capacity inequalities were also included. Every run reached its
time budget without a positive incumbent; the remaining numerical bounds
were loose or unavailable. These are unresolved searches, not exclusions.

The complete 49,772-point pool was also searched with individual additions,
candidate conflicts computed directly, and shared-deletion reference counts.
A 36-second continuation run made 5,348,646 proposals and reached 871
compatible new heads requiring 915 old deletions; these population counts
are trajectory diagnostics, not an improvement. A separate ordinary
temperature schedule made 260,802,367 proposals in 36 seconds, reaching at
most 45 simultaneous new heads. Neither schedule found a positive exchange;
the best nonempty gain was minus three. Their small reports and C++
implementation are retained. The search does not impose antipodality or
complete fibers. Its negative result supplies no upper bound for this pool.

The C++ input is a little-endian binary stream: three int64 values
`(candidate_count, victim_count, 14)`, all candidate coordinates as int16,
then four int32 victim indices per candidate. Victim indices are obtained
by sorting the distinct `blockers` values in the priced NPZ pool. Its command
arguments are `input.bin seconds rng_seed output.json mode`, with mode
`continuation` or `anneal`. It evaluates new/new conflicts using exact integer
products. The stored searches use seeds 20261113 and 20261114 respectively.

## Reproduction

Run from the repository root, leaving a CPU free and setting an appropriate
thread/affinity budget. The original executor used CPUs 6–7 with two BLAS
threads.

```bash
taskset -c 6-7 env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
  python research/reclaim_26_31_20261113/packingstar_audit/bw14_search.py \
  ../PackingStar/25D-31D/new_bounds_20261111/D30_247992_coordinates.npz \
  --out /tmp/bw14-replay --iterations 200000 --seconds 180

taskset -c 6-7 env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
  python research/reclaim_26_31_20261113/packingstar_audit/extend_bw14_to31.py \
  ../PackingStar/25D-31D/new_bounds_20261111/D31_251254_coordinates.npz \
  --out /tmp/bw14-extensions31-replay

taskset -c 6-7 env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
  python research/reclaim_26_31_20261113/packingstar_audit/verify_tail60_ceiling.py
```

For the positive result the complete direct changed-pair check is decisive;
the search objective and the finite matching calculations are not substituted
for it. The external-source full-pair replay remains a separately stated
verification dependency.
