# Verified midpoint constructions in dimensions 29 and 30

Retained construction: the dated comparison and investigation results remain historical; the dated columns retain their original values; the 28 October column records the [asymmetric exchanges](../../dim26-30-asymmetric-exchanges/) and the [31D module exchange](../../dim31-normal-module-exchange/).

Date: 2026-10-25. These are lower bounds on kissing numbers, not exact kissing
numbers or upper bounds. The comparison targets were checked against the
[Felpix-Studios table](https://github.com/Felpix-Studios/kissing-number), including
its current 241,079 entry in dimension 31.

| Dimension | Published comparison | New verified construction | Improvement | Repo bound (26 October) | Repo bound (27 October) | Repo bound (28 October) | Repo bound (29 October) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 29 | 210,495 | **216,568** | **6,073** | — | 226,520 | 231,638 | 231,697 |
| 30 | 222,073 | **224,280** | **2,207** | — | 238,808 | 247,960 | 247,972 |
| 31 | 241,079 | No improvement in this investigation | — | 241,742 | 246,638 | 251,110 | 251,188 |

The final certificates are `certificate29_216568.npz` and
`certificate30_224280.npz` in `certificates/`. Earlier successful certificates
are retained as intermediate results. The portable bundle contains only the
final certificates, their verification receipts, the independent verifier,
the coordinate materializer, the expanded exact coordinate files, and this
document. Nothing has been published
externally.

## The structural change

Use squared norm 32 and maximum off-diagonal inner product 16. Begin with the
196,560 minimal vectors of the Leech lattice. Choose an octad normal space of
dimension 8 and its orthogonal residual space of dimension 16. There are 4,320
relevant residual fibres. At each residual r the 16 original points are
r +/- p_i, where the eight p_i are orthogonal and have squared norm 16;
r also has squared norm 16.

Choose two normal axes p,q. Delete the four points r +/- p and r +/- q. For
orthogonal tail vectors a,b of squared norm 8, insert the following eight:

```
(r +/- (p+q)/2, +/- a)
(r +/- (p-q)/2, +/- b)
```

Signs are independent within each line. Every inserted point has squared norm
32. Its 24-dimensional head is the midpoint of two deleted Leech points.
Consequently its inner product with every retained Leech point is at most 16.
This is an exact convexity argument, not a numerical screen.

An elementary trade gains four points. The search groups residuals into 135
classes of 32, sharing a normal frame. One uniform piece applies the same trade
to all 32 residuals and gains 128 points. Multiple pieces in one group are
allowed when their normal and tail planes are orthogonal.

The extension that matters is allowing all orthogonal pairs of tail root lines,
rather than only the two diagonals of a coordinate plane. The final 29D result
uses D5 tail roots (20 lines, 70 orthogonal pairs); the final 30D result uses D6
tail roots (30 lines, 195 orthogonal pairs). D6 is not a subsystem of E6: moving
to the larger E6 root catalogue had discarded the successful 30D seed. Keeping
that seed and enlarging its D6 pair choices worked better in this experiment.

Separately append a pure tail kissing configuration of squared norm 32: 40 D5
roots in 29D, or 72 E6 roots in 30D. Its inner products with inserted points are
at most sqrt(32*8)=16 by Cauchy-Schwarz. The pure tail configuration need not
belong to the same root system as the inserted tails.

| Dimension | Uniform pieces | Deleted Leech points | Inserted midpoint points | Pure tail points | Total |
|---|---:|---:|---:|---:|---:|
| 29 | 156 | 19,968 | 39,936 | 40 | 216,568 | 231,697 |
| 30 | 216 | 27,648 | 55,296 | 72 | 224,280 | 247,972 |

Thus the totals are 196560 + 128*156 + 40 and 196560 + 128*216 + 72.

## Independent verification and reconstruction

The verifier uses Python and NumPy and imports no search code or search caches.
It regenerates the standard Leech shell from the binary Golay code, checks that
all deleted points are distinct shell members, checks distinctness and exact
norms of every inserted point, and checks the two deleted endpoint witnesses
for every midpoint. It directly checks **every pair of inserted points**.
The standard Leech-shell bound and midpoint argument certify all remaining
bulk interactions. The pure tail roots are independently generated and checked.

The default dot products use float64 BLAS on small integers. Every product and
partial sum is exact: the sum of absolute coordinate products is at most 32 by
weighted Cauchy-Schwarz. No tolerance is used. `--integer` selects integer
matrix multiplication as a slower alternative. Both final verification receipts
report `PASS EXACT` and maximum inserted inner product 16.

From the extracted bundle, run:

```powershell
python verify.py certificate29_216568.npz --output verification29_216568.json
python verify.py certificate30_224280.npz --output verification30_224280.json
python materialize.py certificate29_216568.npz points29_216568.npz
python materialize.py certificate30_224280.npz points30_224280.npz
```

Each materialized file contains integer arrays `a` and `b`; all coordinates
are exactly `a + b*sqrt(3)`. All points have squared norm 32. Divide by sqrt(32)
for a unit-sphere configuration. In 29D the array `b` is zero. In 30D only some
pure E6 roots use the radical coefficient. The compact certificate files contain
explicit inserted coordinates, explicit deleted Leech coordinates, local
endpoint indices, and metric weights; both final files have all weights one.

SHA-256 of the final certificate files:

```
29D 6380f42e112a7340be1a9164f28241fc8883e405fbd7afe4e4e9d5b09382e61e
30D 60261bc68ad9123636005a8a8c99493c31e63d9393ff1a59973b77f92fd6ec30
```

## What happened in 31D

The strongest E7 midpoint candidate in these short probes has 287 uniform
pieces, hence 233,422 points including 126 pure E7 roots. It passes exact
reduced compatibility checks but was not exported as a claimed record; it is
7,657 below 241,079. The trials establish neither an optimum nor an obstruction.

The following distinct representations were actually tested:

- Coordinate tail planes: 245 pieces, 228,046 points.
- All orthogonal E7 root pairs: 287 pieces, 233,422 points after short exchanges.
- Single-axis square trades: 510 uniform axis choices, 229,326 points.
- Directed head choices, dropping compulsory normal antipodes: 1,023 choices,
  229,422 points.
- A mixed catalogue of 141 tail lines, allowing tail inner products larger than
  4 when compensated by smaller head products: no gain over the 287-piece seed.
  The tail catalogue uses exact arithmetic in Q(sqrt(2),sqrt(3)); its failure in
  a short probe is not evidence that correlated tail catalogues are impossible.

The fixed-subspace nonlinear 544-head proposal from the earlier investigation
is not revived here: its original whole-coset model has a proved 528 ceiling.
That ceiling does not apply to all nonlinear head constructions.

## A concrete next change for 31D: split the 32-fibre choices

The current search forces 32 local choices to be identical. This is a large
restriction, and `fibre_census.py` measures its exact effect rather than merely
speculating about it. For two distinct residual groups, the largest residual
product is either 4 or 8. The complete census gives:

| Group relation | Other groups per group | Positive maxima per residual, among 32 |
|---|---:|---:|
| Maximum residual product 4 | 64 | 16 |
| Maximum residual product 8 | 70 | 4 |
| Same group, residual product 16 | 1 | 1 |

For E7 midpoint choices, a forbidden pair across groups of the second type
therefore produces a 4-regular bipartite conflict graph on their 32+32 elementary
trades. The uniform search rejects the entire pair of pieces. The elementary
search can instead make different choices at different residuals. This does
**not** imply an automatic gain: the bipartite conflict graph has a perfect
matching, so gains require coordinated trades involving several alternatives.

A useful next algorithm is a restricted master problem on elementary trades:

1. Expand the current 287 pieces to 9,184 elementary trades.
2. Generate new columns around coarse choices blocked by only a few incumbent
   pieces, using the exact residual Gram matrix to enumerate actual conflicts.
3. Optimize small overlapping neighbourhoods with the exterior fixed. Permit
   coordinated removals and insertions across several residual groups; optimize
   the net number of elementary trades, not the number of complete 32-blocks.
4. Reject neighbourhoods whose relaxation cannot give a positive gain. Export
   and verify explicit coordinates whenever a new count matters.

The first 31D record in this family requires 11,099 elementary trades, yielding
241,082 points. That is 1,915 more than the present candidate. The local
three-plane capacity is 12,960 elementary trades; reaching the target would use
about 85.6% of that capacity. This capacity is not an existence proof, but it
provides a precise target and a falsifiable direction that changes the model.

## Research files

`midpoint_planes.py`, `root_pairs.py`, and `d_roots.py` implement the three tail
models. `augment_roots.py` performs bounded one-to-two and two-to-three exchanges.
`export_coordinates.py` writes midpoint endpoint certificates. The other probe
scripts and JSON files preserve the unsuccessful 31D experiments and the exact
fibre census. Search and export use geometry caches from
`../creative2527_20261025/`; certificate verification and materialization do not.

Re-export the final selections in the existing workspace with:

```powershell
python research/creative2931_20261025/export_coordinates.py 29 --roots
python research/creative2931_20261025/export_coordinates.py 30 --d-roots
```

Searches were short structural probes and small exchanges, not optimality
computations. The new lower bounds are supported by the exported certificates,
independently of whether the search heuristic can be reproduced byte for byte.
