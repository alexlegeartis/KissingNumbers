# 223,528 kissing points in dimension 28

This package gives an exact finite certificate for

\[
\tau_{28}\ge 223528.
\]

It improves Rustem Takhanov's **223,496**-point configuration by **32** and
the preceding repository certificate, **222,120**, by **1,408**. The comparison
is with [Takhanov's version1 release](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1),
published 17 November 2026. This is an explicit lower bound and a comparison
with a pinned source, without an exhaustive worldwide-priority claim.

## Construction

Keep the source outside the twelve-dimensional subspace consisting of the
octad `I = [0,2,4,5,6,10,11,23]` and the four added tail coordinates. Replace
the 680 points in that subspace with 712 points. Against the full source, the
actual exchange is **320 additions and 288 removals**; 223,208 points remain.

The search simultaneously changed retained E8 roots, lifted midpoints, and
auxiliary head frames. Fixing the retained E8 roots had prevented this gain.
The pure E8 midpoint layer doubles from 64 to 128 points at a cost of 32
additional root deletions; the auxiliary and pole counts stay the same.
The new midpoint endpoints form sixteen disjoint four-cycles on 64 deleted
E8 roots, with two tail lifts per edge. The source used eight such cycles.
Exactly 192 of the old auxiliary points are replaced by 192 different points.

| Layer | Takhanov | This certificate |
| --- | ---: | ---: |
| Retained Leech shell | 171,184 | 171,152 |
| Fibre midpoints | 50,688 | 50,688 |
| Row midpoints | 384 | 384 |
| E8 midpoints | 64 | 128 |
| E8 auxiliary points | 384 | 384 |
| Barnes-Wall auxiliary points | 768 | 768 |
| Pure tail poles | 24 | 24 |
| **Total** | **223,496** | **223,528** |

All stored coordinates are integers. They satisfy `x.x = 128` and
`x.y <= 64` for distinct points. Divide by `sqrt(128)` for unit vectors.
The inherited `coordinate_denominator=2` metadata is a lattice-scale label,
not the unit-vector normalization.

## Complete independent verification

Only Python and NumPy are needed; no solver, search output, external checkout,
network access, or unverified Leech-lattice premise is required:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python negative_controls.py
python final.py
```

Reserve at least one logical CPU when running concurrent computations.

`verify.py` checks file identities, source-exchange reconstruction, cardinality,
distinctness, every norm, and **all 24,982,271,628 unordered final pairs**.
It first performs a separate signed-int64 check of all **71,477,600** changed
unordered pairs. No unchanged-source geometry is assumed by the complete
final check. All maximum inner products are exactly 64.

The complete checker first identifies actual antipodal pairs after rejecting
duplicates. Here there are 111,640 pairs and 248 unpaired vectors. Absolute
products are used only for comparisons involving a complete antipodal pair;
unpaired-to-unpaired comparisons retain their actual signs. The pair-coverage
identity is checked explicitly.

The fast kernel is exact binary32 arithmetic on small integers. Every product
and every possible partial sum has absolute value at most
`28*8*8 = 1792 < 2**24`; acceptance uses the exact inequality, with no tolerance.
`python verify.py --direct` replaces that kernel with signed-int64 arithmetic
for a slower replay. Pass `--report PATH` to save a fresh receipt. An additional full binary64 run using the source's
independent geometry checker agrees with the standalone proof.

`negative_controls.py` verifies a valid small control and rejects duplicates,
incorrect norms, a norm-correct geometric conflict, noninteger input, a false
exchange reconstruction, and removal of a point absent from the source.

## Files and provenance

- `data/coordinates28_223528.npz`: the complete final coordinate witness.
- `data/takhanov_d28_223496.npz`: Takhanov's unmodified source coordinate file,
  retained for attribution and exact reconstruction of the exchange.
- `data/exchange.npz`: the 320 added and 288 removed integer rows.
- `manifest.json`: source revision, counts, normalization, file SHA-256 values
  and canonical coordinate SHA-256.
- `verify.py`: standalone full geometry and independent changed-pair proof.
- `negative_controls.py`: six semantic rejection controls.
- `final.py`: read the frozen certificate for the repository registry.
- `verification.json` and `verification.log`: the completed standalone replay.

The canonical final-coordinate SHA-256 (numerically lexicographic rows,
little-endian int16, row-major) is

`19aed25d116134936541f5390be1085f12893478d174bbc3d2ac197318b29077`.

Source credit: Rustem Takhanov, *Lower bounds for kissing numbers between
dimensions 25 and 31* (2026), repository `k-nic/KissingNumbers`, commit
`dadcb2940f7c3c3bd478b0543931d3eef9f4bd65`,
`version1/certificates/data/d28_coordinates.npz`. The new configuration is a
continuation of that construction. The source and the development history
are distinguished throughout.
