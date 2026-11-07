# 370,408 points in dimension 33

**7 November 2026.** This self-contained package proves a kissing lower bound
of **370,408** in dimension 33, **32 above** the actual published repository
certificate of 370,376 and 168 above the pinned 370,240-point Lin construction.
The 370,392-point witness is an **unpublished intermediate**, not the previous
repository record. All historical certificates and the paper remain unchanged.
There is no worldwide-priority, global-optimality, source-union optimality or
finite-pool completeness claim.

## Exact point convention

`certificate33_370408_explicit.json` has SHA-256
`5686dfcd2b891bf6d9be6b01013d8205774aa1a19c02a2c45bdfd86194314a60`
and format `signed_octad_star_points_explicit_v1`. It specifies:

- 131,072 dense sign vectors from a rank-17 binary code of minimum distance 8.
  **Only coordinates 0 through 31** carry these signs; coordinate 32 is zero.
- 1,837 full-even octad supports, each with all 128 even-negative sign patterns
  of magnitude 2, giving 235,136 points.
- 2,064 individual exceptional octads. Each `[support, negative_mask]` pair
  specifies one row of eight entries of magnitude 2. No antipode is implied.
- The complete 2,112 signed roots, with two entries of magnitude 4, minus the
  **156 individual rows** in `removed_root_signed_points`, leaving 1,956 roots.
- **180 individual rows** in `star_signed_points`, each with one entry of
  magnitude 4 and four entries of magnitude 2. These are not antipodal lines
  or complete sign families.

Thus 131,072 + 1,837 × 128 + 2,064 + 1,956 + 180 = **370,408**.
Every row is an integer vector in R^33 with squared norm 32; distinct rows
have **ordinary inner product at most 16**. Dividing by sqrt(32) yields unit
vectors separated by at least 60 degrees. The exact rank is 33, proved by
rank 33 modulo 101. Global uniqueness is independently checked.

Exactly 500 points lack antipodes: 372 exceptional octads, 12 retained roots
and 116 stars. Valid pairs with product -24 occur. Testing absolute products
would wrongly reject this configuration. Adding missing antipodes can create
forbidden +24 products; completing the deletion list would remove 168 roots
rather than 156. Neither operation is part of this construction.

## Reproduction

Python 3.10+ and NumPy suffice. From any working directory, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python /path/to/nonantipodal33/replay.py
```

Use an available CPU affinity that leaves at least one logical CPU free, as
required by the repository's `AGENTS.md`. The replay verifies a complete
SHA-256 manifest and runs in a temporary copy, leaving the source unchanged.
It uses no network, solver, search executable, cached proof or generated
coordinate archive. Do not use Python's `-O` option.

Individual proof entry points are:

```sh
python verify.py
python verify_coordinates.py
python independent_core.py --output /tmp/verify33-intermediate
python independent_audit.py --output /tmp/verify33-final
```

`verify.py` is the producer's standard-library integer proof.
`verify_coordinates.py` independently materializes coordinates and uses int64
products. The separate mathematical auditor, `independent_core.py` with
`independent_audit.py`, imports neither producer verifier nor any search,
matching, graph or optimization function. Its two original source files and
original reports are frozen in `evidence/` for provenance, not replay entry
points. Portable adaptations only relocate inputs, change report destinations,
omit saving redundant row arrays, and transparently decompress the provenance
pool; the finite-geometry proof is unchanged.
The arrays are fully regenerated in memory and their hashes must match the
original independent audit.

The intermediate and final audits each run 16 deliberately corrupted-input or
semantic rejection controls. The replay also mutates the final certificate's
bytes and requires the independent auditor to reject its hash. Controls cover
individual rows, mixed line/point schemas, malformed magnitudes, masks,
counts, dense rank and placement, invalid restored roots and antipodal
completion. The intermediate additionally checks dense minimum distance and
full-layer/star parity. These checks guard interpretation, not only file hashes.

## Complete pair coverage

All five categories have exact norm and shape checks. The auditor covers:

1. Dense/dense through complete codeword enumeration, rank 17 and minimum
   distance 8, yielding maximum product 16.
2. Full/full through all support intersections (at most 4) and all 128 even
   sign patterns. Same-support sign distance is at least 2.
3. Exceptional/exceptional through **2,129,016 unordered pairs**, each tested
   by both Python-integer masks and int64 coordinate products. Their extrema
   agree. Exceptional/full checks every support overlap.
4. Root/root through int64 coordinates. Root/octad products are at most 16,
   dense/root at most 8, and dense/octad at most 16 by verified magnitudes.
5. Every star/dense, star/full, star/exceptional, star/root and distinct
   star/star product by int64 coordinates. The respective maxima are
   12, 16, 16, 16 and 16.

The final independent audit checks 70,698,156 unordered or cross-category
coordinate products in the nonautomatic categories, in addition to the
scalar exceptional checks and full-support proofs. The scalar producer uses
an exact parity-maximization formula for star/full interactions as another
proof route. No floating-point tolerance is used.

## Provenance and exact gains

`sources/` contains byte-identical copies of the published 352,992-point 32D
source, the classic 370,368-point 33D baseline and the actual 370,376-point
root-star certificate. Their original paths and hashes remain in the frozen
delta records, and the portable auditor maps them to local inputs.

The first trade embeds the 32D sparse coordinates into coordinates 1 through
32 while preserving the old 33D dense support at 0 through 31. The independent
source reconstruction checks all 3,280 source IDs. Relative to the actual
370,376 certificate, `delta_370392.json` adds 1,648 octads and removes 1,632,
gaining 16. All old dense points, roots and stars are retained at this stage.

The final `delta_370408.json` keeps that intermediate's dense and octad rows,
adds 148 stars, removes 32 stars and removes 100 additional roots, gaining
another 16. Relative to the published 370,376 control, the final set has
1,796 added and 1,764 removed rows, a net gain of **32**. All set differences
are recomputed from coordinates. The independently regenerated complete root
universe and 336 selected pool IDs also reproduce the final certificate.
`pool.json.gz` is deterministic gzip of the original `pool.json`, whose exact
decompressed SHA-256 is
`7827321d836a37965b1d3d0dd884cf7d5b8fb26532d37b1c63a4dc4dbbbf3fa9`.
The reader checks that hash after decompression. The original audit evidence
retains the raw filename and hash. This is only compact discovery provenance;
its completeness and any optimizer claims are unnecessary and untrusted for the lower bound.

`evidence/public_comparison.json` records the 7 November 13:14 UTC read of
main `0d06328e82838e2eeb32590943a19aa572fb337e`, whose RESULTS table still had
370,376. The publication preparation fresh-read the same head. This is a
repository comparison, not a new comprehensive literature or priority audit.
The Lin comparison and dense-code attribution are retained in the parent
package's `source/public/` and its README.

The underlying public support codes are due to Qihong Lin, the binary
[32,17,8] generator to Cheng–Sloane, and the layered construction to
Edel–Rains–Sloane. The new claim concerns the explicit finite configuration
and the exact improvement over the stated controls.
