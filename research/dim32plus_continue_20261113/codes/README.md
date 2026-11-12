# Compatible sign-family continuation

This lane extends the signed weight-sixteen construction in dimensions 61–67.
The complete proof is in
[`verifications/improved/dim61-67-compatible-sign-families`](../../../verifications/improved/dim61-67-compatible-sign-families/).
The search scores and intermediate families below are historical discovery
records; the registered counts come from that package's exact verifier.

## What changed

The published global family had 1,835 added supports; shortening that globally
chosen family left 588 added supports in 61D/62D. Optimizing after shortening
removes constraints from old supports which no longer exist. This produced 868
added supports with the old search rules, then 937 with the new compatibility
rule. Both 61D and 62D gain 714,752 points over the preceding certificate.

The new rule permits two added supports to meet in ten positions when their
Hamming duals share an eight-position check there, with opposite affine
parities. The eight old/new checks on a candidate span only a rank-four dual.
The last check completing a rank-five Hamming dual is a search variable. This
freedom makes new/new repairs possible while retaining all 2,048 signs.

The global search first reached 1,873 added supports. Filtering against the 66D
tail and optimizing again reached 1,938. A final global refill reached 1,940,
so the full sparse layer has exactly 32,768 support classes and 2^26 points.
The 67D full-block search retained 1,792; a separate proper-affine-subcode lane adds 4,624 points and
is documented in [`../middim/`](../middim/).

## Frozen full-class searches

| Prefix | Candidates | Added supports | Budget or status |
|---|---:|---:|---|
| `short61` |34,272|868|600-second direct-shortening search|
| `global64` |112,896|1,840|310-second initial fill/correct continuation|
| `global_repaired` |112,896|1,873|600 seconds with new/new repairs|
| `short_repaired` |34,272|937|600 seconds with new/new repairs|
| `expanded66` |62,720|1,938|480 seconds|
| `expanded67` |18,816|1,792|480 seconds; no full-block gain|
| `global_final` |112,896|1,940|240 seconds from the 66D family|

These are finite searches, with no claim that any family or dimension is
optimal. The nonzero-tail filters initially forbid every intersection above
eight; they are discovery restrictions rather than a general obstruction.

## Candidate construction and reproducibility

`build_affine_pool.py` reconstructs the affine-plane incidence structure from
the 72 supplied dual octads, derives its field coordinates, and enumerates the
images of a supplied support under affine transformations. The 225,792 images
give 112,896 distinct candidates. Every candidate is explicitly checked for
membership in the actual 64-bit code dual; no ambient code automorphism is
assumed merely from the plane description. `affine_plane.json` records the
field table, point labels and finite orbit data.

`build_records.cpp` computes every original/candidate intersection. It checked
3,480,357,888 pairs and recorded the eight ten-overlaps per candidate. These
large scratch pools are reproducible and ignored by git. Compile with C++17
and OpenMP, then run `prepare_seed.py` before the search. The final exact
certificates do not require rebuilding this pool.

The search engine accepts

```text
seconds candidates.bin prefix old_supports.bin removed_mask seed odd repair
```

For the shortened family the removed mask is exactly `45040394336993280` and
`odd=1`; the full 64D mask and parity flag are zero. `repair=1` enables new/new
ten-overlap repairs. Search logs preserve budgets, seeds and accepted counts.
Each worker was pinned to a single CPU; at least one CPU remained reserved.

## Files

- `build_affine_pool.py`, `affine_plane.json` and `build_records.cpp` reconstruct
  the finite candidate pool and its exact original-support conflicts.
- `prepare_seed.py` converts the published explicit sign codes into search
  seeds and solves the shortened odd-phase equations.
- `shared_search.cpp` is the initial direct-shortening fill/correct search.
- `shared_search_adaptive.cpp` adds adaptive dual completion and new/new repairs.
- `filter_expanded.cpp` and `prepare_expanded.py` construct the 66D/67D strict
  extra-support filters and retain the stronger compatible starting family.
- `export_family.py` discards search scores and indices and writes actual
  physical supports, parity checks and affine shifts.
- `verify_family.py` and `check_family.cpp` are the finite proof prototypes;
  the latter freshly enumerates all 2^28 parent words and checks every pair.
- `geometry_counts.py` proves the 61D/62D geometry and layer counts.
- `section63.py` checks all 72 octad sections and directly enumerates sign words
  on the chosen section.
- `expanded_geometry.py` and `check_expanded.cpp` verify the tail extension,
  including direct scalar checks of all retained parent-support pairs.
- The `*_certificate.json`, `*_verification.json`, `*_geometry.json` and
  `*_section63.json` files preserve intermediate explicit witnesses and replays.
  `seed_preparation.json` and `expanded_preparation.json` record the pools.
  The `*_search.log` files record accepted improvements and the finite stopping
  points. The final proof package, including its manifest, is authoritative.

The counts 61:88,975,146;62:122,539,020;63:184,393,988;64:335,711,104;
65:335,779,330;66:342,271,364 are strictly above the starting repository
certificates. Literature priority is limited to the dated comparison recorded
in the campaign's `source_review.json`.
