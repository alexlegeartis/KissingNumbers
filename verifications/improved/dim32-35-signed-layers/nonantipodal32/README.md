# Non-antipodal signed octads: 352,992 points in dimension 32

**PASS, 7 November 2026.** Both frozen 352,992-point certificates and the published
352,832-point baseline reconstruct as distinct norm-squared-32 integer vectors
whose ordinary pairwise inner products are at most 16. Therefore both new
certificates establish the constructive lower bound τ(32) ≥ 352,992, a gain of
160 over that exact published baseline. This is not an optimality or worldwide
priority claim.

## Reproduce

This directory is self-contained. With Python 3.10+ and NumPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python replay.py
```

Pin to one available CPU if sharing a machine. The replay first checks the complete
`SHA256SUMS.txt` manifest, then uses an isolated temporary copy. It runs the
coordinate checker, the independent audit, the fiber classification check, and
an additional certificate-hash corruption control. No solver, search output,
network access or parent package is required. Do not run with `python -O`.

`verify_coordinates.py` is the producer's independently expressed coordinate
checker. `independent_audit.py` imports no producer, search, or earlier checker
code. It independently enumerates both certificates and the prior baseline,
checks all pair categories, compares scalar-mask and int64 exceptional Gram
histograms, reconstructs all coordinate rows and runs 14 semantic controls.
The original independent audit and receipt are frozen in `evidence/`; the
portable audit changes only its input paths. `report.json` is the portable
replay receipt. Replays leave all checked-in inputs and receipts unchanged.

## Inputs and populations

The immutable input paths and SHA-256 hashes are pinned in `independent_audit.py` and
`report.json`. Every input is hashed before work and again after controls.

- Fixed-body certificate: SHA-256
  `ab06963895d074d8dd6162bea2298796391dc876330ae997a5324adfb6638c5e`.
  Its 131,072 dense points, 218,112 full-layer points on 1,704 supports,
  1,824 explicit points on 29 supports, and 1,984 roots total **352,992**.
- Original discovery certificate: SHA-256
  `a21b8ee4086cdbbc90651bcab378b888d0c3f10643aa0d63b5d208cafe471d1a`.
  Its 131,072 dense points, 217,728 full-layer points on 1,701 supports,
  2,208 explicit points on 32 supports, and 1,984 roots also total **352,992**.
- Published baseline: SHA-256
  `69c3f368a4aae5956ca75a8569f761224148831e8a85edd5972a0b3b057af646`.
  The baseline has 1,664 exceptional points and totals **352,832**.

The two new certificates are not identical coordinate sets. Normalization
recovers three full supports, and changes the 128 odd-parity patterns on support
1099442200 to its 128 even-parity patterns. The complete octad point sets differ
by exactly those 128 removals and 128 additions.

The fixed-body certificate leaves the dense generator and all 1,704 full supports
unchanged from baseline. Its exact exceptional delta is **1,316 removed,
1,476 added, 348 retained, net +160**. It retains 15 old exceptional supports,
introduces 14 and drops 10. Every removed, added and unpaired point matches the
frozen `fixed_body_delta.json` exactly.

## Complete geometric proof scope

All 131,072 binary dense words are enumerated independently in Gray order. Binary
elimination gives generator rank 17; enumeration gives distinctness and minimum
nonzero weight 8. Thus the dense code has parameters [32,17,8]. Every support and
negative mask is checked for type, range, weight and containment.

All full-support pairs and all full/exceptional support pairs have intersection
at most four. The 128 even sign patterns are reconstructed and their full Gram
matrix checked. Every exceptional unordered pair is checked twice: once by
Python integer mask arithmetic and once by an independently reconstructed int64
coordinate Gram matrix. This means 1,662,576 pairs for the fixed-body certificate,
2,436,528 for the discovery certificate, and 1,383,616 for baseline. The two exact
histograms agree in each case. All 1,967,136 unordered root pairs are checked in
int64 coordinates.

All ten unordered pair categories are covered:

1. Dense/dense: 32 − 2d ≤ 16 from the enumerated linear-code differences.
2. Dense/full and dense/exceptional: eight contributions of magnitude two.
3. Dense/root: two contributions of magnitude four, giving at most eight.
4. Full/full: exact same-support sign Gram; different-support intersection ≤ 4.
5. Full/exceptional: every support intersection ≤ 4, each contribution magnitude four.
6. Full/root and exceptional/root: at most two contributions of magnitude eight.
7. Exceptional/exceptional and root/root: complete exact pair checks above.

The full coordinate sets are materialized and all global norms, distinctness and
antipode counts checked. Each new configuration has exactly **380 points lacking
antipodes**. The fixed-body exceptional set has 1,444 paired points, or 722 pairs.
The maximum distinct-pair inner product is exactly 16 and the minimum squared
Euclidean distance is 32. Scaling vectors by 1/√8 gives sphere centers at distance
two from the origin and mutual distance at least two.

This is a complete proof by pair categories, not a full 62-billion-pair Gram
sweep. All decisive arithmetic is exact Python integer or NumPy int64; products
have magnitude at most 32. No unenumerated lattice minimum or population premise
is used.

## Why non-antipodal semantics matter

Each entry in `exceptional_signed_points` means **one point**, not an antipodal
line. There are **768 valid unordered exceptional pairs with dot product −24**
in each new certificate. These are not antipodal pairs. An absolute-product
threshold would incorrectly reject them.

For example, the existing points (support, negative mask)
`[142637192, 8200]` and `[142637192, 142626944]` have dot product −24.
The first point's missing antipode is `[142637192, 142628992]`; adding it produces
a forbidden +24 product with the second point. The checker verifies this
explicitly and rejects the altered certificate.

The 14 controls comprise 13 rejected corruptions and one explicit false-rejection
demonstration for the absolute-product rule. They cover duplicate explicit points,
wrong count, deficient dense generator, weight-one generator, invalid masks,
duplicate or overlapping full layers, mixed-layer conflicts, an actual
nonduplicate sign corruption with product +24, and adding a missing antipode.

## Fresh read-only comparison

Saved snapshots use GitHub connector reads only; no remote writes were made.

- `evidence/main_head_snapshot.json` resolves `alexlegeartis/KissingNumbers` main directly
  to `2017d94a57a564382ec8549f3f4592944f6bd5a6`.
- `evidence/public_source_snapshot.json` contains fresh main/default-branch README,
  RESULTS and baseline certificate reads. Both tables still report 352,832 in
  dimension 32. The remotely read certificate hashes to the exact baseline
  SHA-256 above. Read time: 2026-11-07T12:44:03Z.
- `evidence/lin_source_snapshot.json` records a fresh Lin README at 352,320 and latest
  commit search at `f367a9c5a051c0d6821ce71324b9dfa1ead76c2e`. The new count is
  672 above that named comparison.
- Cohn's leaderboard was not re-fetched in this audit. These dated source checks
  do not establish exhaustive worldwide novelty or priority.

Public links:

- https://github.com/alexlegeartis/KissingNumbers/commit/2017d94a57a564382ec8549f3f4592944f6bd5a6
- https://github.com/alexlegeartis/KissingNumbers/blob/main/RESULTS.md
- https://github.com/QH-Lin/kissing-number-32/tree/f367a9c5a051c0d6821ce71324b9dfa1ead76c2e

## Explicit-sign schema

Masks are unsigned 32-bit integers, with bit j referring to coordinate j. Each
entry `[support, negatives]` in `exceptional_signed_points` encodes one vector:
coordinate j is zero outside support, -2 on negatives, and +2 otherwise.
Every support has weight eight and negatives is a subset of support. This
format never supplies an implicit antipode. The legacy baseline's
`exceptional_antipodal_lines` does supply both signs and is decoded separately.

The dense layer contains the 131,072 sign vectors from all binary combinations
of the 17 `dense_generator` masks, with bit one encoded as -1. Each full-even
support supplies all 128 even-parity sign assignments with entries ±2. Roots
are all vectors ±4e_i ±4e_j for i<j. There are no other layers.

The immutable fixed-body certificate includes the metadata field
`equivalent_frozen_candidate_sha256`. Here “equivalent” means the same certified
point count only. The two configurations are not identical coordinate sets;
the exact difference is described and independently checked above.

## Fiber structure and limits

`verify_fibers.py` exactly checks the frozen classification: 14 affine sign
fibers and 15 nonlinear fibers. All 15 nonlinear fibers have partial
antipodality. For each support, it enumerates the 256 translations in the binary
sign cube and computes the translation kernel. A nonempty fiber is affine
exactly when its kernel has the same size as the fiber. These observations are
properties of this witness, not proofs that every construction with this many
points must be nonlinear or non-antipodal.

The 1,664-point capacity proof in the preceding `exceptional32` package concerns
its fixed 25-support set only. It does not constrain the new 29-support set,
which removes 10 old supports and introduces 14. No optimality of the new
support pool, sign selection, or global kissing number is claimed.
