# November 2026 sextet exchanges in dimensions 28 and 29

**Verified lower bounds: 222,120 in 28D and 231,776 in 29D.**
These exceed PackingStar's November 11 release by **16 and 2** points.
They are modest improvements; the campaign has not obtained a substantial
increase across all dimensions 26–31.

| Dimension | Previous repository certificate | PackingStar, November 11 | This construction | Removed from PackingStar | Added |
| --- | ---: | ---: | ---: | ---: | ---: |
| 28 | 221,988 | 222,104 | **222,120** | 40 | 56 |
| 29 | 231,697 | 231,774 | **231,776** | 6 | 8 |

The retained points are credited to PackingStar. Its inputs are pinned to
commit [`2f3453b`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
Our contribution is the explicit exchange, discovered with the existing
corrected sextet isometries and bipartite union method. This is a continuation
of that method, not a new construction principle or a claim of global optimality.

## Exact geometry

The package supplies every coordinate. Every integer numerator has squared
norm **128**, and every distinct pair has product at most **64**. Division by
sqrt(128) therefore gives a kissing configuration. The denominator-two field
only records the scale relative to the source files.

The complete verifier checks norms, distinctness, actual opposite pairs,
paired-to-unpaired products and all remaining signed products. Antipodality
is not assumed. With P actual opposite pairs and U unpaired points, it covers

    4*binom(P,2) + P + 2*P*U + binom(U,2) = binom(2*P+U,2).

Every coefficient has absolute value at most 16 and the dimension is at most
31. Every coordinate product and partial dot-product sum is an integer with
absolute value at most 31*16^2 = 7,936 < 2^24. Thus the binary32 matrix
multiplication used here is exact integer arithmetic, without a tolerance.
The checker is inherited verbatim from the earlier joint-refill package.

`verify_delta.py` uses a separate direct-int64 coordinate sum for every new
point against every final point. It also reconstructs the complete pinned
source from the final configuration and the explicit removals/insertions,
then checks its normalized canonical SHA-256. Neither proof assumes the
PackingStar input is geometrically valid; the full final check establishes
all retained and changed inequalities.

## Replay

Python 3.9+ and NumPy suffice. No optimizer, network, GPU or search state is
needed. From this directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python final.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python negative_controls.py
```

An optional dimension argument to `verify.py` selects 28 or 29. Reserve at
least one available CPU as described in the repository's `AGENTS.md`.

The four controls must reject a wrong norm, a duplicate, an illegal pair of
unpaired points, and an illegal interaction involving an actual opposite pair.
Recorded wall times are evidence of the executed replay, not mathematical inputs.

## Discovery and remaining scope

A Golay sextet partitions the first 24 coordinates into six tetrads. On
each tetrad apply J/2−I, then negate the first tetrad; the tail coordinates
are fixed. The corrected transformation is orthogonal and preserves the
Leech shell. Both the current configuration and its image are valid, so
the exclusive parts have a bipartite conflict graph. A maximum matching
and alternating paths select a largest compatible union of the admitted
vertices.

In 28D, the second frame deletes 40 points and inserts 48; the third frame
adds eight further points. All 35 frames in the selected normal-octad
family were tested with a 32-blocker candidate cutoff. In 29D, the eighth
frame gives the six-for-eight trade. The cutoff and frame family restrict
the search only; they do not restrict the complete coordinate proof.

The dated source review and wider experiments are recorded in
[`../../../research/reclaim_26_31_20261112/`](../../../research/reclaim_26_31_20261112/).
They distinguish exact findings, numerical timeouts and conjectural next
steps. No worldwide-priority or unrestricted maximality result is asserted.

## File inventory

| File | Purpose |
| --- | --- |
| `final.py` | Derive the registered counts from the complete arrays. |
| `verify.py` | Check the manifest and replay both proof routes. |
| `verify_geometry.py` | Complete exact finite geometry. |
| `verify_delta.py` | Independent direct-integer changed pairs and source reconstruction. |
| `negative_controls.py` | Four semantic rejection checks. |
| `data/coordinates*.npz` | Complete norm-128 integer configurations. |
| `data/exchange*.npz` | Explicit additions and removals. |
| `data/provenance.json` | Source commit, hashes, counts and integer replay results. |
| `data/verification*.json` | Executed complete-geometry reports. |
| `data/manifest.json` | SHA-256 for every other file in this package. |
