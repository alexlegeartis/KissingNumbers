# Independent November sextet proof review

This read-only review covers the 222,120-point 28D and 231,776-point 29D
certificates in
[`dim28-29-november-sextets`](../../../verifications/improved/dim28-29-november-sextets/).
It found no blocking error in the complete geometry verifier or the independent
direct-integer delta check. The certificate package remains the source of the
complete large-configuration verification; this review does not claim a second
large all-pairs brute-force run.

## Mathematical and implementation review

The geometry verifier first checks integer coordinates, squared norm 128,
uniqueness and coefficient bounds. Each binary32 multiplication and every
possible partial sum is an integer of absolute value at most
`31 * 16^2 = 7936 < 2^24`. Consequently these operations are exact on standard
CPU binary32 arithmetic, independent of summation order. The same argument
applies to binary64. There is no angle tolerance.

The antipodal reduction uses a line twice only after row uniqueness proves
that both signs actually occur. Absolute products cover four point pairs
between two paired lines, and two pairs between a paired line and an unpaired
point. Products between two unpaired points retain their signs. The triangular
block loops cover each unordered pair, and the independent coverage identity
equals `N(N-1)/2`. The negative product fixture specifically checks that an
unpaired product of -96 is accepted.

The delta replay reconstructs the pinned source set from the final points,
additions and removals. It checks every added point against every other final
point by direct int64 coordinate sums, without BLAS. This supplies an
independent check of every changed inequality. Unchanged pairs are covered by
the complete geometry verifier; the delta check by itself is not a complete
new proof of the retained source.

## Reproducible checks

`review.py` independently loads the original PackingStar source files,
checks their full file SHA256 hashes, normalizes them from squared norm 32
to 128, and uses an independent lexicographic array sort to reproduce the
source-set hashes. It then invokes and checks the direct-integer delta replay.
The verified counts are 12,438,664 changed comparisons in 28D and 1,854,200
in 29D, with maximum product 64 in each case.

Nine synthetic fixtures cover valid paired, unpaired and mixed sets; a valid
large negative product; a bad norm; a duplicate; and violations among paired
lines, paired/unpaired points and two unpaired points. A separate direct
Python-integer all-pairs predicate supplies the expected result. Testing both
binary32 and binary64 at block sizes 1, 2, 3, 7 and 64 gives 90 comparisons.
All agree. `results.json` also pins the reviewed geometry and delta modules.

The review originally identified one wrapper issue: an unsupported dimension
could select no records and return success. The package author added an
explicit unsupported-dimension guard. This review script does not mutate or
refresh the certificate package's manifest. Run that package's top-level
verification command for its manifest and command-line checks.

From this directory, with NumPy installed and the pinned PackingStar files
available, run:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  taskset -c 0-1 python review.py \
  --source-dir /path/to/PackingStar/25D-31D/new_bounds_20261111 \
  --output /tmp/november-sextet-review.json
```

The script uses repository-relative paths for the proof package, or accepts
`--package /path/to/dim28-29-november-sextets`. Source files are required from
PackingStar commit `2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`; its two URLs and
expected hashes are in the certificate package's `data/provenance.json`.
No source coordinates are duplicated here.
