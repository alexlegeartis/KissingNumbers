# 26D: 201,201 kissing points from a 27D section

This gives **tau(26) >= 201,201**, 9 above the 201,192-point configuration in
[R. Takhanov’s released certificates](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1).
It is an explicit lower bound, not a proof of optimality.

## Construction

All points have squared norm 128 and distinct points have inner product at
most 64. Slice Takhanov’s 27-dimensional construction at coordinate 25 = 0
(zero-based, the middle of the three tail coordinates) and drop that
coordinate. The resulting 9,028 off-equator rows allow 192,144 Leech roots
after all unblocked roots are restored.

Let I = {0, 2, 4, 5, 6, 10, 11, 23} and let J be the other 16 Leech coordinates.
Keep the 200,864 rows with nonzero J projection. Replace the remaining 308
rows with a new 337-point code in the ten coordinates I plus the two tails.
The frozen rows project to only 2,601 distinct linear constraints on the
new code. The replacement has 250 integer rows and 87 nonintegral rational
rows. The rational directions were found by maximizing linear forms over the
constraint polytope, rationalizing outward directions on the sphere, and
selecting compatible subsets by exact conflict-graph exchanges. Verification
does not require rerunning the search. The latest kernel adds 6 points over the first published 331-point kernel, by reoptimizing neighborhoods of exact conflicts among the rational candidates.

## Replay

Run from this directory with Python 3 and NumPy:

```sh
python verify.py --full
```

The extra `--full` checks every pair in the 200,864-row integer exterior;
without it the verifier uses the standard Leech minimum-distance bound for
the shell-shell pairs, but still checks every other constraint. All
integer-row dot products computed using binary32 GEMM are exact, since each
entry and partial sum is an integer of absolute value at most 128.
Rational pair checks use integer arithmetic, including exact numerator and
denominator comparisons.

`section.i8.gz` stores the 9,028 selected parent rows as gzip-compressed
row-major signed bytes, each row with 26 coordinates.
`kernel.txt` stores ten integer numerators followed by a positive denominator
per row. Its other 16 coordinates are zero.

For an explicit parent provenance check, optionally pass the original
[upstream certificate](https://github.com/k-nic/KissingNumbers/blob/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1/certificates/data/d27_certificate.npz):
`python verify.py d27_certificate.npz --full`. This additionally checks
the upstream SHA-256 and equality of the extracted rows with the unmodified
parent.
