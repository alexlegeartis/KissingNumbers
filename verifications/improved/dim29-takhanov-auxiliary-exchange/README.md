# 232,154 points in dimension 29

The integer coordinate witness proves `τ(29) ≥ 232,154`, one above
[Takhanov’s alternate 232,153-point witness at revision dadcb29](https://github.com/k-nic/KissingNumbers/blob/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1/certificates/data/d29_alt_coordinates.npz).

The exchange is supported on the eight octad coordinates
`[0, 2, 4, 5, 6, 10, 11, 23]` and the five tail coordinates. It removes two
E8 roots, inserts four lifted E8 midpoints, and replaces 19 auxiliary points
with 18. The four midpoint additions force six auxiliary deletions. A search
over the auxiliary points recovered a compatible placement of 18 points
within the 19-point deletion scope.

The frozen 22-for-21 exchange is checked against the byte-identical source
archive bundled here. Signed int64 replay checks **5,107,135** unordered pairs
involving an inserted point. Complete exact geometry covers all
**26,947,623,781** unordered pairs in the final configuration. The squared
norm is 128 and the largest distinct-point inner product is 64, so scaling
by `1/√128` gives a kissing configuration.

From the repository checkout, run the verifier with NumPy.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verifications/improved/dim29-takhanov-auxiliary-exchange/verify.py
```

`final.py` reads the frozen coordinate count for the repository registry.

[Coordinate and exchange checksums](manifest.json) · [Exact verification report](verification.json)
