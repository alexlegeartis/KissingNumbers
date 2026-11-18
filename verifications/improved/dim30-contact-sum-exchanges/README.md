# 248,351 points in dimension 30

The integer coordinate witness proves `τ(30) ≥ 248,351`. It adds 62 points
to the [248,289-point repository witness at revision 315f783](https://github.com/alexlegeartis/KissingNumbers/tree/315f7838e8e453f14c8715f22917c77f46fa3db3/verifications/improved/dim30-31-midpoint-exchanges/).

For source points `x` and `y` with squared norm 128 and inner product −64,
the contact sum `x + y` also has squared norm 128. We priced each admitted
candidate against the current witness and searched for mutually compatible
sets whose shared deletions cost fewer points than their insertions. The first
round also reused earlier midpoint and auxiliary candidates. The second round generated contact sums from points introduced in the first round
and repriced the earlier catalogue against the new witness.

| Configuration | Points |
| --- | ---: |
| Pinned source | 248,289 |
| First round | 248,332 |
| Second round | 248,351 |

The frozen exchange inserts **188 points** and removes **126** from the pinned
source. The verifier reconstructs the final set, checks **46,672,222** pairs
involving an inserted point in signed int64, and checks all **30,838,985,425**
unordered pairs of final points. Every point has squared norm 128 and the
largest distinct-point inner product is 64. Dividing all coordinates by
`√128` gives the unit-sphere configuration.

From the repository checkout, run the verifier with NumPy.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verifications/improved/dim30-contact-sum-exchanges/verify.py
```

`final.py` reads the frozen coordinate count for the repository registry.

[Coordinate and exchange checksums](manifest.json) · [Exact verification report](verification.json)
