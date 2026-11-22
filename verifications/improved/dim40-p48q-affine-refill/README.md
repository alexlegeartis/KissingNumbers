# 40D affine refill in the P48q complement

This construction gives **1,126,351** kissing points, an increase of 24
over the preceding affine construction. It uses the rank-40 complement
of \(\sqrt{3}E_8\) in \(P_{48q}\).
[Preceding certificate](../dim40-affine-centre-refill/).

| Coefficient squared norm | Core or complete short families | Refill |
|---:|---:|---:|
| 54 | 1,084,590 | 0 |
| 36 | 40,590 | 30 |
| 30 | 1,134 | 0 |
| 6 | 0 | 7 |

The search enumerated the complete shortest dual shell. In the
\(P_{48q}\) complement, 16 discriminant classes whose shortest normal
representatives have squared norm \(8/3\) contain no bodies of squared
norm \(10/3\). If a candidate of coefficient squared norm 36 has one of
those classes as its difference from the core class, a forbidden
core contact would force a difference of coefficient squared norm 30.
The empty shell makes the candidate compatible with the entire core.
The search then selected the refill vectors jointly against the
short families and the available norm-six points.

Dorofeev, Sun and Wang supplied the rank-40 lattice. The basis used here
is an integer re-reduction of the published Gram matrix.
[Paper](https://arxiv.org/html/2607.20359v4).
[Pinned P48q section data](https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/packings_out/records/dim40-delta41/P48q-91ec3c80.cert.txt).

From this directory, run

```sh
OPENBLAS_NUM_THREADS=16 python verify.py
python final.py
```

The verifier proves positive definiteness and enumerates all lattice
vectors up to squared norm four, establishing a minimum of at least
six. It then enumerates the complete norm-54 affine core from its
coefficient residue modulo three. The lattice minimum bounds every
core–core product by 27. All remaining unordered pairs are checked
using binary32 matrix products with an explicit bound below \(2^{24}\)
on every integer accumulation. The proof uses the shipped Gram form,
residue, enclosure and explicit short layers.
