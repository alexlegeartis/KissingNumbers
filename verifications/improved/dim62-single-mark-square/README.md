# Single-mark signed shells in 62D and 63D

| Dimension | Verified lower bound | Gain over preceding bound |
|---:|---:|---:|
| 62 | 166,921,228 | 16,359,424 (+10.87%) |
| 63 | 226,280,708 | 9,236,480 (+4.26%) |

Start from the [self-sparse Reed–Muller construction](../dim62-67-rm2-selfsparse/). Every weight-sixteen word of its 64-coordinate sign code is orthogonal over `F₂` to the entire dense sign stock. This follows because quadratic words are orthogonal to degree-three words, and the eighty cubic shifts themselves are pairwise orthogonal. The replay also checks the needed parity directly on every added support.

**The rule.** Take a weight-sixteen support meeting the original marked coordinates exactly once. Keep its fifteen ordinary coordinates with amplitude two, and replace the marked coordinate by a fixed directed head of length two. Attach a shortened extended-Hamming coset of 1,024 body signings. Choose its parity opposite the marked bit of any dense fibre that can meet that head too closely. This forces at least one sign disagreement, which pays for the extra tail product. At a common marked site, two old supports overlap on at most seven ordinary sites, so the same head is safe. At different marked sites they may overlap on eight, so assign heads with nonpositive dot products.

In **62D**, use the four square directions for the four marked sites, giving 4,032 additional supports per direction. The four marked sites are 0, 9, 20 and 41 in the six-bit coordinate space. The only dangerous dense hexagon banks share the same bit at the corresponding marked site. This uses a different set of marked sites from the earlier lower bound; it reconstructs the dense hexagon stock, the 11,929 unmarked weight-sixteen supports, 9,145 affine planes, and the poles on the new coordinates.

In **63D**, use the two opposite head directions along the collapsed axis in the exact equality section. There are 7,136 supports of each marked type. Keep the dense section, the unmarked weight-sixteen shell, the signed affine planes and the poles. Omit the old weight-sixteen shell with both marked sites. The new construction exceeds the stronger preceding 63D bound of 217,044,228; it does not assume the new vectors fit that separate construction unchanged.

All vectors have squared norm 64. No optimization is needed for replay. Run `python3 verify.py` with a C++17 compiler available; it first replays the hash-pinned parent counts and length-64 distance with independent enumerators, then checks the exact parity, head and Hamming-coset conditions. These are lower bounds, not optimality claims.
