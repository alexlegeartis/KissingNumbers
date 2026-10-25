# Three affine lattice shells in dimension 40

Checkpoint: **26 October 2026**. Exact lower bound: **1,126,314**.
Previous method: Dorofeev–Sun–Wang, [Table 2, arXiv:2607.20359v4](https://arxiv.org/html/2607.20359v4), **1,092,000** points; gain **34,314**.
The inherited lattice basis comes from [wcgbg/sphere-packing](https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/sphere_packing/lattices/p48nm.cc).

1,084,590 norm-54 coefficient rows + 40,590 norm-36 rows + 1,134 norm-30
rows = **1,126,314** normalized points in the positive definite rank-40
Gram metric of determinant 6,561. `source/PROOF.txt` supplies the proof;
`source/README.txt` supplies construction details and provenance.

The core is one affine class modulo 3. An independent exact lattice
enumeration rules out nonzero vectors of squared norm at most 4. Evenness
therefore gives minimum at least 6, proving every core product bound.
All five remaining pair types are explicitly checked. Each integer
accumulation is proved below 2^24 before exact representation in IEEE
binary32; no tolerance is used. The standalone verifier also offers
`--integer-pairs` for a slower int64 replay. The lattice enumeration uses
bounded C++ __int128 arithmetic and an exact rational enclosure.

```bash
python final.py
python verify.py
```

Requires Python, NumPy, and a C++17 g++ compiler supporting __int128.
The wrapper verifies the complete SHA-256 manifest and replays every check
in a temporary directory, including compilation and the lattice minimum
proof. All coefficient rows, the Gram matrix, enclosure, C++ source,
standalone verifier and original receipt are supplied. No parent lattice
minimum assumption, optimizer, download, or search cache is needed.
