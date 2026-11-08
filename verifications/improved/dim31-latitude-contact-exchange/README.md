# Latitude-contact exchange: 251,200 points in dimension 31

8 November 2026. An explicit four-for-eight exchange improves the preceding
repository certificate from **251,196 to 251,200**. This is four more points
than the 4 November algebraic completion, and 31 more than PackingStar's
latest checked 251,169-point release. The preceding certificate is included
byte for byte. Earlier constructions and verification packages are preserved.

## Exact construction

Use zero-based coordinate indices and squared-norm 128 integer numerators:

```
A = 4(-e0 + e2 + e6 - e23)
B = 2(-e8 + e9 + e13 + e14 - e15 + e17 + e18 + e22)
C = 4(-e24 + e25)
D = 4(e27 - e29).
```

The four blocks are mutually orthogonal, with squared norms 64, 32, 32, 32.
Delete the four points **±A ±(C+D)**, and insert the eight points
**±A ±B ±C**, with independent signs. The eight insertions have squared norm 128
and mutual products at most 64. The complete proof below checks their
products with every retained point too.

The four roots ±2B and ±2C remain in the configuration. Together with the
insertions they form a scaled A3 root system. The exchange converts four
norm 64-tail latitude points into eight norm 32-tail midpoint points. It is
not a free insertion and is not an antipodal-symmetry assumption: the final
configuration has 125,570 antipodal pairs and 60 unpaired points.

The search enumerated negative-contact sums anchored at **all 6,592** points
of the old norm 64-tail latitude. It generated 1,243,881 distinct new sums;
52 had at most two incumbent blockers. Joint selection found this exchange.
The optimizer is not part of the proof. A second complete latitude pass on
the new source found no additional gain in its stated two-blocker pool.

## Reconstruct and verify

Python, NumPy, and g++ on a little-endian host suffice. No optimizer, network retrieval,
search cache, lattice minimum, or unenumerated population premise is needed.

```sh
OPENBLAS_NUM_THREADS=1 python3 reproduce.py
OPENBLAS_NUM_THREADS=1 python3 verify.py
python3 negative_controls.py
```

- `final.py` supplies the dimension and cardinality to the repository registry.
- `construction.py` reconstructs the entire witness from the pinned source
  and the four-block formula. `data/exchange31.npz` also records the delta.
- `verify_geometry.py` checks every norm, uniqueness condition, opposite
  pairing, and all unordered pairs. Its bounded integer products are exact
  in binary32; there is no tolerance. `--precision 64` selects an alternative complete binary64 replay.
- `verify.py` separately checks every inserted point against every final
  point by direct int64 coordinate sums, without BLAS.
- `independent_audit31.py` reconstructs the source/target
  difference and checks all **31,550,594,400** unordered target pairs with
  a separate C++ integer kernel, without floating arithmetic or antipodal
  compression. Its negative controls reject wrong norms, duplicate points,
  and a same-norm pair with forbidden product 112.

The independent pair counts are 31,548,584,836 retained/retained,
2,009,536 retained/added, and 28 added/added. Every category has maximum 64
and zero violations. The producer's separate inserted-versus-final pass
checks 2,009,592 directed comparisons, including each added/added pair twice.
These are complete constructive lower-bound proofs, not optimality claims.

The file manifest covers all package files except itself and Python caches.
Ordinary verification does not rewrite evidence. The wrapper `verify.py` runs
the independent checker in an isolated temporary output directory. It
automatically uses the scalar integer kernel when AVX-512F/BW is unavailable.

## Comparison and scope

The [fresh November 8 primary-source comparison](../../../common/source_review_20261108_31.json)
checks PackingStar main at 477df34, Takhanov–Yun at b376831, and Li at b823640.
Their reviewed 31D values are 251,169, 238,666, and 241,079. Cohn's leaderboard
returned 502, so the dated comparison does **not** claim exhaustive worldwide
priority. The manuscript is deliberately unchanged.

The broader campaign also tests whole exceptional-tail replacements, complete or explicitly partial
core-tetrad sweeps, every core six-support, and complete latitude/unpaired contact
closures. Its negative results are restricted to their stated families and
budgets. See `research/exceptional_tail_replacement_20261108` for research
programs and reports; none is needed to verify this certificate.
