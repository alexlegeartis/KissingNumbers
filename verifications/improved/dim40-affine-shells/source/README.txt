1,126,314 KISSING POINTS IN DIMENSION 40
26 October 2026

This explicit construction proves tau(40) >= 1,126,314.

It improves the earlier local construction of 1,092,072 by 34,242 points,
and the published benchmark of 1,092,000 by 34,314 points.

The published comparison is Table 2 of Dorofeev, Sun and Wang,
"Optimal Extensions of Cross-Sections: Sphere Packings in Dimensions
38 to 43", arXiv:2607.20359v4, 29 August 2026, checked on 26 October:
https://arxiv.org/html/2607.20359v4

REPLAY

Requires Python 3.9+, NumPy, and g++ with C++17 and __int128 support.

    python verify.py

The program compiles the minimum checker itself, checks the supplied data,
checks every cross-family pair and every pair of auxiliary points, and writes
replay/verification.json. It requires no network or external lattice data.
It uses moderate matrix tiles to limit memory. The completed default replay
took 193.4 seconds in the development environment; hardware and BLAS affect
timing. Its full output is in verification.recorded.json.

The pair products are evaluated as EXACT bounded integers represented in
IEEE binary32. The verifier proves an absolute accumulation bound below
2^24 before each matrix multiplication. Every operand, product and partial
dot-product sum is consequently represented exactly. There is no angular
tolerance, approximate threshold, or rounded geometric acceptance test.
All norms, coset checks, and the minimum proof use integer or rational
arithmetic. See PROOF.txt for the arithmetic justification.

An alternative, substantially slower mode uses int64 for all pair products:

    python verify.py --integer-pairs

CONSTRUCTION

Let Q be gram.json and use <a,b>_Q = a Q b^T as the Euclidean inner product
on R^40. The unit vectors, specified by integer coefficient rows, are:

    source                         points       Q-norm   normalization
    core.npz: coefficients       1,084,590          54    divide by sqrt(54)
    auxiliary.npz: norm36            40,590          36    divide by 6
    auxiliary.npz: norm30             1,134          30    divide by sqrt(30)
    TOTAL                        1,126,314

For ordinary orthonormal coordinates choose R with R R^T=Q and multiply
each row by R before the indicated normalization. Existence follows from
the exact positive-definiteness check; computing R numerically is unnecessary
for the proof.

The core coefficients all lie in one affine class modulo 3. Their pairwise
differences are three times nonzero integer lattice vectors, which proves
their mutual angle bound from the minimum of Q. All other pair types are
checked directly, including every auxiliary vector against every core row.

The gain comes from replacing the original zero-class minimum shell with
a shifted shell. This loses 7,410 core points relative to 1,092,000, but
admits 40,590 norm-four points and 1,134 short-dual points, yielding a net
gain of 34,314 over that benchmark.

DISCOVERY AND PROVENANCE

Q is the same rank-40 P48n / sqrt(3) E8 complement used in the earlier
1,092,072 certificate. The construction uses affine class 199 in the
workspace's E8 mod-3 labeling. The 40,590-point layer consists of 82
495-point families: the center class and 81 root-neighbor classes.
The 1,134-point layer uses 72 further root-neighbor classes.
discovery_families.npz records those labels and doubled E8 head offsets.
These labels explain the search; the explicit coefficient arrays define
the construction and suffice for verification.

The shifted core was freshly enumerated with exact branch bounds. Eighteen
previously missing 495-point families were freshly enumerated; the remaining
64 were copied from archived explicit arrays and checked anew. The final
proof needs only the supplied vectors, not completeness of these searches.
It does not assume the extremality or published minimum of P48n.

The source lattice data ultimately derive from kP48nBasis at
https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/sphere_packing/lattices/p48nm.cc

FILES

    gram.json                   exact 40 by 40 Gram matrix
    core.npz                    all 1,084,590 core coefficient vectors
    auxiliary.npz               all 41,724 auxiliary coefficient vectors
    enclosure.npz               exact minimum-enumeration enclosure
    exact_minimum.cpp           integer-only minimum enumerator
    verify.py                   complete independent replay
    PROOF.txt                   proof and arithmetic details
    verification.recorded.json  completed packaged-verifier output
    discovery_cross_check.json  earlier independent cross-check output
    discovery_families.npz      optional geometric discovery metadata
    requirements.txt           Python dependency

This establishes a kissing-number lower bound, not an optimality result
or a sphere-packing density improvement.
