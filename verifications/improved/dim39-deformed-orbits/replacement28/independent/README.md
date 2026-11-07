# Independent exact audit: 769,820 points in dimension 39

This package reconstructs and audits the frozen integer certificate with SHA256
`dbc293625abdbc7f4377ee94e77344967877864ea7717482e34479ce794085df`.
It contains a separately written Python reconstruction and C++17 direct integer
verifier. No code from the producer, old orbit verifiers, NumPy, transforms,
optimizer, solver, or saved pair-product tables is used.

## Reproduction

Requirements: Python 3.10 or later and a C++17 compiler such as GCC or Clang.
Only Python's standard library is used. There are no downloads. Run from any
working directory:

    python3 /path/to/this/package/verify_independent.py

Set `CXX=clang++` or use `--compiler clang++` to choose a compiler. The verifier
uses one compute thread; when a process has multiple CPUs available it excludes
the highest numbered one, following the repository's resource guidance. It
builds a temporary executable and regenerated text input in `_build/`, then
writes `independent_report.json`, `direct_results.json` and `direct_progress.log`.
The original complete run is retained in `run.log`.

All calculations affecting acceptance are integer calculations. Python uses
unbounded integers. C++ uses standard signed 64-bit integers, with a checked
coordinate bound of 31,000,000. Every dense dot is bounded in absolute value by
39 × 31,000,000² < 2^63. Sparse cross dots have explicit overflow guards before
squaring. No dense dot is squared, and no compiler-specific 128-bit type is needed.

## Configuration and count

The sparse layer is reconstructed from the original 3,424 binary support rows,
the four checked support replacements, and the four heavy stars in
`inputs/base_certificate.json`. Its count is:

- 3,424 distinct supports, each carrying all 128 even-parity sign choices at
  magnitude 2: **438,272** vectors.
- All four sign choices of magnitude 4 on every coordinate pair except the 14
  distinct center/leaf pairs removed by the stars: **2,908** roots.
- Four stars, each with center sign ±4 and all eight even-parity leaf patterns at
  magnitude 2: **64** heavy vectors.

Every sparse vector has squared norm 32; sparse count is **441,244**.

The dense layer consists of **10,268** integer representatives and all **32**
coordinate-sign translations in the explicit XOR group, or **328,576** vectors.
Each dense vector is normalized by its own exact norm. All dense coordinates
are nonzero. Total: **769,820**.

## Complete pair coverage

**Sparse/sparse.** Every unordered support pair is checked directly for
intersection size at most four. Within one support all 128 even patterns are
pairwise checked for dot at most 16. Roots are generated distinctly, and every
root/root pair is checked directly. Root/middle dots are bounded by two shared
coordinates, hence at most 2 × 4 × 2 = 16. All heavy/heavy and heavy/root pairs
are checked directly. For every heavy vector and every middle support, every
one of the 128 even sign patterns is evaluated directly; all products are at
most 16. Thus sparse/sparse normalized products are at most 1/2.

**Dense/middle.** On a support, a sign-group element changes the sparse even
parity to the parity of its projection onto that support. The verifier checks
all group projections. If any is odd, the permitted union is all 256 patterns;
otherwise it is exactly the 128 even patterns. For each representative and
support it explicitly generates all 256 eight-sign sums by Gray-code additions,
accepting exactly the permitted parity classes. This is direct enumeration,
not the producer's sum-of-absolute-values/parity-correction formula.
If M is the largest accepted sign sum, the unnormalized sparse product is 2M.
The verifier requires M² ≤ 2L², with L = 99,999,996, and independently checks
that every dense norm is at least L. Therefore every dense/middle normalized
product is at most 1/2.

**Dense/roots.** For every representative and every retained root coordinate
pair, all four signs have exact maximum 4(|v_i|+|v_j|), invariant under the
sign group. The verifier checks its square against 8L².

**Dense/heavy.** Every representative, every sign-group translation, and all
64 actual heavy vectors are tested directly. Every positive dot d satisfies
d² ≤ 8L². Negative dots already satisfy the required inequality.

**Dense/dense.** XOR closure, identity, range and uniqueness of the group masks
are independently checked. By common sign action, every pair reduces to
(v_i, g v_j). The C++ verifier evaluates the direct 39-coordinate integer dot
for every i < j and every g, plus every i = j and every nonidentity g. It skips
only an actual vector paired with itself. This covers

32 × 10,268 × 10,267 / 2 + 31 × 10,268 = **1,687,063,204** relative products.

Each product is required to be at most floor(L²/2). Since each vector's norm
is at least L, this proves every normalized dot is at most 1/2. XOR group
reduction is exact; no sampling or distance estimate is used.

## Distinctness

Within the middle layer, distinct supports and distinct sign patterns give
unique vectors. Root and heavy uniqueness is explicitly checked. Nonzero
support sizes 8, 2 and 5 separate the three sparse families. All dense points
have 39 nonzero coordinates, so cannot equal any normalized sparse point.
Every nonidentity sign translation changes a nonzero dense coordinate; more
generally the complete dense pair inequality rules out all repeated normalized
rays, whether or not the integer representatives have different norms. The
count therefore counts distinct unit vectors.

## Comparison, attribution, and limitations

`remote_baseline.json` records the independent check of origin/main
`daa3a9a54ccbe450a08a4de03bed6ef1c7a43cd6`: its README, RESULTS and active 39D
verification all give **769,436**. This candidate is **384 above that pinned
remote**, and **128 above the local 769,692 predecessor** named in its provenance.
The remote can change after the observation; publication must check it again.
This audit is not a worldwide literature search, an optimality proof, a claim of
external acceptance, or evidence that the result has already been published.

The geometric proof is a complete finite-certificate proof without an inherited
lattice minimum or unenumerated population assumption. Ordinary trust in the
Python interpreter, C++ compiler, integer hardware, and file-hash implementation
remains. This is not a proof-assistant formalization.

The original support data are Qihong Lin's, used under CC BY 4.0; original
notices are preserved in `inputs/THIRD_PARTY.txt`. Only the support data from
that notice are needed for this independent certificate, since the old binary
dense code is entirely replaced. Four support modifications, sparse root/star
replacement and the new integer dense directions are explicitly supplied.

## Independent witnesses and rejection controls

After the full replay, run without Python's `-O` option:

    python3 check_witnesses.py
    python3 negative_controls.py

`check_witnesses.py` resolves all four reported dense-related maxima to actual
legal point pairs and verifies their angular inequalities again using Python's
arbitrary-precision integers. It does not substitute for global enumeration.

`negative_controls.py` tests the direct C++ engine past the outer checksum gate:
zero norm, duplicate normalized dense direction, a norm-valid dense point that
collides with an actual middle support, and a false representative count must
all be rejected. Mutants are created only in temporary files. The original
certificate and all audited inputs remain unchanged. Reports are
`witness_checks.json` and `negative_controls.json`.

## Recorded result

The complete direct audit passed. Its maximum dense dot was
**4,999,995,789,376,822**, below the permitted **4,999,999,600,000,008**, at
representative indices 188 and 10,261 and relative group mask 277,388,594,648.
All **8,858,408,960 permitted middle sign sums** passed; their maximum was
141,421,309 and their conservative square slack was **11,760,726,551**.
All **21,028,864 dense/heavy direct checks** passed. The reported root and heavy
maxima are realized by legal point pairs, independently confirmed in Python.
See the JSON reports for every count, extremum, witness and rejection control.
