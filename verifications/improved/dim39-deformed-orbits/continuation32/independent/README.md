# Independent complete audit: 769,436 points in dimension 39

PASS. This audits the frozen certificate at `../orbit_certificate.json`.

SHA-256: `e9f2dd59bf5169bf915be8e7d622a7ff14481471e2e0a96e0f9d90b46a5c6ad7`.

The point count is 328,192 dense + 438,272 middle + 2,908 roots + 64 stars = **769,436**. No producer verification implementation was imported or executed. No floating-point arithmetic or numerical tolerance is used for geometric acceptance.

## Dense layer

`prepare_dense.py` checks the frozen input hash, integer coordinates, and literal XOR closure of the sign-mask group. `dense_direct.cpp` computes every unordered representative pair against every explicit relative sign mask using direct signed-int64 dot products. It does not use the producer's Walsh-transform algorithm.

There are 1,683,122,416 tested relative products. Only a representative paired with itself under the identity mask is excluded. Closure of the diagonal sign group makes these tests cover every distinct pair of actual dense points. The maximum integer product is 4,999,992,742,814,284, below the sufficient threshold 4,999,999,600,000,008. Every exact squared norm is at least 99,999,996 squared, so this threshold proves normalized cosine at most one-half. The same checks rule out duplicates and nontrivial stabilizers.

Each coefficient has absolute value at most 100,000,000. The absolute accumulation bound 39 times 100,000,000 squared is below 2^63. AVX2 integer multiplication is only an acceleration; a scalar fallback is included, and 8,192 deterministic direct scalar comparisons are checked. The extremal witness is also reconstructed using Python arbitrary-precision integers in `dense_result.json`.

## Sparse and mixed pairs

`audit_sparse39.py` independently reconstructs the four actual support replacements, all retained middle sign masks, the root deletions and all heavy-star points. It checks all 5,860,176 support intersections and every sparse pair family, including actual same-support sign classes.

For every ordered support it explicitly forms the union of retained middle masks XOR projected dense-group masks. It does not assume uniform sign parity or use the producer's parity shortcut. It checks 8,848,056,320 actual signed coordinate sums. The maximum is 141,421,267; its squared margin against twice the uniform dense norm lower bound squared is 23,640,114,743. Actual retained roots and group-transformed heavy stars are independently checked too.

The explicit support and sign-mask reconstruction is stored in `actual_sparse_patterns39.json`. `sparse_dense_sparse_audit39.json` gives all counts, witnesses, arithmetic bounds and source hashes. Sparse within-family distinctness follows from unique supports/sign masks/root vectors/star vectors. The four sparse and dense shapes are disjoint; all dense coordinates are nonzero.

## Reproduction

From this directory, run:

    python replay.py

Or run `python verify.py` from the package root. Python 3, NumPy and g++ with OpenMP are required. The replay regenerates the binary input, compiles the direct verifier, runs the literal sparse/sign-pattern audit, reconstructs extreme-point witnesses and checks two corrupted binary inputs. It never imports producer geometry code or trusts archived caches.

This certifies the stated finite configuration, not optimality or worldwide priority. All dense representatives were jointly deformed; numerical count changes should not be described as literal insertion into an unchanged old dense set.

## Portable package replay

Run `python replay.py` from this directory, or `python verify.py` from the package root. Python 3, NumPy and g++ with OpenMP are required. The C++ binary and all derived inputs/pattern tables are regenerated; no archived executable or cache is trusted. Only input/output paths and the explicit no-`-O` guard in `prepare_dense.py` were adjusted in the independent source. `replay.py` adds orchestration, fresh integer witness checks, receipt comparison, and two binary-level rejection controls. `original_audit.tar.gz` preserves the complete original audit, including all source, logs, derived tables and original checksums.
