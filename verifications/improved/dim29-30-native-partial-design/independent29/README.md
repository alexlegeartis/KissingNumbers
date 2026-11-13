# Independent audit: K(29) ≥231814

This separately written verifier freezes the proposed coordinate artifact and checks it without assuming the source configuration is valid and without antipodal grouping.

Run `OPENBLAS_NUM_THREADS=2 python audit.py`. It verifies:

*231814 distinct integer vectors in29dimensions;
*every squared norm equals128;
*every unordered pair has inner product≤64;
*exact exchange from the bundled231776-point source:492 genuinely new points added,454 existing points removed;
*the rebuilt source exchange equals the full candidate set exactly.

The normalization metadata describes coordinates divided by2. Thus the geometric squared norm is32 and the allowed inner product16; equivalently the integer array uses128 and64. Normalizing gives unit vectors with pair products≤1/2, proving a kissing-number lower bound.

Arithmetic: squared norms and structural checks use int64. Pair checks use ordinary CPU binary64 matrix products on integer inputs. Each individual product and any partial sum is bounded by29 times the maximum absolute coordinate squared, which is far below2^53. Every multiply/add is therefore exact; no floating tolerance is used. The checker visits all n(n−1)/2 unordered pairs directly.

The receipt records exact file hashes, pair count, maximum product and elapsed time. This audit establishes the construction, not worldwide priority or optimality.

Data are read from ../data; no coordinate arrays are duplicated here. The verifier is independent29/audit.py. Run from any directory using its path.

Frozen full coordinate file SHA256:
c237076120b06c6d9f928ef91da410e972a9299f9cf8de5cb986015a463d8b78

Source: KissingNumbers dim28-29-november-sextets coordinates29_231776.npz. Proposed construction and delta supplied by the native partial-design search in section30_newdesign. No source files were modified by this audit.
