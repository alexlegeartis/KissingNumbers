Independent sparse and dense/sparse audit: 769436 points in dimension 39

STATUS: PASS for the scope listed below. Dense/dense separation and dense orbit distinctness are covered by the separate dense audit in this directory, not by audit_sparse39.py.

Reproduce from the workspace root:
  OPENBLAS_NUM_THREADS=1 python research/independent_audit39_769436_20261106/audit_sparse39.py
  python research/independent_audit39_769436_20261106/audit_sparse_witnesses39.py --size 769436

Independent implementation:
  audit_sparse39.py imports no producer code. It reads certificate JSON and source text directly. All acceptance arithmetic is exact: int64 dot products with formally checked overflow bounds, and arbitrary-precision Python squares. No floating-point or numerical tolerance is used for acceptance.
  audit_sparse_witnesses39.py resolves each reported extremal dense/sparse product to an original sparse vector and a group-transformed representative, recomputing products with pure Python integers.

Complete scope:
  - Frozen orbit certificate hash, source hashes, dimension, integer coordinates, group closure and representative norm bounds.
  - Original source support decoding; all four support edits verified against their original rows; support uniqueness.
  - Every actual retained middle mask, projected group mask and reachable-mask union, with exact parity counts.
  - Every one of the 5,860,176 distinct support intersections. Maximum intersection is four.
  - Every actual same-support mask class. Minimum sign Hamming distance is two.
  - All 2908 roots after deleting all four signs on exactly 14 center-leaf pairs, and all 64 actual heavy-star vectors.
  - All root/root, star/star and root/star dot products; all actual middle/star patterns; exact middle/root maxima from root-sign completeness. Every sparse pair family has inner product at most 16, with squared norms 32.
  - Every dense/middle pair by literal reachable-mask enumeration: 8848056320 exact signed sums.
  - Every dense/root pair by the sign-invariant full root-sign set: 7456112 exact coordinate-pair checks.
  - Every dense/star pair by the literal union of group-transformed star vectors: 1312768 exact products.
  - Sparse distinctness and dense/sparse distinctness. Dense orbit distinctness remains the separately documented dense audit's responsibility.

Exact results, using L = 99999996:
  Middle M = 141421267; M^2 = 19999974759885289 <= 2 L^2 = 19999998400000032; slack 23640114743.
  Root M = 49454881; 2 M^2 = 4891570509448322 <= L^2 = 9999999200000016; slack 5108428690551694.
  Half-star M = 128757197; M^2 = 16578415779296809 <= 2 L^2 = 19999998400000032; slack 3421582620703223.

Counts:
  Middle 438272; roots 2908; stars 64; sparse total 441244.
  Dense representatives 10256; sign group order 32; dense total 328192 provided dense distinctness passes; total 769436.
  Retained-mask count distribution: {'128': 3424}.
  Reachable-mask count distribution: {'128': 108, '256': 3316}.

Files:
  sparse_dense_sparse_audit39.json: scope, complete numerical results, and SHA-256 hashes.
  actual_sparse_patterns39.json: every actual support, retained mask list, projected mask list, reachable mask list, all roots and stars, and per-support dense maxima.
  actual_extreme_pair_witnesses39.json: independently reconstructed actual extremal point pairs.

Reachable masks are only a dense/sparse verification device. They do not add sparse points or replace the actual retained sparse mask lists. The source binary dense kernel is hashed for provenance only; it is not used to prove the new dense configuration.
