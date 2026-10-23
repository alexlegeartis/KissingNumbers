36D kissing number: an exact-verified 501,147-point construction
=============================================================

This package supplies a finite construction giving tau(36) >= 501147.
Against Cohn's 23 October 2026 table entry 484760, the gain is 16387
points (approximately 3.38 percent).

Run with Python 3.9 or newer:

    python verify.py

The verifier uses only Python's standard library. It enumerates the short
Coxeter--Todd lattice shells with rational arithmetic and checks every
pair of the 5038 supplied outer cosets. It requires no search software,
optimization solver, downloads, or floating point geometry.

Files:

  certificate.json    Gram matrix, sublattice, label map, all outer cosets.
  verify.py           Standalone exact verifier.
  verification.json   Verification output, including certificate SHA-256.
  proof.txt           Mathematical proof, discovery method, and sources.
  discovery.json      Algebraic core and the 81 extra cosets.
  vectors.npz         Optional explicit integer coefficients and Gram matrix.
  coordinates_check.json  Checks on the auxiliary coefficient export.

All vectors have squared norm 12; all distinct pairs have squared distance
at least 12, hence inner product at most 6. Normalization gives unit vectors
with pairwise inner products at most 1/2.

Important coordinate convention: vectors.npz stores coefficients in the
supplied Gram basis. For ordinary Cartesian unit coordinates, use

    coefficients @ numpy.linalg.cholesky(gram) / sqrt(12)

The proof and standalone verifier do not depend on this floating conversion
or the optional NumPy coordinate file.

The structured core has 498960 points. Its extra-block problem separates
into 252 small independent-set problems after imposing a + b A + c = 0
over F_3^6. The final certificate adds 81 more blocks of 27 points.

Appending a zero coordinate and adding the two poles immediately gives
tau(37) >= 501149, also above the retrieved 37D table entry 500008.

This is a local computation and proof package, not a claim that an external
leaderboard has been updated or that third-party review has occurred.
