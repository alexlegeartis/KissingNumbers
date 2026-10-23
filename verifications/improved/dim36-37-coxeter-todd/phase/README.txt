Verified kissing configurations: 501633 points in 36D, 505359 in 37D
24 October 2026

Run the exact lower-bound verifier:
    python verify.py

Python 3.9+ is sufficient; no third-party libraries or search are needed.
It writes verification.json. All checks passed in the recorded run.

Run the optional exact section-obstruction check:
    python verify_sections.py

This additionally uses NumPy. It writes section_obstruction.json.

Read proof.txt for the construction and mathematical proof.
Read research_directions.txt for the tested and proposed next steps.

certificate.json defines both configurations.
exact_lattice.py supplies rational exhaustive lattice enumeration.
vectors36.npz and vectors37.npz are optional explicit coefficient lists.
coordinates_check.json records exact norms and distinctness of those lists.

COORDINATE CONVENTION
Rows are integer coefficients in a Gram basis, not Cartesian coordinates.
Their squared norm is x @ gram @ x.T. Unit Cartesian coordinates are
coordinates @ numpy.linalg.cholesky(gram) / sqrt(12).

The constructions are locally verified and ready for independent review.
They have not been submitted, published, or added to an external registry.
Benchmark values refer to the Cohn table snapshot from 23 October 2026.
