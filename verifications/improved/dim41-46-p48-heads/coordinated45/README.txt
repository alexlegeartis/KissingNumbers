45D coordinated family movement, 2026-10-24

Retained count: 7,655,056, improvement 990 over the previous project bound 7,654,066.
The cumulative gain over the starting 7,653,571-point construction is 1,485.
This is a lower-bound construction, not an optimality proof or a literature-wide record audit.

Run: python verify.py
Dependencies: Python 3 and NumPy. No optimizer or network access is used by the verifier.

The configuration contains 168 large families (9,180 points each),
778 medium families (495 points each), 1,200 small families (18,900 points total),
the original 41 core families, and 16 pure-head points. Exact rational heads lie in a
rank-five rational subspace. The body space has dimension 40, so total dimension is 45.

For family (c,r), body coefficients a run through integer rows with
a mod 3 = body_residues[c] mod 3 and a Q a^T = 3r.
Their search-unit Gram products are a Q b^T / 3.
coordinates.py gives the exact coordinate definition and a numerical exporter.

New body classes: [(77, 12), (202, 12), (239, 12)]. added_bodies.npz explicitly lists their body vectors.
The verifier checks their residues, norms, distinctness, population and within-family
body angles. It also reconstructs the parent embedding, the full rank-40 body
complement, the gluing map, all 30,240 small bodies and their population table.
Existing large-family populations use the included population_note.txt. The parent's
published minimum and design property remain mathematical inputs. This bundle does
not freshly enumerate every old large body family or re-prove those published inputs.

For different families, twice the body inner product is bounded by
r+s-18+mu(c-d). Repeated copies of the same body shell use the bound 2r.
Within one family, the bound is 2r-18. After adding head inner products,
the verifier checks each angle using integer squaring with the sign retained.
It also checks pure-head angles, the exact span projector and all family populations.

The search first replaced one medium family by two and moved large and medium heads
together. A second joint optimization restored every small family. Candidate searches
and rejected numerical trials are outside this bundle.
