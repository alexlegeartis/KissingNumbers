EXACT KISSING CONFIGURATIONS IN DIMENSIONS 41, 42 AND 44
Research checkpoint: 22 October 2026 — second pass

RESULTS

 dimension   points       families   exact family-pair conditions
 41          1,460,034       151        11,476
 42          1,989,954       181        16,471
 44          4,939,281     1,034       535,095

These are constructive lower bounds. The full standalone replay completed
successfully, including a fresh exact enumeration of the 25 distinct core
body classes used by these three configurations. Its recorded results are in
reports/results.json. The run took approximately 145 seconds with six workers
in the development environment; other machines may take substantially longer.

REPLAY

Requires Python 3, numpy, python-flint and a C++17 compiler named g++.

  python3 -m pip install -r requirements.txt
  OPENBLAS_NUM_THREADS=1 python3 run_verify.py --jobs 6

Use a smaller --jobs value on machines with fewer processors. Run without
Python's -O flag: the verifier intentionally uses assertions. The optional
--skip-recount performs the algebra checks using recorded core populations;
it is explicitly labelled as a partial replay and is not a new population
verification. All essential input data are included; no network is required
after dependencies are installed.

WHAT IS CONSTRUCTED

Let L be the pinned extremal even unimodular 48-dimensional lattice P48n.
It has minimum squared norm 6. An embedded sqrt(3) E8 lattice has orthogonal
complement K of rank 40 and determinant 3^8. The script geometry.py verifies
the embedding and that the supplied basis is the full integral orthogonal
complement, not a proper sublattice.

Each row of certificate_DIM.npz specifies:
  classes: a discriminant class in K*/K, in the compact E8 mod 3 labelling;
  r3: three times the squared radius of its body vectors;
  head_numerators z and head_denominator D: a head h=z/(D sqrt(3));
  weights: the number of body vectors used in the family.

The head arrays have eight columns, but their exact rational ranks are 1, 2
and 4 respectively. Thus the ambient dimensions really are 40+1, 40+2 and
40+4. The final points are (x,h)/sqrt(||x||^2+||h||^2).

For a fully explicit description of the bodies, let b be the integer 40 by
48 matrix in data/e8_complement_reduced.npz, and q=b b^T/96. The full shell
in a class consists of vectors

  x = w b/(3 sqrt(96)),    w in Z^40,
  w mod 3 = the class label reconstructed by geometry.py,
  w q w^T = 3 r3.

The radius-10/3 bodies are supplied directly as the verified subset in
data/dual40_min.npz, with coordinates H/(3 sqrt(96)). Completeness of this
small subset is not needed: the verifier checks its 30,240 vectors, distinctness,
norms, lattice membership, orthogonality and discriminant labels individually.

POPULATION PROOFS

The radius-6 central family belongs to a NONZERO affine class. Its population
is 1,084,590. It must not be confused with the ordinary zero-class shell,
whose population is 1,092,000.

In the new 41D and 42D constructions, the radius-16/3 outer families have
shortest normal representatives of squared norm 2/3. Their uniform population
116,235 is proved by the replayed moment identities. The 44D outer families
instead have shortest normal representatives of squared norm 8/3, and their
individual populations are enumerated exactly. These two cases must not be
confused.

For the core, run_verify.py regenerates an enclosing rational ellipsoid from
the exact Gram matrix. The C++ program enumerates every integer coefficient
vector in the required congruence class inside that ellipsoid and accepts
only vectors with the exact target norm.

More precisely, with P=2^24, integer lower-triangular M and positive diagonal
D0, let E=P^3 q-M D0 M^T. Put eta=max_i sum_j |E_ij| and T=trace(q^-1).
For w q w^T=t,

  w M D0 M^T w^T <= t (P^3 + eta T).

This follows from the spectral norm bound by eta and ||w||^2 <= T t.
The rounded floating-point Cholesky factors are therefore only used to choose
an efficient ellipsoid; the exact residual certifies enclosure. Recursive
range decisions use integer arithmetic and corrected integer square roots.
The accepted-vector test uses the original integer Gram matrix. Runtime
assertions check the integer size bounds used by the C++ implementation.

The radius-14/3 and radius-4 populations are 9,180 and 495 respectively.
verify_masses.py replays exact degree-at-most-10 polynomial identities for
these counts using the 11-design property of the parent minimal shell. It
also proves the relevant projection domain exhaustively. Abstract E8 Weyl
symmetry is used to transport polynomial identities, not as an assumed
automorphism of P48n.

ANGLE PROOF

For two body classes, let mu be the squared norm of a short normal gluing
representative for their difference. The parent minimum gives

  ||x-y||^2 >= 6-mu,
  <x,y> <= (r_i+r_j-6+mu)/2.

For the same class and radius at two different heads, the verifier instead
uses <x,y> <= r_i, to cover the possibility x=y. Within one family, distinct
bodies have squared separation at least 6, giving <x,y> <= r_i-3.

Write this chosen body inner-product bound as M_ij. Define integer quantities

  N_i = r3_i D^2 + z_i.z_i,
  P_ij = (6 M_ij) D^2 + 2 z_i.z_j.

The normalized inner-product condition is exactly certified by

  P_ij <= 0, or P_ij^2 <= N_i N_j.

All these comparisons use exact integers. Every pair i<=j is checked,
including internal family conditions. Passing them proves a spherical code
with maximum inner product at most 1/2, hence the stated kissing bound.

POPULATION ACCOUNTING

 body squared radius       41D points    42D points    44D points
 6                         1,084,590     1,084,590     1,084,590
 16/3                        232,470       697,410     2,774,520
 14/3                        110,160       165,240       881,280
 4                            31,680        41,580       190,575
 10/3                          1,134         1,134         8,316
 TOTAL                     1,460,034     1,989,954     4,939,281

PUBLISHED INPUTS AND PROVENANCE

This is an exact construction and count replay conditional only on the
identified, published parent lattice facts: P48n has minimum 6, 52,416,000
minimal vectors, and its minimal shell is an 11-design. It does not reprove
the minimum of P48n by a 48-dimensional shortest-vector enumeration.

Parent lattice entry:
https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ48a.html

G. Nebe, Boris Venkov's Theory of Lattices and Spherical Designs,
Theorem 4.1 (all nonempty shells of extremal even unimodular lattices in
dimension divisible by 24 are 11-designs):
https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf

Current comparison table:
https://cohn.mit.edu/kissing-numbers/

The pinned parent embedding, complement basis, small-body vector pool and
polynomial identities were recovered from the user's earlier 43D exact
certificate/research archives. geometry.py and verify_masses.py are adapted
from those verification scripts. The new head configurations, class
relabelings, complete core recounts and standalone assembly were produced in
this research pass. The parent basis and Gram hashes are checked in code.

The older archive names are kissing43_2878155_exact_certificate.zip and
kissing43_structural_research_v4.zip. No unpublished later certificate is
needed to replay these three new constructions.
