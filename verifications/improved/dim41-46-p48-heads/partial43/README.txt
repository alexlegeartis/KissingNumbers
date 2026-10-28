43-dimensional kissing configuration with 2,919,782 points.

Run: python verify.py
Dependency: NumPy. The replay uses integer arithmetic for every angle.

The configuration consists of the complete families in configuration_43.npz,
the explicitly selected tiny-body/head pairs in explicit_points_43.npz, and
the five pure heads. Complete-family populations use population_note.txt.
The published parent minimum 6 and 11-design property remain inputs.

Coordinate definition in 43 dimensions:
Let q be body_complement.npz['q'] and L a real Cholesky factor, L L^T=q.
For each complete family with residue c, r3=r, and head h=z/D, enumerate
all a in Z^40 with a mod 3 = c and a q a^T = 3r. Include
  (a L / sqrt(3), h) / sqrt(r + ||h||^2).
The residue c is the corresponding row of family_residues.
For an explicit tiny body X, recover its integer coefficient a from
  a = X b^T (3 q^-1) / 288,
and use the same formula with r=10 and its supplied head.
For each pure direction p include (0_40, p/||p||).

The explicit tiny layer replaces the previous complete tiny families.
It must not be added on top of those old tiny families.

The short-difference improvement used by the verifier is exact: absence of
a norm-10/3 difference implies squared difference norm at least 16/3 in
that discriminant class. The full 30,240-vector tiny pool is checked by
verify_families.py and its completeness is explained in population_note.txt.

28 October 2026 continuation: 3,447 explicit tiny points, replacing 46
old points with 85 and gaining 39. The previous layer, exact delta,
independent coefficient-coordinate verifier and reproduction are in
../continuation43/. Complete-family population assumptions are unchanged.
