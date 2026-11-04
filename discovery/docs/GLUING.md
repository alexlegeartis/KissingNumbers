# Exact finite Gram gluing

`gluing.union` accepts two independently exact sources, each either explicit
weighted integer coordinates (`spherical_code`) or a rational Gram matrix
(`gram_spherical_code`). Sources can use different representations and different
norms. Implicit code or orbit descriptions, floating-point coordinates, intervals,
and general algebraic-number backends are explicitly unsupported here.

The required `cross_gram` is the full matrix of **unscaled** inner products
between the source point lists in a proposed common base space. It is supplied
construction data. The operator does not infer a relative frame from unrelated
source coordinates or assume that a bound on cross angles is sufficient. It
reconstructs each source Gram and checks the full joint matrix
`[[G_A,C],[transpose(C),G_B]]` for exact positive semidefiniteness. This certifies
a simultaneous real realization of both sources and all specified cross values.

Optional `scales` has two exact nonzero rational values; negative values are
allowed. **Zero scales are rejected as unsupported noninjective maps.** This
backend does not claim to handle all scalings. Optional `heads` contains two
rational vectors of the same width, in a shared space orthogonal to the entire
base realization. Scalars are integers or `[numerator,positive_denominator]`;
floats are rejected. Defaults are unit scales and empty head vectors.

The mapped points are `(s_A*x,h_A)` and `(s_B*y,h_B)`. Their within-source
products are `s_A^2*G_A + dot(h_A,h_A)` and similarly for B; cross products are
`s_A*s_B*C + dot(h_A,h_B)`. The two common squared norms must agree. All signs
are evaluated exactly, including negative scales. A source's input coordinate
width is recorded separately from the union's exact Gram rank. A head or a
nonorthogonal relative frame can increase rank; the operator never substitutes
the maximum input dimension for the computed union dimension.

Joint PSD is checked **before** deduplication. In a PSD common-norm Gram,
`g_ij=norm` is exactly the zero-distance condition. The operator computes its
equivalence classes in deterministic source/list order and keeps one point from
each class. Partial overlaps and identical components therefore contribute only
their exact set-union cardinality. An inconsistent fake overlap cannot hide
behind deletion of a row: the original joint PSD check rejects it first.
The deduplicated Gram then undergoes the shared exact norm, all-pair threshold,
PSD, and rank checker. Every distinct-pair product must be at most half the
common squared norm; no numerical tolerance or unsafe envelope is accepted.

The Python API in `discovery.compatibility` exposes `within(source, context,
scale=..., head=...)`, `between(A,B,cross_gram,context,scales=...,heads=...)`, and
`certified_union(...)`. The first two return independently recomputed
`Verification` reports; the last returns the replayable Gram artifact.

The output embeds its full validated source artifacts, their content digests,
maps, provided cross values, and overlap classes under `data.gluing`. Verification
replays all source mathematics and reconstructs the output Gram, rank and count.
Saved exact-status fields confer no authority. Source provenance is retained as
metadata; historical, priority, or novelty claims are not certified by gluing.
Nested gluing is permitted within the explicit evidence-node limit.

Point count, pre-dedup pair count, Gram-entry count and head-vector budgets are
checked before joint arrays or source work. Integer-coordinate-to-Gram work,
rational bit lengths, and PSD/rank work also have explicit bounds. Input
coordinate width and output rank are separately subject to dimension budgets.
These are finite exact backends, not a claim of a general scalable 96D gluing
solver, an implicit-orbit compatibility oracle, or global optimality.

Run the mixed coordinate/Gram rank-three control with:

```bash
python -m discovery run discovery/recipes/nonorthogonal_gluing.json --out runs/gluing.json
python -m discovery verify runs/gluing.json
```
