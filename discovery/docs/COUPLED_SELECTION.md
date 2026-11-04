# Exact coupled point selection

`selection.isometric_exchange` takes one explicit, independently exact verified
weighted integer spherical code. It realizes a second copy by an arbitrary
rational isometry, removes exact overlaps, and chooses a maximum cardinality
subset of the two-copy union. Both within-copy geometries remain valid; all
remaining conflicts form a bipartite graph. Selection is at individual points,
including accessory points, without orbit or product-fiber constraints.

`selection.coupled_union` performs the same optimization on two independently
verified explicit codes. They may have different sizes and alphabets. They must
use the same coordinate width and positive integral diagonal metric. Norms are
aligned by exact common integer scaling when their ratio is a rational square;
unsupported irrational alignments reject. With no map parameter, the given
coordinate frames are used. The second code can otherwise be placed by the
same rational isometry parameters as the one-input operator.

```json
{
  "schema_version": 1,
  "nodes": [
    {"id": "base", "op": "geometry.integer_code", "params": {
      "coordinates": [[-5, 0], [3, 4], [4, -3]], "norm_squared": 25
    }},
    {"id": "exchange", "op": "selection.isometric_exchange", "inputs": ["base"],
     "params": {"normals": [[1, 2]]}}
  ],
  "outputs": ["exchange"]
}
```

This exact control removes two original points, retains one, and adds all three
image points. The result has four points. It demonstrates a coupled exchange;
it is not a kissing-number record.

Exactly one of `matrix` or `normals` is required by `isometric_exchange`. Matrix
entries and normal entries are exact rational JSON scalars: integers or
`[numerator, positive_denominator]`. Rows act on the right, `image = source @ A`.
With source metric `W = diag(weights)`, the map is checked by `A W A^T = W`.
Each Householder normal `v` applies
`x -> x - 2 * (x W v) / (v W v) * v`. A list of normals applies these maps in
order; an empty list is the identity. Thus general rational block mixing and
compositions are supported, beyond permutations of existing frames. Floats,
zero normals, nonisometric maps and unsupported source representations reject.

The implementation uses a common denominator to produce integer coordinates.
Every cross inner product is exact, with a proven safe NumPy `int64` bound or
Python integer fallback. Duplicate image points are removed against the first
copy before optimization. They can safely stay in the first copy: a shared
point is already compatible with each full valid source.

SciPy finds a maximum bipartite matching. The implementation then constructs a
vertex cover and checks every matching edge, endpoint uniqueness, every cover
edge and equality of matching/cover sizes. This supplies the exact upper bound
on the selected finite union. It then independently checks all output norms,
distinctness, pairs and rank through the existing integer verifier. The output
is a plain `spherical_code`, usable by subsequent operators.

The optimization receipt is diagnostic data in the output verification report.
Ordinary saved-artifact replay independently verifies the output coordinates;
it does not infer union-optimality from an unchecked saved optimization receipt.
The 14D campaign separately includes solver-free witness replay: saturating
disjoint geometric conflict pairs certify its reported no-gain upper bounds.

Both operators preflight the full possible union's materialization, dimension
and pair budgets before source pair checking or matrix work. Row-block products
and stored graph edges are also bounded by `max_gram_entries`, and arithmetic
growth is bounded by `max_rational_bits`. Defaults support the 1,932-point 14D
control. Larger unions may need explicitly larger declared pair budgets.

These operators provide exact optimization for two internally valid explicit
components. They do not search arbitrary unknown alphabets, variable masses at
every product point, new source codes, implicit massive codes, arbitrary
algebraic isometries or joint independent sets on three or more mutually
conflicting components. Those remain gaps in the requested universal coupled
product grammar. They impose no dimension-specific family dispatch and make no
global kissing-number upper-bound or novelty claim.

A separate bounded research driver in
[`research/discovery_campaign_20261104/coupled14`](../../research/discovery_campaign_20261104/coupled14/README.md)
does build joint heterogeneous-mass and multiple-image pools with complete exact
graphs, numerical MILP search, and independently verified final coordinates.
That experiment is not yet a generic registered product-pool operator.
