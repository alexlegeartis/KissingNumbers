# Coupled character families

Two search operators generate and replace complete algebraic families. They
accept a binary generator matrix, integer coordinate shapes and positive
integer metric weights. All shapes must have the same positive squared norm.
There is no dimension-specific construction table or target count argument.

For a message `m` and a coordinate-column character `c_j`, the action sends
`x_j` to `(-1)^(m dot c_j) x_j`. For shapes `a,b`, collect the coefficients
`weight_j * a_j * b_j` by character. An exact Walsh transform supplies

```
dot[a,b,m] = sum_j weight_j * a_j * b_j * (-1)^(m dot c_j).
```

This checks all interactions of whole families without enumerating every
point pair. The generator rows are independently checked over F2. Bounds on
the sum of absolute coefficients protect every transform intermediate;
large coefficients use Python integers. Dense kernel allocation and transform
work have explicit budgets.

`search.character_families` assigns one shape or an empty choice to each coset
of a supplied binary quotient. A selected choice contributes the entire sign
orbit of the quotient kernel. For each pair of cosets the maximum over every
message in their difference coset determines compatibility. Nontrivial sign
stabilizers are rejected in this variant. A checked initial assignment can
be replaced completely. Optional `active_cosets` restrict a joint numerical
MILP to a large neighborhood while freezing its complement. This is a search
restriction, not a theorem about all possible constructions.

`search.character_orbits` assembles any number of whole sign orbits. It allows
sparse shapes and counts their exact stabilizers. Duplicate images within
one orbit are counted once; overlapping distinct proposed orbits conflict.
Orbit sizes are objective weights. Roots, signed supports and dense caps can
all appear as independent selectable families, so every incumbent base and
cap family can be removed in one coupled model. It does not impose one shape
per coset.

For example, this program replaces an entire two-point axis orbit with a
four-point rectangular orbit:

```json
{
  "schema_version": 1,
  "nodes": [{
    "id": "sphere", "op": "search.character_orbits",
    "params": {
      "generators": [[1,0],[0,1]],
      "shapes": [[25,0],[15,20]],
      "initial_shapes": [0], "seconds": 1
    }
  }],
  "outputs": ["sphere"]
}
```

The axis orbit has a nontrivial stabilizer. Every new point conflicts with
an old point, so insertion into a fixed incumbent would fail. Joint assembly
finds the four-point family and the existing coordinate checker independently
replays its norms, distinctness, every pair and rank.

Numerical MILP status and dual bounds remain diagnostics. The operators keep
a verified incumbent if optimization finds no better one, then independently
verify the materialized output with the existing exact integer backend. A
returned exact configuration supplies a lower bound, with no optimality or
novelty conclusion. Finite catalogues, symmetry restrictions and sampling
remain visible in research results.

`discovery.numerical_shapes.orbit_shape_objective` supplies a separate continuous
proposal mechanism. It changes all supplied shapes simultaneously and includes
every nonidentity sign operation and every retained-base interaction. The base
must have exact integer coordinates and its full sign-action closure is
checked before symmetry reduction. The analytic gradient is checked against
finite differences. Its loss, worst cosines and coordinates remain numerical;
none can certify a lower bound.

The first competitive experiments are in
[`research/character_families_20261104`](../../research/character_families_20261104/README.md).
The17D mixed-amplitude orbit pool did not improve5730. Joint19D and21D shape
deformations failed feasibility. A separately verified rational variant of
the known19D12460 construction has a strict auxiliary angular gap; its count
is a control, not a new record. No new SOTA has been found in this batch.
