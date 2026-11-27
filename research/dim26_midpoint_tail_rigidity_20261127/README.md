# A fixed-body obstruction for the 26D Takhanov midpoint bank

**No new lower bound.** We use [Takhanov's exact 201,212-point v0.2 configuration](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2), including the nonintegral section. All rows have squared norm 128 and distinct rows must have inner product at most 64.

**Theorem for this stock.** Keep the equator and all non-midpoint rows fixed. Keep all 4,416 existing midpoint bodies, but allow both tails of every antipodal fibre to rotate arbitrarily in the two-dimensional tail space, with tail squared norm 32. Then no new midpoint body obtained from an edge of the deleted Leech shell can support even one additional ray with tail squared norm 32. This excludes arbitrary real tails in the stated joint deformation, not just the two old lattice directions.

## Candidate bodies and the two-colour picture

A missing body is `b=(u+v)/2` for two deleted shell points with `u·v=64`. It has squared norm 96 and is automatically compatible with the retained equator: each retained Leech row has inner products at most 64 with both endpoints. All 4,480 removed shell rows give 367,270 distinct missing bodies of this type after discarding existing bodies. There are 8,832 old lifts forming exactly 4,416 antipodal pairs. Their original tail lines are `A=(4,4)`, `B=(4,-4)`. For old lines the paired bodies with body dot greater than 32 must have opposite colours; when we keep two orthogonal colours, a legal change can swap colours only on entire connected components of this graph. At frozen inputs it has 18 components, including the forced 4,096-body component.

For candidate bodies, 74,278 survive the necessary comparison with all existing paired bodies, `b·b_j≤64`. Exactly 332 candidate/body-tail lines satisfy the frozen comparisons on both signs; every such attempt requires inconsistent flips. Allowing an arbitrary real single target tail while keeping the old palette (with any permitted component flips) leaves 256 candidates under the giant component's two exact radius inequalities. All 7,626 algebraic boundary checks fail.

## Continuous rotations of the existing bank

When `b_i·b_j=64`, their two tail lines must be orthogonal even after arbitrary real rotations. These edges give 1,082 components, each described by one orientation modulo pi and the inherited bipartition. On 954 of them, the integer frozen one-axis caps and orthogonality force both absolute coordinates of the common line to be at most 4. Squared tail norm 32 therefore forces precisely the two diagonal lines.

For each of the 128 other components, actual frozen comparisons prove that its angle from the original line is strictly less than `theta0=arctan(1/8)`. In a coordinate frame with that original angle zero, a comparison with a frozen row becomes

```
|4 (alpha cos(theta) + beta sin(theta))| <= h.
```

For each sign of rotation, an exact comparison outside the endpoint interval witnesses strict violation both at the appropriately signed `theta0` and at the corresponding right angle. After fixing signs of alpha and beta this expression is positive and concave on the intervening first-quadrant angular interval. Thus an entire half of the complement is excluded by each witness, not sampled. All 128 soft components are in the same giant 4,096-body bank.

The remaining old cross-colour edges have body dot 48. Hard–hard edges still force the same flip, and any hard–soft edge forces the hard component to keep its original orientation: with a mismatched flip the absolute tail dot is larger than 16. Between soft components, possible deviations add to less than 30 degrees, so these old edges are unaffected. Thus all hard positions in the giant bank stay original, and hard-only groups outside it have just the old component flips.

- If a new tail is itself on a diagonal line, restoring every soft line to its old position cannot hurt. Against an old same-colour soft body, the tail dot under any permitted change exceeds 31.5, so even the smallest positive half-integer excess of body dot above 32 cannot be rescued. Against the opposite colour, restoring the old orientation only decreases the absolute tail dot. This reduces to the palette replay above.
- If a new tail is not diagonal, its body dot with every hard body is strictly less than 64. At 64 it would need to be orthogonal to a diagonal and thus diagonal. 9,088 bodies meet that test; the two hard giant caps leave 240. To overestimate the help from all soft motions, relax each soft comparison against its **original** colour by 4 in the raw dot scale: the possible change of a tail dot is bounded by `32 * 2 sin(theta0/2) < 4`. Keeping the exact hard flips and all frozen comparisons, a second circle-boundary check tests 14,800 points. Even this relaxation has no surviving ray.

## Replay

Place the official pinned [`d26_certificate.npz`](https://github.com/k-nic/KissingNumbers/blob/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/data/d26_certificate.npz) and [`common.py`](https://github.com/k-nic/KissingNumbers/blob/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/common.py) beside these scripts. The NPZ SHA-256 is `35fd8143cf0c7e6b0146719485e00c86952cbd5e404c7881c5f6b92fa25b21ec` (checked by the script); the source Git blob for `common.py` is `10d8321d8c0f831d474c3769adae7cc7caf29821`. Its only imported function is the published reconstruction of the Leech shell; [the source release](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2) separately verifies the baseline.

Run `bash replay.sh` with NumPy and SciPy to rebuild the candidate bodies and both radical-boundary checks in the prescribed order. The two final geometric checks never treat arbitrary pre-existing candidate JSON as a complete enumeration. After reducing the old comparisons to a finite set of circle halfspaces, any nonempty closed feasible proper subset has a boundary point where some halfspace binds. The entire-circle case is covered by an explicitly tested old root direction. Boundary candidates have coordinates in a quadratic extension; every inequality uses integer sign and square comparisons. The displayed exact scalar body products have small integer partial sums even when computed through binary32 (by Cauchy–Schwarz).

This theorem says nothing about moving or removing the frozen caps, making some old midpoint fibres one-sided, or changing existing bodies or radii. Those are genuine next directions rather than additional searches of the same chart.
