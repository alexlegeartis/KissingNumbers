# Limits of tail recombination for the 232,230-point 29D witness

This is a scoped obstruction, **not** a new kissing-number estimate. The baseline is [Takhanov's v0.2 232,230-point configuration](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2), not our separately certified 232,157-point configuration. All input below is from his compact certificate, with exact integer squared norm 128 and separation threshold 64.

## Mixing full tail-isometric copies

Write a ray as `(b,t)` with a 24-coordinate body and five-coordinate tail. A tail isometry `(b,t) ↦ (b,Rt)` preserves the code and fixes the equator. The two copies share the equator and any common off-equator rows. Every common row is compatible with every row in either copy.

After removing the common rows from both sides, build a bipartite graph on the two equal-sized remainders, putting an edge on a cross inner product above 64. The two sides are individually compatible. By the bipartite matching/vertex-cover theorem, the largest code in this **union** is the old size plus the deficiency of a maximum matching on either side. In particular, a perfect matching is a short certificate that **no paid joint exchange**, even using all donors simultaneously, gains a point. The global rotation itself is a tie, not a new record.

| Tail operation (indices 0–4) | Rows per asymmetric side | Cross conflicts | Matching deficiency |
| --- | ---: | ---: | ---: |
| Negate all five axes | 866 | 5,290 | 0 |
| Negate axis 0 | 25,820 | 801,482 | 0 |
| Exchange axes 1 and 2 | 32,654 | 961,630 | 0 |
| Quarter-turn in axes 0 and 4 | 46,276 | 1,716,138 | 0 |

The row norms and small coordinates imply that every float32 dot sum in the implementation is exact: all products are integers and every partial sum has magnitude at most 128, by Cauchy–Schwarz. More directly, the replay checks every matched edge again in signed integer arithmetic and checks that the matching is a bijection. We do not infer that other tail isometries or rotations varying with the body fail.

## A real tail cannot freely complete any of the 193 singleton midpoints

There are 64,905 midpoint lifts over 32,549 distinct bodies of norm squared 96. Each fibre holds at most two points, and if it holds two their tails are antipodal: on a common body separation demands `96 + t·u ≤ 64` while `t² = u² = 32`. Thus precisely 193 bodies have just one lift.

Fix such a body and replace its singleton by a pair `(b,±u)`. To keep every other row `(b_j,t_j)`, one needs

```
|u·t_j| ≤ 64 − b·b_j, with u·u = 32.
```

Equatorial constraints are unchanged because the body was already present. Seventeen of the 193 candidates have an unavoidable negative right side. For 42 more, the zero-right-side constraints span the entire tail space. The other 134 live in subspaces of dimension one, two or three. Their slab polytope is bounded: the retained pure-tail signed two-coordinate roots bound every tail coordinate. In a basis of the zero-right-side nullspace, enumerate vertices by intersecting independent bounding planes and checking all remaining exact inequalities. The maximum squared norm, attained at a vertex by convexity, is at most **20**, strictly below the required 32. This covers **arbitrary real** rotated antipodal pairs replacing a single midpoint with no other deletions, not merely lattice tails.

The obstruction suggests changing midpoint bodies or co-designing several interacting fibres and another layer with paid deletions. It does not exclude either.

## Replay

Download the pinned [`d29_certificate.npz`](https://github.com/k-nic/KissingNumbers/blob/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/data/d29_certificate.npz) beside the scripts. Its SHA-256 is `b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824`, matching [Takhanov's certificate manifest](https://github.com/k-nic/KissingNumbers/blob/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/SHA256SUMS.txt).

```sh
OPENBLAS_NUM_THREADS=1 python mirror_search.py selected  # NumPy, SciPy
OPENBLAS_NUM_THREADS=1 python rotate_pair.py             # NumPy, SymPy
```

The second replay outputs one record per singleton and its actual upper radius after all tail slabs. Results rely on the pinned, separately verified baseline. The scripts do not recompute its full all-pairs certificate.
