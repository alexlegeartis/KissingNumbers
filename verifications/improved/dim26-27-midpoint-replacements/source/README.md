# Shared midpoint replacements in dimensions 26 and 27

The current-bound column follows the [29 October joint-refill continuation](../../dim27-31-joint-refills/).

Retained construction: the dated benchmark, construction and gain columns below are historical; the current-bound column uses the [26D asymmetric exchange](../../dim26-30-asymmetric-exchanges/) and [27D quadratic module](../../dim27-quadratic-normal-module/).

Exact lower bounds obtained on 2026-10-25:

| Dimension | Previous benchmark | This construction | Improvement | Current repository bound (29 October 2026) |
|---|---:|---:|---:|---:|
| 26 | 199,868 | **200,918** | **1,050** | 201,136 |
| 27 | 201,589 | **204,636** | **3,047** | 209,807 |

The comparison benchmarks are those published in the [Takhanov–Yun certificate repository](https://github.com/k-nic/kissing_d25-31). The construction here rebuilds from the standard Leech shell and does not use their lifted configurations.

## Verify

Only Python and NumPy are needed. From this directory:

```text
python verify.py certificate26_200918.npz
python verify.py certificate27_204636.npz
```

Both commands have passed. The verifier rebuilds the Leech shell without caches or downloads, identifies deleted points by their coordinates, verifies the midpoint identities, and checks every new/new pair using integer arithmetic. The bulk/new and auxiliary pairs are covered by the elementary arguments below. No optimizer, normal-frame compatibility table, or fibre cache is used in this verifier.

| Component | 26D | 27D |
|---|---:|---:|
| Original Leech points | 196,560 | 196,560 |
| Deleted Leech points | −4,352 | −8,064 |
| Inserted points | +8,704 | +16,128 |
| Pure tail points | +6 | +12 |
| Total | **200,918** | **204,636** |

## Geometry and proof

Use squared norm 32 throughout. The kissing condition is inner product at most 16.

Let `u,v` be Leech minimal vectors with `u·v=16`. Their midpoint `h=(u+v)/2` has squared norm 24. A tail `t` of squared norm 8 gives a point `P=(h,t)` of squared norm 32.

For every original Leech point `z` distinct from both endpoints,

```text
h·z = (u·z + v·z)/2 <= 16.
```

Thus deleting the two endpoints suffices to clear the entire remaining equator. Several midpoint heads can share those deletions. The endpoint indices stored in the certificate are checked against the explicit deleted-coordinate array.

Here is the local replacement used in the search. Take mutually orthogonal `r,p,q`, each of squared norm 16, with all four points `r±p,r±q` in the Leech shell. Delete those four points. Write `a=(2,2)` and `b=(2,-2)` in a tail plane, so `a·b=0` and `|a|²=|b|²=8`. Insert

```text
(r ± (p+q)/2, ±a)  [all four sign combinations]
(r ± (p-q)/2, ±b)  [all four sign combinations].
```

Each foot is a midpoint of two deleted endpoints. These eight points form a translated cross-polytope in the local four-dimensional space. Their maximum mutual inner product is 16. The replacement gains four points.

The two tail directions correspond to two matchings on the same four deleted vertices. This is the source of the gain: both matchings use the same deletion set. Midpoint heads within one matching cannot share an endpoint while meeting the required same-direction separation, so this sharing uses different tail directions.

An octad split of the Leech shell supplies 4,320 suitable residual fibres. The eight normal axes of each fibre form one of 135 frames; every frame occurs with 32 residuals. A compatible selection of one normal-axis pair per frame activates a group of 32 replacements, gaining 128 points. The 26D certificate uses 34 groups. In 27D the tail plane may be `xy`, `xz`, or `yz`; the certificate uses 63 groups.

Every inserted coordinate is an integer. All inserted tails have squared norm 8. Therefore any pure tail point of squared norm 32 clears every inserted point by Cauchy–Schwarz. Add a regular hexagon in 26D, or the twelve cuboctahedron vectors `(±4,±4,0)` and permutations in 27D. These pure points also satisfy the kissing condition among themselves. They are distinct from the inserted points because every inserted foot has squared norm 24.

All new/new inequalities, including inequalities between different groups and tail planes, are checked directly. The largest is exactly 16. Bulk/bulk inequalities are those of the standard Leech minimal shell generated from the extended binary Golay code.

## Certificate format

Each NPZ contains:

- `new_points`: explicit integer coordinates, excluding pure tail points;
- `removed_leech`: explicit coordinates of deleted original Leech points;
- `midpoint_owner_indices`: two indices into `removed_leech` for each inserted point.

No external shell ordering is assumed.

SHA-256:

```text
6a1fcba8c37657cd6111bb2954c354e5a2fc3668ade9c6551c5c13db3156ca43  certificate26_200918.npz
592ba5d22deebeb43f15ab3600d8ab75335c3b0591ab5ebbd7a098a0fa1e54ad  certificate27_204636.npz
```

## Search provenance and limits

The search scripts are in the parent directory: `two_axis_codelets.py`, `augment_flags26.py`, and `multi_plane27.py`. Discovery used short greedy probes and small exchanges in the normal-frame graph. Verification does not depend on the search.

These are achieved lower bounds, not optimality results for the configurations or their construction family. No improvement in dimension 25 is claimed here.
