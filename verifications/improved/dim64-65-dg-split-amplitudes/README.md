# Split-amplitude Delsarte–Goethals replacements in 64 and 65 dimensions

| Dimension | Previous exact lower bound | New exact lower bound | Increase |
|---:|---:|---:|---:|
| 64 | 418,478,976 | **422,673,152** | 4,194,176 |
| 65 | 419,000,274 | **422,678,354** | 3,678,080 |

In both dimensions a selected dense sign vector is replaced by two unequal-amplitude vectors. The support stock and populations of the sparse layer are unchanged, but its signs are rephased. Remove the 128 ordinary coordinate poles in both dimensions and all 516,096 affine-hyperplane head points in 65D. Retain the [optimal 65D four-site packing](../dim65-optimal-quadruple-layer/) and the two tail-axis poles. As for that predecessor, the 65D count uses the established existence theorem for its four-site packing; the rest is a finite explicit construction.

## Two donors per selected owner

The parent dense body is the 80 listed cubic cosets of RM(2,6), of size 335,544,320 and binary distance at least 16. Enumerate coordinates in binary truth-table order and pair adjacent coordinates into 32 disjoint pairs. Use the first 32 distinct cubic representatives, one per pair. In each assigned cubic coset designate 131,072 owners: its cubic translate of the full Delsarte–Goethals (DG) family of 1,024 quadratic parts and 128 affine shifts. The included exact generator computes the quadratic parts over GF(32) and checks that the minimum distance within each such DG family is 24, by testing all 15,903 distinct nonzero quadratic differences against all 128 affine parts. Its monomial convention is the same as the parent's.

For every selected owner, replace its two ordinary magnitudes in the assigned pair by `(a,b)` and `(b,a)`, respectively, where `a=28/5` and `b=29/20`. Flip the owner's sign only at the large entry; the other 62 magnitudes remain one, and the tail coordinate (if present) is zero. Their common squared norm is `N=62+a²+b²=38185/400`.

The complete exact dot caps between donors are:

* Same owner: `62−2ab`.
* Distinct owners of the same pair, same profile: `16+a²+b²−2`; swapped profiles: `16+2(ab+1)`. Here the unweighted owner dot is at most 16 by the DG minimum 24.
* Distinct pairs: `32+2(a+b)`. Here the unweighted owner dot is at most 32 by the parent minimum 16, and the sites being modified are disjoint.

With a retained dense point the dot is at most `32+a+b`, since its owner alone, not any of its neighbors, was retired. Any signed four-site point has dot at most `2+a+b`, including in the 65D design. Donors are perpendicular to the two retained tail poles. Every angular inequality is checked by exact squared cross multiplication against the original squared norms.

## Rephase the sparse sign stock

Use the exact 40,414 zero-tail weight-16 masks from [the 64D predecessor](../dim63-68-affine-hyperplane-heads/): 37,996 original cubic-shell supports, 2,354 signed-overlap additions and 64 zero-tail joint additions. All were already present in the 65D predecessor as well. The verifier independently reconstructs every parent word type, the sorted weight-16 shell, and every dangerous sparse pair. Reassign the marked parity guards from the short hexadecimal witness while requiring total parity even on every mask.

Every intersection of size nine or ten has a shared eight-site check of opposite sign parity. Every bank's imposed signed checks are consistent, their homogeneous span contains the full mask, and all other nonzero rows have weight eight. For a basis modulo the full mask, append balanced binary columns by repeatedly splitting the equally sized current column classes. Equal sizes follow by binary Fourier inversion from the weight-eight checks. The completed dual has five rows; all its 16 column labels are distinct and start with the overall parity bit. Thus its affine primal coset still has 2,048 words and distance at least four. The shared guards give a sign mismatch in every size-nine or size-ten intersection, so no sparse population is lost.

Every mask is orthogonal modulo two to all 22 RM(2,6) generators and all 80 cubic representatives. If a mask contains a donor's large spike, the donor restriction has odd sign parity and therefore disagrees at least once with any old sparse sign word. With both spikes or the large spike only the dot bounds are `12+a+b` and `13+a`, respectively. With the small spike only or neither, even without a sign disagreement the bounds are `15+b` and 16. These satisfy the exact angular bound against squared norm 16.

## Counts and replay

In both dimensions retire 4,194,304 owners and insert 8,388,608 donors, a net 4,194,304 point increase before the pole and atlas deletions. The 64D predecessor includes 10,416 affine 2-flat supports with all sixteen signs. In 65D retain 10,741 four-site supports by the [Bao–Ji quadruple packing theorem used in the pinned predecessor](../dim65-optimal-quadruple-layer/); each still accepts all signs because the new compatibility bound uses only its size.

Run `python3 -B verify.py` in this proof directory with a C++20 compiler. It SHA-pins the published input witnesses and claim rows and independently replays the base palette, all 816,625,491 support-pair comparisons, the 14,723 guards, the DG count and distance, all exact mixed gates, and both counts. No floating-point decisions, search logs or network access are needed.
