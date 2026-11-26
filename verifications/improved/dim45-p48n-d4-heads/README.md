# 7,655,249 kissing points in dimension 45 from P48n

The [preceding P48m construction](https://github.com/alexlegeartis/KissingNumbers/tree/6d5f57eb0454618b28c109c486417e266ffda1aa/verifications/improved/dim45-p48m-d4-heads) gives 7,655,213 points. Changing the ambient extremal lattice **from P48m to P48n**, not rerunning that P48m frame, gives **7,655,249** (+36). The 53 projected classes, independent scalar heads, affine count identities and angle criterion are inherited from that construction. There are still exactly two pure poles. This is a lower bound, not an optimality certificate.

## What changes

Use the Gram matrix `G` for P48n in [the pinned parent embedding](https://github.com/alexlegeartis/KissingNumbers/blob/6d5f57eb0454618b28c109c486417e266ffda1aa/verifications/improved/dim41-46-p48-heads/partial43/parent_embedding.npz), and the four integer normal rows in `section.json`. They span a copy of √3 D4. P48n is the known extremal even unimodular lattice; its minimum and the modular-design identities for its shells are the same mathematical inputs as in the parent packages. The new calibrations in the 1,656-point slice are (0,14,8,136) instead of (0,12,8,130). Applying the **same** [D4 affine moment identities](https://github.com/alexlegeartis/KissingNumbers/tree/53069eac091a87601d1d570c8b5fbab933b12499/verifications/improved/dim44-d4-section) gives 7,655,247 body-and-head points plus the two poles. No assumption that the P48m and P48n complements or slices are isometric is made.

## Independent completeness and distance checks

Write a for the first normal and b for the negative of the second. Both have squared norm six and a·b=3. The stored 1,656 distinct integer vectors x all have squared norm six, a·x=3 and b·x=0. For completeness put u=a·x, t=b·x, P(u)=(u+3)(u+2)(u+1)u(u−1)(u−2), and Q(t)=−(t−1)(t−2)(t−3). The universal 11-design moments give ΣP(u)=26,611,200 and ΣP(u)Q(t)=9,616,320. Unless x=±a the first polynomial selects u=3. On u=3, t lies between zero and three except for x=b and x=a−b: apply the lattice minimum to x and to a−x. Subtract the four exceptional contributions: (u,t)=(6,3),(−6,−3),(3,6),(3,−3). Dividing the resulting second sum by P(3)Q(0)=4,320 proves that exactly 1,656 such vectors exist. All four other normal profiles can therefore be counted from this complete, supplied slice.

For each witnessed family let R be twelve times the squared projected radius and H the integer scalar-head numerator. A distinct body difference from two classes lifts to a nonzero P48n lattice vector using the shortest normal representative in its difference class. Thus its squared distance in body units is at least six minus that representative's squared norm. In a repeated class and radius, distinct lifts might have coincident body, and the cross-family check uses zero instead. The verifier handles all these cases, the within-family distances, all 1,378 cross-family pairs, and both poles without a floating-point angle check. The exact populations of the norm-six and norm-eight layers are universal affine functions of the four calibrated counts, whose identities are replayed with `--moments`.

From repository root (assuming this package at `verifications/improved/dim45-p48n-d4-heads/`):

```sh
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim45-p48n-d4-heads/verify.py --moments
```

The stated improvement compares the two explicit repository witnesses. It is not a claim to have exhausted all ambient lattices or all 45-dimensional configurations.
