# A 7,655,091-point kissing configuration in 45 dimensions

This raises the preceding repository bound of 7,655,056 by 35. The construction uses 51 complete affine families from a 44-dimensional projection of the extremal lattice P48m, independently chosen rational scalar coordinates, and both pure poles on the extra axis. The published P48m input is credited to [the Dorofeev–Sun–Wang certificate](https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/packings_out/records/dim45-kissing7379838/P48m-d914e7d8.cert.txt); this is a different construction and does not claim their result. It is a lower bound, not an optimality proof or a claim of exhaustive priority.

## Construction and proof

Use the first two normals of the [pinned 44D exact certificate](https://github.com/alexlegeartis/KissingNumbers/tree/53069eac091a87601d1d570c8b5fbab933b12499/verifications/improved/dim44-d4-section) and the negatives of its section points with zero-based indices 329 and 588. Their Gram matrix is again that of the same scaled D4 normal lattice. Recounting the four calibrated profiles in the certified complete 1,656-point slice now gives (0, 10, 4, 144). The affine count formulas in the 44D verifier were proved for every embedding with this Gram matrix using spherical moments and harmonic theta series, so they apply with these new counts. The resulting 51 populations sum to 7,655,089.

A row of the witness contains a dual class, an integer R equal to 12 times the projected squared norm, and an integer H specifying the extra coordinate H/(8 sqrt(12)). Normalize the resulting vector. For noncoincident bodies in classes c and c', the projected squared distance is at least six minus the minimum squared norm in the normal difference class: the difference can be lifted to a nonzero parent vector using a shortest normal representative. If the two rows have equal class and radius, coincident bodies at different heights must also be allowed, and the checker instead uses zero as the distance lower bound. Different rows in the witness specify different class/radius/height triples. There are 324 normal dual classes and their shortest representatives have squared norm at most three. As proved in the parent package, this also ensures primitivity and disjointness of different projected classes.

Let D be 12 times the appropriate distance lower bound. The squared norm of a raw lifted vector is (64R+H²)/768. For a pair from different rows, the checker verifies that the integer 64(R+R'-D)+2HH' is no greater than the square root of (64R+H²)(64R'+H'²), checking its sign before squaring. This proves their normalized dot is at most 1/2. The within-family test uses minimum squared body distance six, and both poles are valid precisely when 4H² is no larger than 64R+H². The replay covers all 1,275 cross-family inequalities, 51 within-family bounds, and both poles.

## Replay

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim45-p48m-d4-heads/verify.py --moments
```

The run verifies the affine count identities and geometric inequalities exactly. It inherits the unchanged minimum-six lattice certificate and completeness of the supplied coordinate slice from the pinned parent package; run its `verify.py` to recheck those premises. The included [independent integer replay receipt](verification.recorded.json) gives all 51 populations. No fresh enumeration of millions of points is needed.
