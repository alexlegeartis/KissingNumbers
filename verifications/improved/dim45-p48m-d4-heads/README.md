# A 7,655,213-point kissing configuration in 45 dimensions

This raises the preceding repository bound of [7,655,091](https://github.com/alexlegeartis/KissingNumbers/tree/91c1bef4ebe42f0b774adb0f29b069af2697d481/verifications/improved/dim45-p48m-d4-heads/) by 122. The construction uses 53 complete affine families from a 44-dimensional projection of the extremal lattice P48m, independently chosen rational scalar coordinates, and both pure poles on the extra axis. The published P48m input is credited to [the Dorofeev–Sun–Wang certificate](https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/packings_out/records/dim45-kissing7379838/P48m-d914e7d8.cert.txt); this is a different construction and does not claim their result. It is a lower bound, not an optimality proof or a claim of exhaustive priority.

## Construction and proof

Use the four integral normal rows in [section.json](section.json), in the same P48m integral basis as the [pinned 44D exact certificate](https://github.com/alexlegeartis/KissingNumbers/tree/53069eac091a87601d1d570c8b5fbab933b12499/verifications/improved/dim44-d4-section). The new normal Gram is the same scaled D4. The included complete 1,656-point slice (compressed signed-byte rows) yields calibration counts (0, 12, 8, 130). Its completeness follows from the same 11-design moment argument for any two normal roots of these norms and mutual product. The affine count formulas of the pinned parent apply to every such D4 embedding. Together they give 7,655,211 points in 53 complete families, plus two poles. Choosing the normal frame from a larger shell rather than holding its first two rows fixed improved the populations; two radius-36 families contribute 12 and 8 points. No sharper lower bound on inter-family distance is needed.

A row of the witness contains a dual class, an integer R equal to 12 times the projected squared norm, and an integer H specifying the extra coordinate H/(8 sqrt(12)). Normalize the resulting vector. For noncoincident bodies in classes c and c', the projected squared distance is at least six minus the minimum squared norm in the normal difference class: the difference can be lifted to a nonzero parent vector using a shortest normal representative. If the two rows have equal class and radius, coincident bodies at different heights must also be allowed, and the checker instead uses zero as the distance lower bound. Different rows in the witness specify different class/radius/height triples. There are 324 normal dual classes and their shortest representatives have squared norm at most three. As proved in the parent package, this also ensures primitivity and disjointness of different projected classes.

Let D be 12 times the appropriate distance lower bound. The squared norm of a raw lifted vector is (64R+H²)/768. For a pair from different rows, the checker verifies that the integer 64(R+R'-D)+2HH' is no greater than the square root of (64R+H²)(64R'+H'²), checking its sign before squaring. This proves their normalized dot is at most 1/2. The within-family test uses minimum squared body distance six, and both poles are valid precisely when 4H² is no larger than 64R+H². The replay covers all 1,378 cross-family inequalities, 53 within-family bounds, and both poles.

## Replay

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim45-p48m-d4-heads/verify.py --moments
```

The run verifies the affine count identities and geometric inequalities exactly. It inherits the unchanged minimum-six lattice certificate and completeness of the supplied coordinate slice from the pinned parent package; run its `verify.py` to recheck those premises. The included [exact replay receipt](verification.recorded.json) gives all 53 populations, independently rechecked from the full integer slice. No fresh enumeration of millions of points is needed.
