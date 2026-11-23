# Additional norm-four tail layers in 84D and 85D

| Dimension | Signed-mask predecessor | Added points | Final lower bound |
|---:|---:|---:|---:|
| 84 | 6,292,536,896 | 107,136 | **6,292,644,032** |
| 85 | 6,297,207,062 | 114,048 | **6,297,321,110** |

The underlying lattice is Γ72. The signed-mask predecessor for 84D is
included; its first two steps improve the earlier public 6,292,388,112
using the same 840-point tail and smaller individual horizontal packets.
For 85D the script replays the [signed norm-4+8 predecessor](../dim85-gamma72-signed-apron/).
The final extension does not remove any signed-mask points.

The new tails have norm squared four. They use 168 (D_{12}) roots
after deleting those on 24 pairs, 384 even-signed octads, and 192
additional vectors with one double and four unit coordinates, all
in coordinates rescaled by (1/sqrt2). The 85D extension adds
48 more vectors involving the last coordinate. The vectors are partitioned
by distinct nonzero B-only horizontal phase labels. A tail at each such
phase is paired with 144 exact horizontal positions of squared norm four:
the 48 signed raw Leech norm-eight vectors in any of the three blocks.

Two new tails in the same phase have nonpositive dot product; between
phases the tail dot product is at most two, and the horizontal
phase clears a further four in squared distance. Same-tail horizontal
vectors are pairwise orthogonal or antipodal. Old mixed tail dot products
are at most two, and their different first color clears at least six
horizontally. Old pure tails are twice the old mixed tails, and Γ72
equatorial points clear by the phase norm and new tail norm. The
finite inner products and counting are checked without a search solver.

Run `OPENBLAS_NUM_THREADS=1 python3 -B verify.py` with NumPy and SymPy;
pass the shared parent of the predecessor verification folders to run
outside the repository. signed84.py is the included first-stage replay
and norm4.py checks the new tail configurations in both dimensions.
