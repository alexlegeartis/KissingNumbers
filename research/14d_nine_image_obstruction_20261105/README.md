# Exact obstruction for nine 14D group images

The supplied `baseline.bin` is the 1,932-point 14D configuration reproduced in the repository's group benchmark. Its SHA-256 is `25fd242c2c5b52afb8f065fc4651d466d19f7b4016809eb847e17f1df5099e6b`. It contains signed int8 numerator coordinates in row-major order (1,932 rows by 14 columns); each row has squared norm 8. This is a known lower bound, not a new record.

The two generators in `discovery/recipes/group14_strict504.json` close to a rational orthogonal group of order 504. Acting on the baseline yields exactly nine distinct configurations. They share 1,764 points and their union has 3,276 points. Each of the remaining 1,512 points is compatible with the entire common core.

The graph on those 1,512 points has an edge exactly when their inner product exceeds one half. It is 32-regular. The verifier computes, with exact integer arithmetic,

`(A-32I)(A-24I)(A-16I)(A-12I)(A-6I)(A-4I)(A-2I)A(A+4I) = 0`.

Since A is symmetric, its eigenvalues belong to these nine real roots and are at least -4. Hoffman's independent-set bound therefore gives 1,512*4/(32+4) = 168. A known image attains 168 outside the common core, so the largest kissing configuration inside the entire nine-image union has exactly 1,932 points. This is a restricted-family result: it does not rule out additions outside the union or removing core points and replacing them by other points.

Run from the repository root, with NumPy and SciPy installed:

`OPENBLAS_NUM_THREADS=1 python research/14d_nine_image_obstruction_20261105/verify.py`

The script checks the source checksum, the baseline and group generators, closure of the group, the nine images, all cardinalities, the exact conflict graph, and the matrix polynomial. It bounds intermediate integer values before each sparse matrix product.

## Extending the obstruction by an outside orbit

The verifier additionally constructs an orbit of 56 points outside the nine-image union, with an order-nine stabilizer. These points are mutually compatible and compatible with the common core. Each conflicts with 27 vertices in one of the two 756-vertex components of the variable graph; these 56 neighborhoods cover that component twice. Every neighborhood induces a 27-vertex 16-regular Schläfli graph, whose independence number is at most three by its exact quadratic adjacency identity.

If a compatible set contains `m` old vertices from that component and `k` of the new orbit points, summing the 56 neighborhood bounds yields `2m + 3k <= 168`. Thus these choices contribute at most 84 points in total. The other component contributes at most 84, and the common core at most 1,764. The expanded finite union also cannot improve 1,932. This remains a statement about the specified finite orbit, not other subgroup orbits or deformations.
