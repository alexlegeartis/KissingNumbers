# 36D Coxeter–Todd triangle switches and phase refill

This certificate gives **502,902** kissing points, 81 above the
[502,821-point construction](../dim36-ct-nonlocal/). It uses three
orthogonal copies of the Coxeter–Todd lattice. The vectors have squared
norm 12 and distinct vectors have inner product at most 6.

The two-block layer contains cells with block minima `(0,4,8)` up to
permutation. Each of the 81 norm-eight vectors in a cell is paired with
the triangle of norm-four vectors in the small block. A switch moves
the negative of that triangle into the previously empty block and
keeps the norm-eight vector in place. The certificate switches 63
triangles across seven cells. The two-block layer retains its cardinality.

The three-block layer uses a single norm-four root in each block.
The search jointly selected triangle switches and individual root
triples in a Boolean conflict model. A candidate forces the switches
it needs and pays for every old root triple it excludes. A switch
also excludes any retained or added triple it would hit. The final
delta inserts 1,390 root triples and removes 1,309.

`certificate.json` pins the parent data by SHA-256. A root triple lists
its three quotient labels followed by its three phases. Phases are
obtained by iterating the order-three isometry given by the identity
minus the parent's similarity matrix. A switch lists the parent cell
index and the norm-eight vector in the parent's integral basis.

Both replays use Python and NumPy.

```sh
python verify.py
OPENBLAS_NUM_THREADS=8 python verify_expanded.py
```

`verify.py` enumerates the short shells, compares the possible two-block
placements in the quotient, checks their interfaces with the axial
and three-block layers, and checks all three-block pairs that share a
root. If no root is shared, each block contributes at most 2.

`verify_expanded.py` reconstructs all 502,902 distinct points and
checks every pair involving any of the 1,579 changed points. The
792,834,848 products are computed in binary32 with an absolute
integer partial-sum bound of 864, so every operation is exact.
Pairs of retained points inherit the parent certificate.

`python final.py` recomputes the point count from the pinned data.
