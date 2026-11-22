# 40D affine centre and minimum-shell refill

This construction gives **1,126,327** kissing points, five above the
[previous affine-shell certificate](../dim40-affine-shells/). It uses centre
267 in the ternary discriminant group of the rank-40 P48n complement.
The 1,084,590-point core is accompanied by 40,590 complete-family points of
coefficient squared norm 36 and 1,134 of squared norm 30. The refill adds
ten vectors of squared norm 6 and three of squared norm 36.

The search enumerated all 2,240 nonzero isotropic centres and the entire
1,092,000-point rank-40 minimum shell. At each centre it tested the minimum
shell against the affine core and the short families about that centre.
Forty-eight centres admitted ten vectors. The previous centre admitted
eight.

For the core test, write `Q` for the rank-40 Gram matrix, `x` for a
core row of squared norm 54, and `y` for a refill of squared norm 6.
A forbidden product `x Q yᵀ` is at least 12 and divisible by
three. The row `x - 3y` then has squared norm at most 36 in the same
centre class. The gluing minimum gives squared norm at least 36, so every
potential blocker is detected among the 495 norm-36 bodies in that class.

The three further points belong to class 3006. Their only possible
conflicts with the affine core would give a norm-30 difference in class
6252. The replay reconstructs this short difference shell and excludes
every such contact.

```sh
python verify.py
OPENBLAS_NUM_THREADS=16 python verify_expanded.py
```

The first replay checks the P48n complement and gluing map, the complete
short families, and all refill interactions using integer arithmetic.
It uses the P48n minimum and the shell population identities from the
[body-family certificate](../dim41-46-p48-heads/).

The expanded replay compiles `enumerate_core.cpp` with g++ and uses
C++17 and `__int128`. It checks positive definiteness and the rank-40
lattice minimum, enumerates the entire new affine core, and replays
every cross-layer product. Its matrix products have exact integer
accumulations in binary32.

`certificate.json` contains the centre index and thirteen refill coefficient
rows. `auxiliary.npz` contains both explicit short layers.
`python final.py` recomputes the ledger total.
