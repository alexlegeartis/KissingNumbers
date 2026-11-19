# 36D: 502,821 by opening several CT fibres jointly

This is an exact certificate of **502,821** points, up **54** from the
[earlier 502,767](../dim36-37-ct-cell-rebuild/). That earlier package
also verifies its associated 37D construction; it is unchanged.

The 403,704-point axial and two-block part stays fixed. In the ternary
quotient of three 12-dimensional Coxeter–Todd lattices, the new layer
selects 3,671 cells with block minima (4,4,4), instead of 3,669. Each
cell has 27 points. This selection was found by exact clique searches
in larger neighborhoods: an outside challenger opens its incident
fibres in all three coordinates, so a change can reach several old
fibres together. The certificate does not depend on the search's
completeness or numerical assumptions.

Both verifiers independently enumerate the short lattice shells and
check the 5,100 outer cosets. Their minimum pairwise squared quotient
distance is 12. The kernel lattice also has squared minimum 12, by the
exhibited integral similarity. Thus distinct points within a coset
or between cosets obey the kissing-angle condition at squared norm
12. This is a lower bound, not a bound on the optimum of the model.

Run from this directory:

```sh
python verify.py                 # standard library
python verify_independent.py     # NumPy
```
