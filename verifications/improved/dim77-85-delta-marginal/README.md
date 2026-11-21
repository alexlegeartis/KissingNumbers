# A third exact shell marginal in dimensions 77–85

The same [shaped norm-ten lattices and tails](../dim76-85-mixed-shape-norm10/) give slightly stronger existence bounds by reducing the upper bound on the number of triangles through a graph vertex. These are not explicit coordinate lists or sphere-packing density bounds.

| Dimension | Previous | Lower bound | Gain |
|---:|---:|---:|---:|
| 77 | 6,218,669,616 | 6,218,669,870 | 254 |
| 78 | 6,218,991,624 | 6,218,992,206 | 582 |
| 79 | 6,219,415,042 | 6,219,416,084 | 1,042 |
| 80 | 6,219,960,348 | 6,219,962,106 | 1,758 |
| 81 | 6,220,773,344 | 6,220,774,236 | 892 |
| 82 | 6,222,087,772 | 6,222,089,254 | 1,482 |
| 83 | 6,223,765,754 | 6,223,768,006 | 2,252 |
| 84 | 6,225,860,206 | 6,225,863,462 | 3,256 |
| 85 | 6,228,414,172 | 6,228,418,680 | 4,508 |

For norm-ten (v,w,u) from Nebe's `Gamma72`, put (a=u\cdot v), (b=u\cdot w) and (j=v\cdot w). In addition to the joint moments and the separate exact marginals of (a,b), use the exact population of (a-b=u\cdot(v-w)) whenever (j=5) or (6). At (j=5) the axis (v-w) itself has squared norm ten, so its population distribution is already known; exclude the exceptional `±(v−w)` heads from the residual LP, and account for them explicitly. At (j=6) the axis is a norm-eight vector, and the marginal about it follows from its fixed-angle moments or, equivalently, double counting the norm-eight population about a norm-ten vector. The case (j=4) retains the preceding two-marginal proof.

The verifier rebuilds the exact shell arithmetic, checks every rational dual inequality, and uses the validated blocked Boßecker–Rautenbach recurrence to turn the local triangle bound into an independent-set lower bound. It relies on [Nebe's lattice](https://arxiv.org/abs/1008.2862) and the inherited [shell-design theorem](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf), not on numerical solver output. The exact replay only needs NumPy:

```sh
python verifications/improved/dim77-85-delta-marginal/verify.py
```
