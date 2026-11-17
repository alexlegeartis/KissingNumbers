# Gamma72 norm-ten layers with exact joint moments, marginals and triangle-sensitive independence bounds

Eight verified existence bounds. Gains over the preceding repository constructions are small: 144,804–3,345,876 points (under 0.054%). No exhaustive worldwide-priority or optimality claim is made. The selected independent sets are proved to exist, not supplied as full coordinate lists. These are general kissing-number lower bounds, not sphere-packing-density or lattice-kissing claims.

| Dimension | New lower bound | Previous repository bound | Added points |
|---:|---:|---:|---:|
| 79 | 6,219,067,454 | 6,218,922,650 | 144,804 |
| 80 | 6,219,714,034 | 6,219,474,300 | 239,734 |
| 81 | 6,220,589,740 | 6,219,649,576 | 940,164 |
| 82 | 6,221,706,936 | 6,220,170,666 | 1,536,270 |
| 83 | 6,223,067,174 | 6,221,008,658 | 2,058,516 |
| 84 | 6,224,663,766 | 6,222,025,404 | 2,638,362 |
| 85 | 6,226,484,034 | 6,223,222,688 | 3,261,346 |
| 86 | 6,228,511,388 | 6,225,165,512 | 3,345,876 |

## Method and premises

The source is [Nebe’s even unimodular Gamma72 with squared minimum8](https://arxiv.org/abs/1008.2862). All its nonempty shells have design strength11 by [Venkov’s theorem, presented by Nebe](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf). These published existence, minimum and design theorems remain explicit premises; this package does not reconstruct the lattice. `inherited_source.py` replays the source theta coefficients and exact one-dimensional product distributions.

For tail rank k=n−72, use all signed weight-four tails. Joint spherical moments and exact one-dimensional marginal populations give rational LP dual upper bounds for common neighbors in the incompatibility graph. A uniform integer triangle count feeds the [Boßecker–Rautenbach theorem](https://doi.org/10.1016/j.disc.2010.05.026), evaluated by forty conservative rational blocks. This theorem permits triangles; the triangle-free-only Shearer theorem is not being misapplied.

The stored dual coefficients are in the UNscaled row basis:36 normalized monomial moments (even total degree≤10), followed by26 indicator rows (first coordinate −6…6, then second coordinate −6…6). Exact verification checks every profile. Neither floating LP feasibility nor optimality is trusted.

See [the full derivation](PROOF.md) and [independent review](INDEPENDENT_AUDIT.md). The detailed proof begins with82D; replace k,M,D,T and block endpoints by the dimension-specific frozen files. All eight final marginal variants passed separate complete independent replays.

## Replay

With Python, NumPy and SymPy installed:

```sh
python verifications/improved/dim79-86-norm10-triangle-layers/verify.py
```

The checker verifies hashes, replays source arithmetic, and independently regenerates each tail orbit, closed-form mixed moment, marginal population, rational dual inequality, triangle bound, recurrence and geometric compatibility. Frozen candidates, exact parameters and replay logs are included.

External comparisons in the main table remain the attributed Gamma72 plus tail direct-sum values. The previous repository values above are separately identified internal history. The manuscript remains its earlier explicitly deferred snapshot.

Verifier entry points: `verify.py` runs the full chain; `verify_one.py` independently checks each dimension; `check_recurrence.py` supplies a second symbolic/rational recurrence review; `inherited_source.py` provides the replayed source arithmetic; `final.py` exports the checked `claims.json` to the repository registry.
