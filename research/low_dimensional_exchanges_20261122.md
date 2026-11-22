# Discrete exchanges in dimensions 25 through 31

On 22 November 2026 the search raised the 30D lower bound from 248,416 to **248,506**. The other six dimensions retained their starting benchmarks. The [30D construction](../verifications/improved/dim30-row-midpoint-packing/) contains the integer coordinates, independent geometric replay, candidate generator, and finite packing instance.

| Dimension | Starting benchmark | Best verified count | Gain |
| ---: | ---: | ---: | ---: |
| 25 | 197,586 | 197,586 | 0 |
| 26 | 201,204 | 201,204 | 0 |
| 27 | 210,108 | 210,108 | 0 |
| 28 | 223,528 | 223,528 | 0 |
| 29 | 232,157 | 232,157 | 0 |
| 30 | 248,416 | **248,506** | **90** |
| 31 | 251,524 | 251,524 | 0 |

In dimension 27 the benchmark is [Takhanov’s certificate](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1). The [previous discrete search](dim27_discrete_20261121/) produced another configuration of the same size. The other starting counts come from the [repository table at the parent revision](https://github.com/alexlegeartis/KissingNumbers/blob/fde77c3e8b9124640ec65f8915337960029f30d4/README.md).

## Search results

In 25D, matching on the conflict graph of two configurations gives the largest compatible subset of their union. We tested sign changes of the Leech coordinates, ranked candidate isometries by how many of the 1,999 deleted roots they preserved, and tested unions with the two preceding repository configurations. Searches over chord lifts and circle triples also reached the original count. We also tested 1,158 oriented one-dimensional sections of the tail coordinates of the 26D through 31D constructions against the 25D incumbent. Numerical screening in this dimension served only to find possible improvements.

In 26D, each proposed kernel point can be priced against the distinct projections of the exterior. Of the projections with squared norm above 32, 64 occur only twice. Releasing these constraints admits new rational directions at a cost of two exterior points per projection. We searched 40,000 local sets of deletions linking these pairs to points of the existing 340-point kernel. This yielded no verified improvement over 201,204. The [26D parent](../verifications/improved/dim26-section-rational-kernel/) records the exact rational kernel.

In dimensions 27 through 31, we pooled candidate points from contact sums, midpoint lifts, and existing configurations. Retaining candidates with small sets of incumbent conflicts gave pools of 58,928, 198,288, 143,110, 201,454, and 538,944 distinct points, respectively. We used bipartite matching to locate sets of candidates with fewer incumbent blockers than candidates. We then enforced the conflicts between candidates in those sets. Searches covered selected components and sets with a common octad, complement, or tail projection.

For 29D we enumerated all 1,104,376,728 even integer vectors of squared norm 128 in the octad and five tail coordinates. Pricing against the current 232,157-point parent yielded 736 additional candidates with at most six blockers beyond the joint pool. In 31D the corresponding 15-coordinate section contains 17,520,853,380 even vectors of squared norm 128. It supplied five additional candidates, each with six blockers. The searches on both expanded pools retained the starting counts.

In 30D the first twelve-for-eight exchange made additional trades possible. Recomputing conflicts after successful exchanges, exploring equal-size configurations, and combining trades against the original parent gave the final 340 insertions and 250 deletions. The exact optimizer recovers the gain of 90 on the supplied 1,403-candidate instance. The complete pair check certifies the lower bound independently of choices made during search. The reported ties in other dimensions leave their global optima open.
