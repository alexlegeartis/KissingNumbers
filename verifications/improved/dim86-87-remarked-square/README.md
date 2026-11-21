# Re-marked QR48, axial Pless and square layers in 86–87D

| Dimension | Dense | Axial Pless | Squares | Full root shell | Lower bound | Gain |
|---|---:|---:|---:|---:|---:|---:|
| 86 | 7,346,640,384 | 661,606,290 | 284,544 | 14,620 | **8,008,545,838** | +339,408 |
| 87 | 7,346,640,384 | 814,042,080 | 284,544 | 14,964 | **8,160,981,972** | +349,026 |

Run `python verifications/improved/dim86-87-remarked-square/verify.py` with NumPy and assertions enabled. The exact replay checks the actual marked QR48 projection, pinned Pless orbit and phase counts, Steiner triple-incidence certificate, and inequalities for the two extra layers. The standard QR48 and Pless minimum-distance-15 facts are inherited from the [axial Pless parent proof](../dim86-87-axial-pless/).

## First two layers

Re-mark the chord coordinates separately in each dimension. The zero-based chord positions in 86D are **11, 15, 27, 29, 32, 37, 42, 45**, with coordinate 0 fixed to zero and discarded in both dense and sparse. In 87D they are **6, 8, 11, 26, 29, 41, 42, 44**, with the extra axis representing coordinate 0 in the sparse layer and identically zero in the dense layer. Each nine-column QR projection has rank nine, so the two-sided dense count is unchanged. The marked Pless histograms and exact integer weighted phase counts are in `claims.json`. Each point in the first two layers has squared norm 45, and every full plane or single axial block has squared norm at most 3. The re-marking replaces the old point sets rather than adding a second copy of either.

## Square and root layers

Use the square-vector construction from the [90–91D A2 square and root proof](../dim90-91-a2-rotated-square/). In each full plane rotate the orthonormal basis so the four unit square corners are at least 15 degrees from every signed A2 root. For every four-plane support take all even-parity signs on its eight coordinates, each of absolute value \(1/\sqrt{2}\), giving 128 vectors of squared norm 4. On a support, distinct signs differ in at least two positions and hence have dot at most 2. Distinct supports intersect in at most two planes, which gives the same bound. A square has dot at most 4 with dense and at most \(4\sqrt{3}\cos(15^\circ)\) with sparse. The latter squared is \(24+12\sqrt{3}<45\); both are below the required threshold \(\sqrt{45}\).

Use the 2,223 quadruples remaining after removing one point from a 40-point Steiner quadruple system. Its exact construction uses [Hanani's \(3n-2\) recipe](https://doi.org/10.4153/CJM-1960-013-3) on the [La Jolla](https://dmgordon.org/ljcr/) 14-point system, in [Sage's relabeling](https://github.com/sagemath/sage/blob/develop/src/sage/combinat/designs/steiner_quadruple_systems.py). The 91 explicit input blocks are in `sqs14.txt`; replay checks all triple incidences at both sizes. The support count is optimal within the triple-disjoint quadruple model on 39 labels: each pair lies in 37 triples, and a packing covers an even number of them, so it must leave at least \(\binom{39}{2}/3=247\) triples uncovered in total.

Replace the old axial-only root tail by the full \(D_n\) shell in this orthonormal basis. Distinct signed two-coordinate roots have dot at most 1; each has squared norm 2. A root meets at most two original blocks of norm squared at most 3, so its mixed dot squared is at most \(2\cdot6<2\cdot45/4\). With a square its dot is at most \(2/\sqrt{2}\), the required threshold. Occupancy separates the layers: all 39 full planes for dense; at least six for sparse; exactly four for square; at most two for roots. Normalize each point to unit norm after combining the four layers.
