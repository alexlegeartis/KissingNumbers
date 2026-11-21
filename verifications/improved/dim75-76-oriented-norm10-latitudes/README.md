# Oriented polyhedral norm-ten latitudes in dimensions 75 and 76

| Dimension | Tail directions | Selected shell points | Pure tail roots | Kissing lower bound | Gain |
|---:|---:|---:|---:|---:|---:|
| 75 | 60 | 100,561 | 12 | **6,218,276,173** | +30,113 |
| 76 | 144 | 224,512 | 24 | **6,218,400,136** | +29,086 |

The 76D comparison is with the immediately preceding shaped-tail bound of 6,218,371,050. The selected shell points are proved to exist by a graph bound; the verification does not materialize their coordinate list.

For 75D take the vertices from cyclic permutations of \((0,\pm1,\pm\varphi)\), where \(\varphi^2=\varphi+1\) and \(\varphi>0\). For each directed adjacent pair \((u,v)\), take the direction of \(32u+12v+7(u\times v)\). There are 60 directions of equal norm. Exact comparisons in \(\mathbb Q[\varphi]\) show that all off-diagonal unit tail products are at most \(9/10\). For 76D take all 144 signed permutations of \((3,1,0,0)\) and \((2,2,1,1)\), the norm-ten shell in \(D_4\), and normalize. Any two distinct such vectors have squared distance at least 2, so again the off-diagonal unit products are at most \(9/10\).

For every full oriented norm-ten vector of Nebe's \(\Gamma_{72}\) and every normalized tail direction, scale the body by \(4/5\) and the unit tail by \(\sqrt{8/5}\). The result has squared norm 8. If the signed body product is at most 3, even the largest tail contribution cannot violate the dot cap 4. Otherwise for distinct shell lines the only possible positive body products are \(j=4,5,6\), with tail threshold \((25-4j)/10\). At \(j=4\) this requires an identical tail. Equal bodies conflict; opposite bodies do not. Include both body orientations as vertices and count each selected vertex once, without antipodal doubling.

The [verified parent populations, joint shell moments and rational recurrence](../dim76-85-mixed-shape-norm10/) yield the uniform local triangle caps and independence estimates. Keep all 6,218,175,600 norm-eight equatorial lattice points (their unscaled dots with a shell body are at most 5 in magnitude); add all \(D_3\) or \(D_4\) pure tail roots of squared norm 8 (cross squared dot at most \(8\cdot8/5<16\)). The layers are disjoint.

Run `python verifications/improved/dim75-76-oriented-norm10-latitudes/verify.py` with NumPy and assertions enabled. The script replays exact geometry, source arithmetic, rational dual inequalities, graph degrees, triangle caps and final integer ceilings. It inherits Nebe's minimum-8 and the Venkov shell-design premise, not solver output.
