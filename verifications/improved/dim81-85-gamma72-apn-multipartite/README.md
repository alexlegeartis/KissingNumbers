# APN-coloured multipartite D_m fibres in 81D–85D

| Dimension | Exact lower bound | Increase over preceding bound |
|---:|---:|---:|
| 81 | 6,232,188,256 | 3,483,134 |
| 82 | 6,235,144,100 | 5,892,158 |
| 83 | 6,238,576,140 | 8,503,968 |
| 84 | 6,242,540,232 | 10,827,600 |
| 85 | 6,246,072,200 | 10,078,350 |

The construction keeps all 6,218,175,600 norm-eight contacts of Nebe's Γ72. Its orthogonal tail lies in D_m, where m=9 through 13. This is a Γ72-based result, not a non-Nebe construction or an optimality claim.

## Three phase colours with a common shift space

In the three-Leech model of Γ72, an integral-block half-phase has labels (a,b,t), a,b in the first polar subspace and t in the second. Set a=23214, b=668618, c=a+b. Use three first-label colours (a,b), (b,c), (c,a), together with zero. They form a two-dimensional additive group. The witness supplies a seven-dimensional space H of second labels for which each nonzero colour at each t in H has exact depth six and 136,704 shortest horizontal points. The verifier checks the root-class list exclusions and orthogonality required for the pinned [depth-six theorem](../dim73-82-gamma72-depth6-lattice/).

Label the m tail vertices by distinct x in F16. The Sidon map x↦(x,x³) gives eight bits, split into a two-bit colour from x and a six-bit second-label address. For two differently coloured vertices i,j, use both rays e_i−e_j and e_i+e_j. Assign each the difference of the vertex labels, toggling the seventh B bit for one sign. This produces distinct phases, all of depth six. The certificate checks distinctness directly.

## Why only singleton differences need clipping

The base Γ72 norm-eight shell is compatible with every tilted fibre by depth six plus tail squared length two. Distinct tail directions have absolute inner product zero or one. A pair of different phase colours has horizontal difference depth six. A same-colour pair has a nonzero B-only difference and therefore horizontal depth at least four, enough when its tails are orthogonal.

If a same-colour pair shares one tail coordinate, use one directed clip against the entire opposing fibre. The B-only short difference has first labels zero, so the three local A labels coincide. If any block vanishes, that common A label is zero, and every nonzero local class is singular with minimum at least eight. If all three blocks are occupied, the raw norm is at least twelve. Therefore the only possible raw norm below twelve is a single occupied block of norm eight. In a clipped fibre retain only raw norm-twelve points split as 4+4+4, with all signs. Such a point x conflicts with the opposing full fibre exactly when at least one of its local B labels, shifted by the difference phase, still labels a norm-four Leech class at its unchanged A label. Indeed differences in that nonzero singular local class have minimum eight, so two local roots in these classes have dot zero and yield a norm-eight difference. This test is independent of local signs. The verifier checks every directed clip and uses factor eight for the signs of the three roots. Cyclic permutation of the three Leech blocks permutes the three nonzero colours and preserves the second label, so the same count table applies to each colour. Same-ray antipodes already have vertical squared separation eight.

The pure D_m roots, scaled to squared norm eight, supply 2m(m−1) points and are compatible with all tilted fibres by Cauchy–Schwarz.

Run `OPENBLAS_NUM_THREADS=1 python verify.py`. Requires NumPy and SymPy. The replay uses only the pinned parent, literal vertex labels and clip lists, and exact integer counts.
