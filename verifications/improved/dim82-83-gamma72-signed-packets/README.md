# Independently signed Γ72 augmentations in 82D and 83D

| Dimension | Verified lower bound | Previous repository bound | Gain |
|---:|---:|---:|---:|
| 82 | 6,269,084,390 | 6,269,068,246 | 16,144 |
| 83 | 6,278,756,940 | 6,278,746,180 | 10,760 |

We keep the 6,218,175,600 contacts of Nebe's Γ72 and the 510-point 10D [Ganzhinov](https://arxiv.org/abs/2207.08266) and 604-point 11D [Bianchi et al.](https://arxiv.org/abs/2606.10402) tail codes from the previous release. The tail code by itself stays antipodal. Its lifted horizontal choices now depend on the **signed vector**, not only on its antipodal ray. We independently select which complete horizontal sign packet is retained at each tail sign. Pure tail vectors are unchanged.

All mixed points have horizontal and tail squared norm six and two, respectively. The pinned Γ72 package supplies exact horizontal packet coordinates and half-phase separation: different first colors clear squared distance six; the difference of different second phases of the same color can have squared distance four in exactly the singleton raw norm-eight block case. Other distances for such pairs are at least six. The only adverse signed-angle case is therefore *positive* tail inner product. Equal phases require different horizontal points unless the tails are antipodal. Because the latter have squared separation eight, their horizontal choices need not be disjoint.

Every potentially conflicting pair is replayed exactly by packed class intersections and integer Leech pairings, including each within-profile conflict. A singleton raw-block change preserves the sorted norm profile, so the 4+4+4, 4+8, and 6+6 packets compose. The 4+4+4 packets contribute eight points at a single signed tail; the two others contribute four. Each mixed point clears the Γ72 equator by its depth-six phase and its tail squared norm two. Each pure point has tail equal to twice a norm-two signed tail vector, so its dot with a mixed point is at most four by Cauchy–Schwarz, and pure-pure dots are at most four by the tail kissing bound.

Install NumPy and SymPy. Place this directory under `verifications/improved/` next to the pinned Γ72 siblings, then run `OPENBLAS_NUM_THREADS=1 python verify.py` with assertions enabled. An explicit sibling path is also accepted as the first argument. The local witness arrays are literal packet masks; the verifier neither needs the discovery solver nor relies on optimality of the search.
