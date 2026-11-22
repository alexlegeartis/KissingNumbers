# Joint signed refill in dimension 33

The 33D lower bound is **370,528**, up **112** over the preceding repository certificate and 288 over Qihong Lin's 370,240. [Previous repository proof](../dim33-octad-switches/) · [Lin's construction](https://arxiv.org/abs/2611.03682)

| Layer | Signed points |
| --- | ---: |
| Dense binary sign code | 131,072 |
| Complete sign families on 1,837 octads | 235,136 |
| Individual signed octads | 2,216 |
| Signed stars | 172 |
| Signed roots | 1,932 |
| **Total** | **370,528** |

The dense code and the set of complete octad supports come from the 370,416-point parent. The exceptional octad layer grows by 152 points and the star layer shrinks by 40. The same 1,932 roots are present in the final configuration.

## Sparse search

Enumerate the eight-coordinate supports meeting every complete octad in at most four coordinates. Generate both parities of sign words on these supports and keep individual points as vertices. Generate each signed star with magnitude four at its center and magnitude two at four leaves. A star can be used on the fixed octad supports precisely when no complete octad contains the center and at least three of its leaves.

A candidate is priced by the number of its conflicts with the explicit part of the parent. The frozen search instance uses a cutoff of fourteen. Vertices include each eligible new point and every incumbent point affected by a candidate. Put an edge between two vertices exactly when their ordinary inner product exceeds 16. This allows stars, signed octad points and roots to participate in one independent set calculation. The first 44,647-vertex instance includes 2,400 points of the incumbent and 42,247 challenges. The first solution selects 2,504 vertices and retains 1,808 incumbent points outside the graph, giving 370,520. Repricing against that configuration gives a second graph of 17,730 vertices, with 2,466 affected incumbent points. Its solution selects 2,474 vertices and retains 1,846 points outside the graph, increasing the count by eight.

The certificate is independent of the optimizer. Rebuild the two graphs with `python search_graph.py /tmp/dim33-initial.graph` and `python search_graph.py /tmp/dim33-repeat.graph --stage repeat`. The archive stores graph coordinates in vertex order as `v`, the number of affected incumbent points as `start`, and the unaffected points as `extra`. The discovery used ReduMIS with seed 1703. [Initial search instance](data/initial_search.npz) · [Repeat search instance](data/repeat_search.npz)

## Verify

```sh
OPENBLAS_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=4 python verify.py --all-pairs
```

The first replay checks the complete geometry by pair categories. It enumerates the dense binary differences, checks the exact support intersections of the complete octads, and computes all products among the explicit points in signed 64-bit arithmetic. For a complete parity family on an octad, the maximum product with an explicit point is twice the sum of magnitudes on that support. If the maximizing signs have the wrong parity, subtract four times the smallest magnitude. This includes every mixed pair with a complete octad. Each incoming sparse point has total absolute magnitude at most 16, giving the dense cross-layer bound.

The optional expanded replay covers all **68,645,314,128 unordered pairs** and finds maximum inner product 16 at common squared norm 32. The registry entry point is `final.py`. [Witness](witness.npz) · [Verification receipt](verification.json)
