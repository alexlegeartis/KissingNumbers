# 62D: a fourth spike prefix via a Kerdock subcode

The new exact kissing-number lower bound is **173,733,268**, an increase of **15,360** over the live [three-prefix predecessor](../dim62-dg-three-prefixes/). Every non-spike layer is unchanged, including the already absent ordinary coordinate poles. Refit the existing three full spike prefixes and add a fourth thinner prefix.

All spikes use the same assigned distinct 60 mother cubics, one per ordinary site. The last column below gives the numerator of the body spike and two real tail coordinates with common denominator 10000.

| Prefix | Quadratic family | (body spike, tail x, tail y) numerator |
|---|---|---|
| 5 | full DG | (51692,3848,−7524) |
| 9 | full DG | (52444,6667,11547) |
| 12 | full DG | (52444,6667,−11547) |
| 13 | 32-word Kerdock subcode | (64829,3721,0) |

The 32 quadratics are the b=0 slice of the explicit 1024-word DG generator in its pinned order. The checker enumerates all affine differences and proves minimum distance 28 within the thin bank, allowing its larger spike; all inter-prefix pairs at a common spike use the DG distance 24 (if the quadratic is identical, a nonzero affine difference has distance at least 32). Different spike sites use distinct mother cubics and thus preserve the parent distance-16 comparison. Each quadratic admits eight affine choices on the marked prefix, so the new population is 60 × 32 × 8 = 15,360. The three heavy banks maintain their previous counts. All dense, guarded sparse, signed plane and tail-axis interfaces are checked exactly using the preceding verified bank categories. Every normalized pair meets the required angle bound.

Run `python verify.py` (Python 3 with assertions enabled). From the repository it automatically replays the pinned parent if present; to require the complete predecessor chain, run `python verify.py --siblings /path/to/verifications/improved`. When unpacked without the sibling tree, the checker performs the full *incremental* exact geometry replay against the pinned, already verified source proof and count.
