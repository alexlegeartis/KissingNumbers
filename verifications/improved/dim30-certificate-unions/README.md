# 30D certificate unions: 248,212

The complete coordinate certificate proves **K(30) >= 248,212**, **32 above**
the preceding [248,180-point native design](../dim29-30-native-partial-design/)
and **220 above** PackingStar's pinned November 11 release of 247,992.
This is a lower bound, without an optimality or exhaustive priority claim.

The improvement combines useful points from older certificates with the newest
construction. First, an exact union with the retained
[Barnes–Wall certificate](../dim30-bw14-module/) restores **24 points without
any deletion**, reaching 248,204. A second union with the historical
[joint refill](../dim27-31-joint-refills/) adds 28 and removes 20, gaining
another eight. Relative to the original 248,180-point source, the final
exchange has **52 additions and 20 removals**.

For two valid codes, only cross-code conflicts need resolving. A maximum
bipartite matching and its alternating paths give a compatible union of
maximum cardinality within those two supplied codes. Shared source deletions
are counted once. No optimizer or matching result is a premise of the final
geometric proof.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python negative_controls.py
```

Python and NumPy suffice. `verify.py` checks the SHA-256 manifest, reconstructs
the exact source exchange, checks that every added point comes from a credited
donor, and independently checks every added-point interaction with int64
coordinate sums. `verify_geometry.py` then checks **all 30,804,474,366 unordered
pairs** of the final configuration, without assuming the source is valid.
All points are rational. Stored integral numerators have squared norm 128,
distinct-pair products at most 64, and coordinate denominator 2. Binary32
matrix products are exact because all absolute partial sums are at most 128.
Six semantic corruption controls in `negative_controls.py` exercise the checker.
`data/verification30.json` records the complete independent verification.

`final.py` reads the frozen count for the registry. Optional
`python reproduce_union.py` repeats both exact exchanges using `union_engine.py`
and the retained source/donor arrays; this discovery replay also requires SciPy.
The broader campaign and bounded unsuccessful searches in 25–31D and 35D are
recorded in `../../../research/joint_orbits_25_35_20261114/`.

PackingStar's original coordinates and the earlier repository constructions
remain credited in their retained packages. The pinned external release is
commit `2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`.
