# A 13D Zinoviev–Ericson tail for the 85D kissing construction

The verified lower bound is **6,297,176,738**, +2,055,936 over 6,295,120,802.

We replace the Barnes–Wall tail in the existing 72+13 dimensional construction by the full 1,154-point antipodal Zinoviev–Ericson code. The 577 rays have exact coordinates in Z[sqrt(3)], at squared norm 16. Source: [Henry Cohn's spherical-code data](https://hdl.handle.net/1721.1/153312), 13-dimensional 1,154-point antipodal entry; see [Zinoviev and Ericson](https://www.mathnet.ru/eng/ppi457). The literal points are in ze99.txt. The verifier checks both exact antipodality and every absolute dot between different rays.

The horizontal core remains Nebe's Γ72. This is a different tail, not a construction independent of Γ72. The change allows more rays to participate in each of the 4+4+4, 4+8 and 6+6 horizontal profiles. The first 384 assigned rays use all different colour/phase pairs; the remaining rays reuse these phases, filling compatible places. A finite component search reallocates the 4+4+4 groups and explicitly verifies that every selected pair is compatible.

The pinned sibling packages regenerate the exact horizontal groups and establish the separation conditions. For an antipodal tail, only the distinction between zero and nonzero tail inner product is needed by this conservative check. Different colours clear the depth-six condition; within one colour, coincident horizontal points are forbidden, and different phases are checked against precisely the pairs whose horizontal difference has raw norm eight. These offending differences are contained within one profile, so the three profile checks compose. Every tail produces two pure contacts. The verifier only uses the literal assignment; it does not search for it again.

Replay from this directory with Python, NumPy, SymPy and a C++17 compiler on a little-endian system; keep the assertions enabled:

```sh
OPENBLAS_NUM_THREADS=1 python3 verify.py
```
