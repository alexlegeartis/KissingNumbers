# Three omitted prefixes in the 62D Delsarte–Goethals construction

**173,717,908** contacts, a gain of **982,920** over the preceding
[DG spike packing](../dim62-63-dg-spike-banks/).

Keep all old dense words, all rephased sparse sign banks (including
the short-guard additions), all signed planes, and the four axial tail
poles. Delete the 120 ordinary coordinate poles. For each of the 60
spike sites use the same assigned cubic and the same DG sign code as
before. Replace the old donor by the following three families,
with bit indices in (0,9,20,41) order as in the predecessor:

| Prefix as integer | Body spike magnitude | Tail coordinates |
|---:|---:|---:|
| 5 | 5 | (1/5, −3/10) |
| 9 | 29/5 | (1, 7/5) |
| 12 | 29/5 | (1, −7/5) |

Each prefix gives eight affine shifts for each of the 1,024 quadratic
choices, for 1,474,560 donor points in all. They replace the earlier
491,520 donors and the 120 removed poles. No donor is identified by
puncturing, by the original distance bound.

The 24-distance DG bound applies to donors assigned to the same
site; the 16-distance parent bound applies at different sites and
against the dense hexagon. In each case the number of lost distances
at the prefix is counted before the exact norm-scaled angle test.
Every old sparse bank uses a fixed total sign parity orthogonal to the
mother code; the replay accounts for the sign disagreement it forces
for each prefix, as well as all possible single- and double-mark head
contributions. The retained planes and tail axes obey direct bounds.
The only square root in the dense acute angle check is bounded by
1732/1000 < sqrt(3) < 1733/1000.

Run python3 verify.py without -O. It hash-pins and fully replays the
predecessor. An optional argument selects an alternate directory
containing sibling certificates.
