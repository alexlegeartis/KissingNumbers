# Shared phases and independent horizontal triples in 84D and 85D

| Dimension | Verified lower bound | Gain over the preceding repository bound |
|---:|---:|---:|
| 84 | **6,292,388,112** | 1,218,424 |
| 85 | **6,295,120,802** | 958,168 |

These constructions still use Nebe's Γ72 core. In 84D the
[Best support tails](../dim81-84-antipodal-tail-codes/) already supplied
420 pure antipodal rays, but only 384 appeared with both a tail and
a horizontal displacement. We now use all 420 for both kinds of
contacts. In 85D we start from the pinned
[Barnes–Wall equality section](../dim82-85-bw16-tail-sections/) and
retain the [same assignment of the 69 leftover rays](../dim85-shared-phase-refill/)
to shared phases.

There are two steps. First, in each added tail packet fill every
compatible signed group in the horizontal profiles 4+4+4, 4+8,
and 6+6. For rays sharing the same colour and phase, only
coincident horizontal points can conflict. For separate phases
of the same colour, an offending pair has nonzero tail dot and
a singleton horizontal difference of raw norm eight. The
other horizontal blocks must stay identical, so contacts from
different profiles cannot conflict. Different colours are safe at
the known depth-six gate. The verifier compares the 4+8 and 6+6
groups against the **actual** occupied masks, testing norm-eight
differences by exact Gram entries; its binary64 products have
integer summands and a checked bound strictly below (2^{53}).

The second step reallocates triple groups. Two rows within one
colour and phase conflict exactly when all three local B labels
are equal, even when their tails are orthogonal. Across distinct
phases of one colour, the conflict instead requires nonzero tail
dot and precisely two equal B labels. In either case the rows
share the A label and the first two B labels modulo the common
phase subspace. These invariants partition the conflict graph.
The bounded clique search provides feasible sets; its full
literal edge replay verifies every selected component. Some
previous triple coordinates are reassigned, so this second step
is not append-only.

In 84D the new packets add 49,188 triple and 399 two-block groups;
reallocation adds a further 26,764 triple groups. In 85D the
previous 86,558 added triples and 138 pure contacts are reproduced
before reallocation; the new profile adds 2,331 two-block groups
and reallocation adds 58,720 triple groups. Triple groups carry
16 points each and two-block groups 8.

Replay with assertions enabled, Python, NumPy, SymPy and a
little-endian C++17 compiler:
```sh
OPENBLAS_NUM_THREADS=1 python3 verify.py 84
OPENBLAS_NUM_THREADS=1 python3 verify.py 85
```
