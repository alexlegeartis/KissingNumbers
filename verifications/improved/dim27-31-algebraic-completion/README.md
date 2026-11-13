# Algebraic completion in dimensions 27 and 31

4 November 2026. Two complete explicit integer configurations improve the
preceding repository bounds:

Retained construction: the 31D certificate is superseded by the
[8 November latitude-contact exchange](../dim31-latitude-contact-exchange/).
The dated constructions and coordinate certificates below remain unchanged.

| Dimension | Preceding repository | New lower bound | Gain | Removed | Added | Current repository bound (13 November 2026) |
|---|---:|---:|---:|---:|---:|---:|
|27|209,807|**209,825**|18|14|32|209,825|
|31|251,188|**251,196**|8|0|8| 251,395 |

The PackingStar release checked for that dated construction has 209,760 and 251,169 points,
respectively. Its retained coordinates and the repository's subsequent
improvements are credited as the baselines. Both complete predecessors are
included byte for byte; their original package remains unchanged.

The 27D discovery breaks signed support orbits into individually selectable
points. Exact support/sign incidence prices the core blockers; complete
integer products find the remaining blockers. The registered
`selection.blocker_exchange` operator jointly selects the new points and
charges shared deletions once. It receives no proposed exchange or target count.

The 31D discovery starts with generated Golay sign changes and signed tail
permutations. The general `selection.isometric_exchange` operator independently
reproduces the first four insertions from a full two-copy union. Negative-contact
sums expose four more. The eight additions have a short formula, in the
norm-32 convention with zero-based indices:

```
±2(e0+e5+e6+e10) ±2(e8+e20) ±2(e24−e25).
```

The three signs vary independently. The orthogonal blocks have squared norms
16, 8 and 8, so all eight vectors have squared norm 32 and mutual product at
most 16. The independent verifier checks each against every retained point.
Stored coordinates are twice these vectors, with squared norm 128.
Together with the four retained points `±4(e8+e20)` and `±4(e24−e25)`,
the eight additions complete a scaled **A3 root system** of twelve points.

## Reproduce and verify

```sh
OPENBLAS_NUM_THREADS=1 python3 reproduce.py
OPENBLAS_NUM_THREADS=1 python3 verify.py
python3 negative_controls.py
```

`final.py` counts the arrays. `verify_geometry.py` contains the full coordinate
checker used by `verify.py`; the latter also runs the independent changed-pair
calculation.

Python and NumPy suffice. No optimizer, network retrieval or search cache is
needed. Every norm, distinctness condition and unordered pair is checked in
the complete explicit configurations. Products are exact bounded integers in
binary32; the proof checks the accumulation bound. A separate direct-int64
calculation checks every inserted point against every final point. No lattice
minimum or unenumerated population is assumed. Search and verification use
different implementations.

The original geometry verifier is preserved here without modification from
`dim27-31-joint-refills`; its hash is recorded. The complete file manifest and
deterministic reports are checked afresh. The manuscript remains its earlier
snapshot. The mathematical proofs establish lower bounds, with no unrestricted
optimality claim. Fresh Cohn retrieval failed; comparison uses the pinned
November 4 capture and the available primary author repositories.

See `research/algebraic_switches_20261104` for executable search programs,
finite pools, complete engine receipts and the unsuccessful 29D, 30D and 40D
experiments. These are seeded construction searches, separate from the
held-out rediscovery benchmark.
