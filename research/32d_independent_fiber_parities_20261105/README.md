# Independent parity choices for the 32–35D signed-octad constructions

This is a new search space, **not a new lower-bound record**. Each full octad
fiber contains all 128 even-parity sign patterns. Any subset of the following
fibers can instead use all 128 odd-parity patterns, independently, without
changing cardinality or violating the kissing constraint:

| Dimension | Certified points | Independently flippable fibers |
|---:|---:|---:|
| 32 | 352,800 | 1,704 of 1,704 |
| 33 | 370,368 | 1,848 of 1,848 |
| 34 | 390,982 | 1,878 of 2,013 |
| 35 | 415,310 | 2,063 of 2,202 |

Thus, for example, the 33D source has 2^1848 distinct coordinate
configurations. No claim is made that they are pairwise nonisometric.

Each fiber has norm squared 32 and, within one parity, distinct sign patterns
differ in at least two positions, giving inner product at most 16. Different
full supports meet in at most four coordinates; every exceptional support
also meets every full support in at most four. The dense sign and two-coordinate
root layers each contribute at most 16 to a full-fiber product. The 34D and
35D certificates include extra heavy heads, so each fiber's allowed parities
must additionally pass the exact normalized product test against all four
heads. The stated counts are the fibers for which both parities pass.

`verify_parity_freedom.py` first replays the existing complete source
certificates, then checks all support intersections, within-fiber parity
inequalities and (where applicable) heavy-head inequalities. The statements
do not depend on numerical optimization or tolerances.

The engine operator `search.signed_fiber_exchanges` can choose parities for
selected full fibers while jointly deleting old antipodal lines and adding
individual signed octad lines. Its output is independently replayed by the
compressed signed-support verifier. This enables searches to break complete
orbits, but does not promise a gain. The current engine's compressed backend
represents the ordinary layers; a source containing extra heavy heads requires
additional compatibility checks or discarding those heads.
