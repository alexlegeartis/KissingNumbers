# Global ownership exchange and post-refill, dimensions 50–60

New certified lower bounds:

| Dimension | Signed ancestor | Net body exchange | Net new tail | Lower bound |
|---:|---:|---:|---:|---:|
| 50 | 52,462,872 | 111 | 0 | 52,462,983 |
| 51 | 52,509,740 | 216 | 0 | 52,509,956 |
| 52 | 52,603,402 | 437 | 0 | 52,603,839 |
| 53 | 52,728,054 | 724 | 0 | 52,728,778 |
| 54 | 52,976,623 | 1,321 | 156 | 52,978,100 |
| 55 | 53,393,546 | 2,324 | 168 | 53,396,038 |
| 56 | 54,263,398 | 4,475 | 0 | 54,267,873 |
| 57 | 54,760,693 | 5,936 | 352 | 54,766,981 |
| 58 | 56,270,076 | 11,619 | 0 | 56,281,695 |
| 59 | 56,952,414 | 14,933 | 224 | 56,967,571 |
| 60 | 58,631,571 | 27,016 | 5,135 | 58,663,722 |

A donor may replace multiple incompatible signed owners in its destination
class. A signed owner appearing in a different class may also be reassigned,
provided its old assignment is deleted. We priced the union of these
requirements with one maximum-closure computation across all classes for each
dimension, using a common bank of 4,000 candidate vectors. Next we refilled
against the resulting retained stock, favoring several compatible donors that
share a small set of surviving blockers. The 60D row includes a third
pass that permits groups sharing up to two such blockers. Verification does not rely on optimization output
beyond the literal inserted coordinates/classes and checks the exact removal
union from scratch.

The parent is the same norm-six shell of the 48-dimensional P48 lattice as in
the pinned signed ancestor. Each world head has at most one owner, and a moved
owner is removed from its original class. New heads in one class have raw
inner product at most 2 with every other new or retained head in that class.
Other classes have the shell bound of 3 on raw inner product for distinct
signed world heads. These are precisely the inherited body geometry conditions;
same class uses the same tail line, and distinct classes obey the prior
absolute tail inner-product bound of one half. All class isometries and integer
bilinear arithmetic used for the new heads are checked.

The 54D minuscule-tail proof is replayed from its pinned sibling. The
tetrahedral tail palettes in 55, 57 and 59 replace the older trine palettes,
and 60 uses the 16-color frame plus four sparse tetrahedra and 20 upgraded
cube subcosets, distributed across three syndromes. Those literal tail
geometry verifiers replay independently against the same signed ancestor.
The composition then checks the union of old and new tail reserve directions,
including both signs, against the full world coordinates of all body inserts
and predecessor owners. The stronger tail replaces, rather than adds to, the
earlier tail ledger.

Replay (Python 3 with NumPy, assertions enabled):
```sh
python verify.py
python verify.py 60
```
Optionally pass `--siblings` for the directory with the signed and pointwise
ancestor packages, and `--tail` for the complementary-tail sibling.
Bounded integer-valued matrix products in the body use float64 only where the
absolute sum of all products is strictly less than `2**53`, so comparisons
are exact. World-head equality uses all 48 signed coordinates.
