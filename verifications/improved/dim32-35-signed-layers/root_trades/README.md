# Shared root-star deletions in dimensions 33–35

28 October 2026. Each new bound has a separate formula proof and an expanded
integer-coordinate proof, freshly replaying its retained baseline first.

| Dimension | Preceding repository bound | New bound | Added star points | Removed roots | Gain |
|---|---:|---:|---:|---:|---:|
| 33 | 370,368 | **370,376** | 64 | 56 | **8** |
| 34 | 390,982 | **390,990** | 64 | 56 | **8** |
| 35 | 415,314 | **415,334** | 128 | 108 | **20** |

For a center a and four leaves B, take the 16 vectors
`+/-4 e_a + 2*sum(s_b e_b)`, with an even number of negative leaf signs.
Their squared norm is 32. Delete all four roots on every center/leaf coordinate
pair. One star alone adds and removes 16 points. Compatible stars can share
deleted root pairs: each shared coordinate pair saves four deletions.
The certificates choose four, four and eight stars, respectively.

The original dense code, middle vectors and heads are retained. For every
new vector, `verify.py` maximizes exactly over each full middle support's sign
patterns, checks all exceptional signed lines, retained roots and normalized
heads, and uses the projected l1 bound against the dense sign code. Every
new/new pair is checked. Retained old/old pairs remain valid by the freshly
verified baseline. The head comparisons use `4*dot^2 <= 32*head_norm`.

`verify_coordinates.py` independently expands each star and all retained base
layers. It checks every new/base coordinate product: **23,699,968** in 33D,
**25,019,264** in 34D and **53,146,368** in 35D. Mixed-norm squares use int64.
Both routes verify distinctness and derive the counts from actual deletions.
All baseline inputs are bound by SHA-256; no optimizer output is trusted.

`negative_controls.py` requires both proofs to reject duplicate stars, undeleted
roots and false counts in each dimension: eighteen independent rejection checks.

From the package root, `verify.py` replays both routes in a temporary copy.
The original `source/classic/` and `source/multi35/` witnesses remain unchanged.
This extends the repository's 39D root-star mechanism; it is not an optimality
claim for these dimensions or for the whole construction family.

`search.py` and `enumerate_stars.cpp` are optional discovery tools. The C++
enumerator visits all center/four-leaf choices compatible with the fixed full
middle supports. Python filters exceptional lines and existing heads, then
uses a binary MILP to select compatible stars and shared root deletions.
This finite pool contains 28, 41 and 31 stars in the three dimensions.
Python, NumPy, SciPy and a C++17 compiler suffice for the optional search:

```bash
python root_trades/search.py 35 --output /tmp/root-trade35.json
```

Solvers may select different witnesses with the same gain. The explicit
shipped certificates, independently checked by both proof routes, determine
the accepted bounds. The proof replays require Python and NumPy only.
