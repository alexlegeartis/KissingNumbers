# Sparse filling and parity-compatible central sections

This continuation improves six registered bounds. Comparisons are with commit
`9d38a38b93fcbbc88d48981f382464569189c3c8`; they are not a new worldwide
literature survey. The parent codes retain their existing distance proofs.

| Dimension | Bound | Gain over preceding certificate |
|---|---:|---:|
| 61 | 72,376,856 | 1,066,368 |
| 62 | 100,437,548 | 66,080 |
| 63 | 182,537,988 | 130,560 |
| 66 | 338,302,340 | 10,240 |
| 67 | 351,420,806 | 20,480 |
| 68 | 371,552,152 | 30,720 |

## The weight-eight bridge in 61D

`bridge61.npz` exhibits 8,331 weight-eight supports on 61 coordinates. Every
support occupies at most five of the fifteen four-cell rows in coordinates
0 through 59. Coordinate 60 is the latitude coordinate. Any two supports
intersect in at most four positions.

Give each support the 128 even-parity sign patterns, with nonzero magnitude
1/sqrt(8). Two signs in the same support differ in at least two positions,
and different supports overlap in at most four, so their inner products
are at most 1/2. The added layer has 8,331*128 = 1,066,368 points.

Each old weight-15 support has one cell in every four-cell row. Its
intersection with a new support is therefore at most five, and
25/(8*15) <= 1/4 proves the mixed-angle bound after squaring. For the old
weight-two layer the corresponding bound is 4/(8*2) = 1/4.

The existing dense compression has head/tail magnitudes (1/8,1/4) or
(sqrt(17)/32,1/16). Without the latitude coordinate, the squared mixed
products are bounded by 8/64 and 136/1024. With it, they are bounded by
81/512 and (837+28sqrt(17))/8192. Each is at most 1/4, checked exactly.
Nonzero-coordinate counts distinguish every new point from the old
61-coordinate dense, 15-coordinate and two-coordinate points.

`bridge61.py` checks all 34,698,615 support pairs using packed integer
population counts and independently using binary incidence coordinates.
The latter products are exact integer sums bounded by 61. It enumerates
the sign code, checks all mixed inequalities with rational/quadratic
arithmetic and rejects eight malformed-certificate controls.

The seed has 2,625 supports: restrict the affine SQS(32) to the thirty
nonzero-row cells, double them, and use three pairings of each four-cell
row. `search/seed.py` constructs and verifies it. `search/fill_correct.cpp`
performs integer-only insertions, one-for-two replacements and neutral
swaps. The recorded run uses seed 2026102804, 150 seconds, 87,544,578
accepted trials and 374,815 neutral swaps. Its elapsed-time stopping rule
does not promise a byte-identical rediscovery on another machine. The
exhibited certificate is verified independently of the discovery search.

```bash
python search/seed.py /tmp/bridge_seed.txt
g++ -O3 -march=native search/fill_correct.cpp -o /tmp/bridge_search
/tmp/bridge_search /tmp/bridge_seed.txt /tmp/bridge_result.txt 150
python bridge61.py
```

## Parity-compatible sections in 62D and 63D

`parent64.json` is an exact copy of the signed-support package's
`data/construction.json`. Its 211 extra supports each conflict with eight
different old supports in ten coordinates. All 1,688 affected old supports
are distinct. Every conflict supplies an eight-coordinate common parity
check. Opposite parities force a sign disagreement and bound the product
by 4*(10-2) = 32 at squared norm 64.

The earlier sections independently optimized the signs on every old support.
The new construction uses the augmented parent. For an affected old support,
put the common eight-set first, then group positions by the disjoint section
supports. Give it the even or odd common parity opposite to its new neighbour.
`section_shifts.json` specifies the coset shift of each new support separately
for 62D and 63D. Different new supports intersect in at most eight positions,
so those shifts introduce no further constraints.

`signed_sections.py` rechecks all added/old incidences, disjoint ownership,
sign distances and parity conditions. It counts every relevant fibre twice,
using bit counts and expanded integer sign coordinates. The dense counts come
from the actual binary projection; the independent four-support layer is
relabelled to make the section sets consecutive affine 3-spaces. Central
sections inherit all norms and angles from the checked parent configuration.
`signed_sections.json` records the derived terms and gains.

The parent package's `verify.py` freshly enumerates 536,870,912 words and
checks every parity repair. It passed again for this continuation. The
new sections preserve its binary-code distance input and add no new
assumption about that code.

## Complete quadruple packings in 66--68D

`netto.py` constructs the orbit of {infinity,0,1,30} on the projective line
over F67 under determinant-one Mobius transformations, identifying matrices
up to sign. The seed satisfies 30^2-30+1 = 0 modulo 67. All 150,348
transformations are evaluated. The resulting 12,529 blocks cover each of
the 50,116 triples exactly once. An independent incidence-coordinate replay
checks every distinct block intersection is at most two.

Restriction to 66, 67 and 68 points gives 11,088, 11,792 and 12,529 supports,
replacing 10,448, 10,512 and 10,609. Each carries all sixteen signs of
magnitude four at squared norm 64. Intersections of at most two bound
products by 32. Mixed products with the unchanged dense, weight-16 and
axial layers are bounded by 16, 32 and 32. Each additional support adds
sixteen points; `netto.json` records all totals.

The 67- and 68-point support counts attain the nested Johnson packing bound
for these four-subsets. This is not an optimality statement about the
complete kissing configurations or other support sizes.

## Scope and replay

From the package root, `python verify_continuation.py` replays all three
certificates in a temporary copy and compares complete receipts. The default
`verify.py` includes these checks as well as the earlier code proofs.
`final.py` derives the new totals while retaining stronger 95--96D inputs.

The finite coordinate-shortening search also tested 170,540, 240,454 and
338,724 starts for the 66D, 67D and 68D weight-16 layers, without improving
their existing support counts. This is not an impossibility result.
No new certificate in this continuation changes 64--65D or 69--72D.
