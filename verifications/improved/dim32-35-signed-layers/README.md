# Signed octad layers and coordinated heads in dimensions 32–35

Checkpoint: **7 November 2026**. All four lower bounds have two independent
exact verification routes.

| Dimension | Previous reviewed bound | This construction | Added points | Previous method |
|---|---:|---:|---:|---|
| 32 | 352,320 | **352,992** | **672** | [Lin's support code and ERS layers](https://github.com/QH-Lin/kissing-number-32) |
| 33 | 370,240 | **370,408** | **168** | [Lin's support code and ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) |
| 34 | 390,980 | **390,992** | **12** | [Lin's support code and ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) |
| 35 | 415,308 | **415,334** | **26** | [Lin's support code and ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) |

The preceding constructions use 1,713, 1,852, 2,013 and 2,202 supports from
[Brouwer's table](https://aeb.win.tue.nl/codes/Andw.html), reviewed on 26 October.
Their counts are 131,072 + 128 M + 2n(n−1). This is the gain over those public
witnesses. A complete new retrieval of Cohn's table was unavailable.

## Non-antipodal 33D octads and individual root-star trades

The current 33D bound is **370,408**, 32 above the actual published 370,376
repository certificate and 168 above the pinned Lin construction. The
[self-contained 33D package](nonantipodal33/README.md) contains 131,072 dense
points, 1,837 full-even octad layers, 2,064 individual exceptional octads,
1,956 retained roots and 180 individual stars. All have squared norm 32,
ordinary products at most 16, and rank 33. Exactly 500 points lack antipodes.
Dense signs occupy coordinates 0 through 31 inside R^33; stars and the 156
root deletions are individual signed vectors, never implicitly doubled.

The unpublished 370,392 intermediate is preserved as provenance. Its octad
exchange gains 16 over 370,376, then the individual star/root trade gains
another 16. Independent reconstruction checks the actual old certificate,
all coordinate-set deltas and all pair categories. Both intermediate and
final audits have 16 corruption controls. `nonantipodal33/replay.py` checks
the complete manifest and reruns every proof in a temporary directory.
No global optimum, source-union optimum or worldwide priority is claimed.

## Non-antipodal signed points: 7 November 32D continuation

The current 32D bound is **352,992**, a gain of 160 over the reviewed 352,832
repository certificate and 672 over the pinned Lin construction. The
[self-contained continuation package](nonantipodal32/README.md) keeps all
1,704 full-even layers and replaces the exceptional set by 1,824 individual
signed points on 29 supports. Exactly 380 points lack antipodes. Two independent
implementations verify every pair category; the audit checks every exceptional
pair both by scalar masks and int64 coordinates. Fourteen semantic controls and
one hash-corruption control protect the explicit-point interpretation.

The preceding fixed-25-support capacity bound does not apply to this new
29-support set. Its 14 affine and 15 nonlinear fibers are properties of this
witness, not a necessity theorem. No worldwide-priority or global-optimality
claim is made. Every preceding certificate and its proof remain unchanged.

## Historical antipodal signed octads: earlier 7 November continuation

The preceding 32D bound was **352,832**, 32 above the preceding repository certificate
and 512 above the pinned Lin construction. [The new self-contained proof](exceptional32/README.md)
replaces 380 exceptional points with 412 on the same 25 supports. Every other
layer remains unchanged. Three exact routes cover all pair categories; the
standalone `independent_audit.py` materializes all 352,832 coordinate rows and
checks every exceptional pair directly. `replay.py` also runs the two original
checkers, `verify_delta.py`, `verify_pool_bound.py` and `negative_controls.py`.
The exact 1,664-point capacity bound applies only to these 25 supports, not to
other supports, coordinate shapes or the kissing number. Fresh source snapshots
and unavailable Cohn retrieval are recorded; worldwide priority is not asserted.
The historical 352,800-point certificate remains byte-for-byte unchanged.

## Rational heavy-head rotation: 6 November continuation

The 34D bound is now **390,992**, two above the preceding repository
certificate. Moving the dense layer to coordinates excluding 0 and 2,
switching one star parity, and rationally rotating one heavy head restore
two roots. [The exact construction](rotated34/README.md) has 390,992 unique
integer vectors after common scaling, all norm squared 5,780,000 and all
pairwise products at most 2,890,000. No numerical premise is used.

`rotated34/verify.py` supplies a symbolic integer proof;
`rotated34/verify_coordinates.py` independently expands every layer,
checks all interaction types, uniqueness and rank, and optionally writes
the full coordinate archive. `rotated34/negative_controls.py` checks ten
corrupt-certificate rejections. The parent `verify.py` runs all three in
an isolated copy; the historical source certificates remain unchanged.

## Shared root deletions: 28 October continuation

The 28 October bounds improve the preceding repository certificates by **8, 8
and 20 points** in 33D, 34D and 35D. Four, four and eight complete root-star
families add 64, 64 and 128 points while removing 56, 56 and 108 roots.
All preceding dense, middle and head points are retained. Shared deleted
root pairs make these trades profitable. The [root-star proof](root_trades/README.md)
states the exact coordinate rules and both independent verification routes.

`root_trades/verify.py` proves the sign maxima symbolically;
`root_trades/verify_coordinates.py` independently expands and checks every
new/base product. `root_trades/negative_controls.py` rejects duplicates,
undeleted roots and false counts in both routes. `root_trades/search.py`
and `root_trades/enumerate_stars.cpp` are optional discovery tools; no solver
is used by the proof. The complete manifest covers every new file.

## Retained 26 October construction and proof

In 32D, 1,704 complete even-sign octad layers and 816 exceptional antipodal
lines replace a small part of the support packing. In 33D there are 1,848
complete layers and 320 exceptional lines. Their counts are
131,072 + 218,112 + 1,632 + 1,984 = **352,800** and
131,072 + 236,544 + 640 + 2,112 = **370,368**.

In 34D a distance-six coset of the dense binary code permits four full-support
heads of squared norm 50. Delete one antipodal root pair, giving
131,072 + 257,664 + 2,242 + 4 = **390,982**.

In 35D, coordinate the parity of all 2,202 octad supports with three head
families. Six antipodal lines give twelve heads of squared norm 51; three
antipodal root lines are removed. The resulting count is
131,072 + 281,856 + 2,374 + 12 = **415,314**.
Exactly 59 support parities are odd. The maximum absolute head products
against dense, middle and retained root vectors are 20, 20 and 16; the
maximum between non-antipodal heads is 17. Normalized comparisons use
4 dot² ≤ norm1 × norm2, with integers throughout.

`source/classic/proof.txt` proves the signed replacements and four-head
extension. `source/multi35/proof.txt` proves the stronger twelve-head witness.
The earlier 415,310-point 35D witness is retained and explicitly checked as
history. No construction optimum or kissing-number optimum is claimed.

## Reproduction

The current 33D package runs `nonantipodal33/replay.py`, which includes the
scalar and coordinate checkers and the independent intermediate/final audits.
The independent replay entry points are `nonantipodal33/independent_core.py`
and `nonantipodal33/independent_audit.py`;
`nonantipodal33/evidence/original_independent_core.py` and
`nonantipodal33/evidence/original_independent_audit.py` are frozen provenance
copies, not replay entry points.
It preserves and rechecks the previous 370,376 certificate as the actual
published control. See its README for the individual root/star convention.

The current 32D package runs `nonantipodal32/replay.py`,
`nonantipodal32/verify_coordinates.py`, `nonantipodal32/independent_audit.py`
and `nonantipodal32/verify_fibers.py`. Its
`nonantipodal32/evidence/original_independent_audit.py` is a frozen provenance
copy, not a replay entry point. See that package's README for its standalone
command, individual-point schema and complete proof.

```bash
python final.py
python verify.py
```

`verify.py` checks the complete SHA-256 manifest, then runs both proof
representations in an isolated temporary copy. Python 3.10+ and NumPy suffice.

- `source/classic/verify.py` uses standard-library masks, exact code enumeration
  and signed-support bounds for the original four certificates.
- `source/classic/verify_coordinates.py` independently expands every middle
  point, root and added head; all exceptional pair products use integer arrays.
- `source/multi35/verify.py` uses standard-library parity formulas and exact
  code enumeration for the twelve-head 35D certificate.
- `source/multi35/verify_coordinates.py` independently checks all 4,983,624
  head/base products and all head/head pairs in integer coordinates.
- `source/multi35/source/multi_head.py` deterministically regenerates the
  stronger certificate when passed the supplied
  `--baseline source/original_certificate35.json --output certificate35.json`.
- `source/public/` contains the downloaded support witnesses, pinned hashes,
  attribution and dense-code provenance. `final.py` checks the witness hashes
  and derives every count. Original proof receipts are retained.

Both checkers regenerate the [32,17,8] Cheng–Sloane code and its minimum
distance. Support intersections and coordinate shapes prove the large base
pair families; every exceptional case is checked. No search, cache, external
download, floating-point tolerance or inherited unverified population is used.
`data/SHA256SUMS.txt` covers every package file.

## Attribution

The public support codes are due to Qihong Lin (October 2026). The binary
generator is Cheng–Sloane's [32,17,8] code, supplied in the archived Qiushi
certificate; `source/public/dense_code.json` preserves that provenance.
The layered construction is [Edel–Rains–Sloane](https://arxiv.org/abs/math/0207291).
Signed replacements also appear in [Yang et al.](https://arxiv.org/html/2609.35051v1).
The new bounds concern the four explicit configurations above.
