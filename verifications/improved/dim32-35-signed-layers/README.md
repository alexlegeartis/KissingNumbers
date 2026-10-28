# Signed octad layers and coordinated heads in dimensions 32–35

Checkpoint: **28 October 2026**. All four lower bounds have two independent
exact verification routes.

| Dimension | Previous reviewed bound | This construction | Added points | Previous method |
|---|---:|---:|---:|---|
| 32 | 352,320 | **352,800** | **480** | [Lin's support code and ERS layers](https://github.com/QH-Lin/kissing-number-32) |
| 33 | 370,240 | **370,376** | **136** | [Lin's support code and ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) |
| 34 | 390,980 | **390,990** | **10** | [Lin's support code and ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) |
| 35 | 415,308 | **415,334** | **26** | [Lin's support code and ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) |

The preceding constructions use 1,713, 1,852, 2,013 and 2,202 supports from
[Brouwer's table](https://aeb.win.tue.nl/codes/Andw.html), reviewed on 26 October.
Their counts are 131,072 + 128 M + 2n(n−1). This is the gain over those public
witnesses. A complete new retrieval of Cohn's table was unavailable.

## Shared root deletions: 28 October continuation

The current bounds improve the preceding repository certificates by **8, 8
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

```bash
python final.py
python verify.py
```

`verify.py` checks the complete SHA-256 manifest, then runs both proof
representations in an isolated temporary copy. Python 3.9+ and NumPy suffice.

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
