# Simultaneous rephasing in the punctured 60D/61D constructions

The guarded version improves the ordinary punctured constructions as follows.

| Dimension | Before both changes | Base punctured | Final | Gain vs previous repository bound |
|---:|---:|---:|---:|---:|
| 60 | 58,665,998 | 66,657,320 | **68,033,576** | 9,367,578 (+15.97%) |
| 61 | 91,038,762 | 116,711,210 | **118,639,402** | 27,600,640 (+30.32%) |

Leave all dense points, coordinate poles, affine four-sets, and previous sparse **supports** unchanged. The earlier signs on some supports are replaced, keeping their sizes and total parity. Take 672 new zero-mark supports in each dimension, and another 539 single-mark supports in 61D. These are drawn from the published signed-overlap support witness; all new geometric and parity premises are rechecked directly, so the earlier witness is not a trusted premise.

The new weight-16 zero-mark supports carry odd Hamming sign cosets of size 2,048. Each new single-mark support uses the same 15-ordinary-site map and shortening rule as in the base 61D construction; its marked site represents the same fixed height as the old supports of that mark. The dangerous dense pairing in each marked direction needs even body sign parity, because the facing large dense label has bit one at that mark. The Hamming sign code shortened at the erased marked site retains 1,024 sign words. In both types the verifier checks even dense-parent restriction on the unpunctured support, using the explicit quadratic rows and cubic coset representatives.

Each potentially bad old/new or new/new overlap comes with a common eight-site guard with opposite sign parity. A deleted marked site, if present in the guard on both banks, is assigned a fixed positive dummy sign (binary zero) in both. It is not allowed to pay for a disagreement. Thus the guard forces a disagreement on an ordinary coordinate, reducing dot by at least two. The raw dot of every bank pair before paying for disagreements is the shared ordinary count plus the product of signed heights, bounded by 10 whenever a guard is needed. The result is at most eight on the squared-norm-16 sphere.

All guard conditions and fixed total parity are enforced simultaneously. A bank's local homogeneous checks span a space of rank at most five containing the all-ones row, with every other nonzero row of weight eight. Quotient by the all-ones row, and take any independent basis of the quotient. Its restriction labels the sixteen columns with equally sized classes (every nonzero character sums to zero). Repeatedly bisect these equally sized classes until all labels are distinct, appending a check at every bisection. The result is the row space of a full extended-Hamming parity check, of rank five with minimum distance four. Its affine syndrome is consistent because the replay checks the assigned parity on every dependency; new rows may be assigned arbitrarily. No singleton is in that full dual, so shortening at a dummy costs exactly one half. Hence each listed guard bank has its full asserted population; the union is not a sum of individually optimized substitutions.

The numerical totals add `672 × 2,048` to the 60D base and `672 × 2,048 + 539 × 1,024` to the 61D base.

Run `python3 verify_extra.py PATH_TO_EXISTING_DIM62_67_RM2_SELFSPARSE_PACKAGE`. This replays the baseline independently, loads the compressed exact witnesses and checks the entire additional pairwise geometry (integer popcounts), every bank's local affine completion, and both totals. `extra60.json.gz.b64` and `extra61.json.gz.b64` are base64-wrapped gzip-compressed JSON with the original 64-bit support masks and one exact opposite-parity assignment per guarded pair. No optimization solver or floating-point tolerance is needed.
