# Independent equal-row adapted-sign audit

Both follow-up claims PASS:

- 93D: 31,044,526,932, gain 2,678,788 over current published 31,041,848,144.
- 94D: 47,698,939,020, gain 630,840 over current published 47,698,308,180.

Baselines agree with the freshly fetched RESULTS.md blob `5ee520eafa6c797f191b03a79bc9d948c541ac74` from the preceding audit pass. The independent checker is `audit_equal.py`; `equal_report.json` pins input hashes and records counts.

## What was independently checked

All 105 stored equal-row pattern records were checked using explicit binary generator matrix multiplication and signed integer matrix products, without importing source proof code. Mixed even/odd row occupancies necessarily have incompatible sum parity and contribute zero. For every other pattern, every group has the prescribed size, all sign positions are distinct and in range, and enumeration of all 16,384 sign words reproduces the count.

For each support, order selected physical rows by increasing occupied-cell count, then by physical row index. Assign stored groups in that order, matching each group's ordered sign positions to ascending occupied physical cells. Empty groups map to empty rows. Assign unused sign positions ascending to all remaining occupied physical cells ascending. The audit instantiates this mapping and checks it is a bijection and pulls the physical equality constraints back to the exact stored group constraints.

This reconstructs 25,018 physical permutations in 93D, omitting 1,004 mixed-parity supports, and 25,261 permutations in 94D, omitting 761 mixed-parity supports. Histograms match exactly. Sparse totals are 139,980,900 and 184,949,152 respectively.

The 93D row set is [0,1,8,16], not the preceding [0,1,4,8]. Uniformity cannot simply be inferred from the old first-twelve-column argument. Complete direct image enumeration verifies its Z8 image has 4,096 vectors and 1,810 qualifying images across both Gray cosets. The 94D row set [0,1,8] has 512 images and 346 qualifying images. Exact uniform fibers give dense counts 30,366,760,960 and 46,439,333,888. Direct sign-code enumeration gives transversal counts 536,870,912 and 1,073,741,824.

Both use the same 913,952-point RS layer, placed outside the selected rows; sufficient outside cells are verified. Pure-root counts are 208 and 204. All total and gain arithmetic is independently checked.

## Complete compatibility and dimension transfer

Within each support, a fixed coordinate bijection preserves binary distance at least six. Distinct additional supports still intersect in at most twelve cells, so arbitrary independently chosen signs give inner product at most 48. The support bytes are the same as those validated by both full parent replays in the preceding audit. Additional/transversal support intersections are also at most twelve by the previously replayed dual-trace obstructions, independently of the additional supports' sign permutations. Thus the equal-row transversal layer remains compatible.

All other bounds are unchanged: dense/weight-24 <=48, weight-24/weight-six <=48, root/weight-24 <=4*sqrt(48)<48, and pole/weight-24 <=2*sqrt(96)<48. Six-support vectors avoid selected rows. Norms are 96, and the unchanged parent and other layer checks cover every remaining pair category.

The equal-row constraint matrices have identity minors of size 3 and 2, so their kernels have dimensions 93 and 94. Roots and poles span the zero-row-sum subspace, while the all-+1 dense point supplies the common nonzero-row-sum direction, proving the actual configurations span their claimed dimensions. Layer support sizes, distinct support families, unique sign words and coordinate bijections prevent duplicates. Orthogonal identification of the section with Euclidean space preserves all norms, products and distinctness.

## Limits

No mathematical defect was found. The source checker still needs an explicit physical-map convention in final documentation; the convention above closes that presentation gap. Integration needs manifests and corruption controls, and should explicitly include additional/transversal compatibility when discussing equal-row cases. Heuristic pattern selection proves a construction count, not an optimum or worldwide priority. No remote files were modified.
