# Independent nonlinear sign-fiber audit

## Verdict

All seven supplied candidates PASS independently against current main. No mathematical blocker was found.

|Dimension|Verified total|Current published baseline|Gain|
|---|---:|---:|---:|
|88|6,967,552,256|6,942,236,928|25,315,328|
|89|9,270,596,604|9,243,338,748|27,257,856|
|90|12,341,351,672|12,316,642,276|24,709,396|
|91|16,435,507,572|16,416,373,024|19,134,548|
|92|21,896,439,920|21,883,990,128|12,449,792|
|93|31,056,799,856|31,044,526,932|12,272,924|
|94|47,703,347,340|47,698,939,020|4,408,320|

Baseline source: freshly fetched remote RESULTS.md, Git blob `6ea1e0cdad095ac9a10f2253ae92ad771b28ef26`; exact relevant rows are in `published_baseline.json`. Every candidate baseline was programmatically compared to these rows.

Reproduction: run `OPENBLAS_NUM_THREADS=1 python research/z8_nonlinear_signs_audit_20261110/audit.py`. This imports no producer proof functions and writes `report.json`; `audit.log` records the completed run. Input hashes are pinned in the receipt.

## Actual NR construction and partition

The audit independently enumerates all 4^4 coefficient tuples of the supplied 4x8 octacode generator, using scalar modular arithmetic and the Gray map [0,1,3,2]. It obtains exactly the 256 supplied distinct 16-bit words and scalar-counts all their pair distances, verifying minimum distance six.

It rebuilds all eight supplied translates, verifies equality with each stored partition component, and checks that their union has 2,048 distinct words. It independently generates the span from the supplied binary basis and verifies that it is exactly this union, contains the original NR code, and has minimum nonzero weight four. Thus the partition is not merely assumed to be a partition of a distance-four code.

For each 4-bit label x, the construction puts opposite bits into each of four head pairs, then appends the NR translate with color (x&7) XOR (7 if x has its high bit set else 0). The audit reconstructs all 4,096 distinct 24-bit words and compares the complete set with nr_pair4.json.

All 8,386,560 unordered pairs were then checked by scalar XOR/popcount, independently of the producer's integer-matrix pair proof. Minimum distance is six. Exact distance histogram:

- 6: 393,216
- 8: 702,464
- 10: 1,966,080
- 12: 2,260,992
- 14: 1,966,080
- 16: 702,464
- 18: 393,216
- 24: 2,048

An additional structural explanation confirms the result: equal labels have NR tail distance at least six; distinct labels with the same color are antipodal labels, giving head distance eight; different colors lie in disjoint components of the distance-four span, so tail distance at least four combines with head distance at least two.

This verifies the supplied finite code without relying on the Nordstrom–Robinson name or a literature assertion about it.

## Fibers, translations and physical assignments

Every support receives exactly one family: either the 16,384-word linear distance-six code or a translated copy of the 4,096-word nonlinear code. It does NOT receive their union. This distinction is essential: compatibility of overlapping fibers on a single support would require another argument.

For every stored zero/equal pattern, the audit verifies family tag, translation range, group sizes, disjoint sign positions and exact count. It independently enumerates translated scalar binary words and computes each group sum by popcount. A binary translation preserves all pair distances. Linear patterns use zero translation. Mixed-parity equal-row occupancies necessarily contribute zero.

Every retained support's physical assignment is instantiated: selected rows are ordered by increasing occupied-cell count, breaking ties by physical row index; their occupied cells are ascending. Stored sign groups are matched in order. Unused sign positions and remaining occupied physical cells are matched ascending. All mappings are checked bijective and their pullbacks are checked to equal the stored pattern groups, including empty groups for equal sections. Thus the row-sum conditions are conditions on actual physical coordinates, not merely abstract sign positions.

The pattern histograms, number of physical permutations and separate linear/NR point totals match every supplied record. One family is selected per support, with all accepted sign words distinct.

## Dense counts, ranks and compatibility

For every new row set, the projected Z4 image (zero sections) or Z8 image (equal sections) is independently enumerated with array group closure, rather than using the producer's tuple-set enumeration or a presumed initial-column projection. Exact uniform fibers reproduce every dense count. Both Gray cosets are tested for equal sections, and transversal sign counts are independently enumerated from the original linear code.

The section constraints have explicit identity minors, proving rank 96-d. For zero sections, A3 roots and outside coordinate poles span the full section. For equal sections, the all-+1 dense point adds the common-row-sum direction to that zero-row-sum span. The resulting configurations therefore span exactly their claimed dimensions. An orthogonal identification with R^d preserves norms, inner products and distinctness.

All vectors have squared norm 96. On one support, distance at least six gives product <=4*(24-12)=48, whether the selected family is linear or nonlinear. On different supports, overlap at most twelve gives product <=48 regardless of family, translations or coordinate permutations. Additional/transversal supports also overlap at most twelve by the previously verified trace obstructions. Dense/sparse products are at most 48, sparse/six-support products at most 48, and the pure-root/pole cross-bounds remain strictly below 48. The six-support layer uses only outside cells, with sufficient space checked for each row set.

The underlying support bytes match both earlier fresh full parent replays exactly (SHA-256 `d0def7d6587ca3aa1b39c9f3e7ccb844310285a7f1794dc8e9492f0fff9fc65d`). Those already checked all 338,559,231 support pairs, the trace obstructions, the classical auxiliary codes, and the Z8 parent distance in two representations. This follow-up reuses those validated unchanged inputs; it does not claim another full parent replay.

Layer support sizes separate dense, weight-24, weight-six, roots and poles. Transversal supports occupy 24 rows and additional supports at most 14, separating the two weight-24 families. Support uniqueness, code uniqueness, translation and bijective placement prevent duplicates within each component.

## Remaining gaps and scope

No proof gap affecting the seven numerical bounds was found. Production integration should pin the NR partition, complete word set, pattern tables, physical-map convention, section row sets and parent inputs, and add candidate-specific corruption controls for all of them. The standalone producer verify.py does not itself replay all parent geometry; its relationship to the parent verification must remain explicit in the final wrapper. The independently rebuilt partition/span check here is stronger than trusting the names or stored receipts.

Pattern and row searches are heuristic or restricted to the specified construction class. No global optimality or worldwide-priority conclusion follows. The classical nonlinear parent and NR code are prior ingredients, not new code discoveries. This audit made no remote changes and does not itself authorize publication.
