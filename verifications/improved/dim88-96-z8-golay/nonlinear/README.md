# Nonlinear Nordstrom–Robinson sign fibers for Z8-Golay sections

## Verified bounds

Against published commit `bede564fdf5ee8a07e79715736f1fe9fa4253750`:

| Dimension | Candidate | Increase |
|---|---:|---:|
| 88 | 6,967,552,256 | 25,315,328 |
| 89 | 9,270,596,604 | 27,257,856 |
| 90 | 12,341,351,672 | 24,709,396 |
| 91 | 16,435,507,572 | 19,134,548 |
| 92 | 21,896,439,920 | 12,449,792 |
| 93 | 31,056,799,856 | 12,272,924 |
| 94 | 47,703,347,340 | 4,408,320 |

These are explicit verified constructions. They use a nonlinear replacement for some restricted sparse sign fibers, while retaining the same Z8 parent and exact same physical supports. No new Nordstrom–Robinson code or global optimality is claimed.

## A 4,096-word binary length-24 fiber with four opposite pairs

`nr_partition.py` constructs the classical Nordstrom–Robinson code N by Gray mapping a displayed four-generator Z4 octacode. It explicitly verifies that N has 256 distinct words of length 16 and minimum distance6.

The linear binary span H has 2,048 words and minimum distance4. Eight stored translates N_i partition H; all words and translations are in `nr_partition.json`. The partition, rather than any unverified existence statement, is checked explicitly.

For every binary label x of length 4, form eight head bits as four opposite pairs: pair j is (x_j,1-x_j). Assign color c(x) equal to its low three bits XOR111 if its fourth bit is1. Thus equal colors identify only x and its bitwise complement.

Append every sixteen-bit word in N_c(x). The result B has 16*256=4,096 distinct length 24 words; every prescribed head pair is opposite.

Distance proof:
- Same label: distance at least6 inside the corresponding NR translate.
- Different labels, different colors: head distance is at least2, while the distinct tails belong to the distance4 code H, totaling at least6.
- Different labels, same color: labels are complementary, so head distance alone is8.

`nr_pair4.json` is the complete explicit B. Every one of its8,386,560 unordered pairs is additionally checked by the discovery script's scalar population counts. The independent verifier reconstructs B and checks all pairs using integer matrix products. Its minimum distance is exactly 6.

## Why this improves restricted fibers

Four prescribed opposite pairs can retain all 4,096 nonlinear words, versus2,048 from the preceding adapted linear code. Three pair rows can also retain4,096 by placing the unused fourth opposite pair in unselected coordinates. Several mixed pair/four-cell patterns improve as well. The nonlinear fiber has fewer total words than the16,384-word linear code, but can have more words satisfying the required row equations.

For each occupancy pattern, the certificate chooses whichever tested family retains more words: the existing linear [24,14,6] code, or B after an explicit binary translation. The stored ordered groups specify a coordinate permutation. Binary translation, permutation and taking a subset all preserve the minimum distance lower bound.

## Complete kissing configuration

The current row sets and exact counts are in `verified.json`. For a given physical support:

1. Inspect selected-row occupancies. Odd occupancies give no zero-row points; mixed odd/even occupancies give no equal-row points.
2. Choose the pattern indexed by `(number of two-cell rows, number of four-cell rows)` for zero sections, or the five counts of occupancies0 through 4 for equal sections.
3. Select its linear or nonlinear family. For an NR pattern XOR its explicit24-bit translation into every B word.
4. Order selected physical rows by occupancy and then row index. Zero sections omit empty rows; equal sections include them. Map the stored sign-position groups to the occupied physical cells in each selected row, in cell order. Map remaining physical support cells to remaining sign positions, both in increasing order.
5. Retain precisely those family words whose selected sums are zero or equal, as appropriate. Place signs times2 on the support.

Every map is a bijection of 24 cells, explicitly reconstructed by `verify.py`. Each support may use a different family, translation and permutation.

Within one support, binary distance>=6 gives product<=4*(24-2*6)=48. Across different supports, the unchanged complete support proof gives intersection<=12, hence product<=4*12=48 regardless of signs or family. Dense/sparse products are<=24*2=48; sparse/RS products are<=6*2*4=48. The root/pole bounds likewise apply to arbitrary signs. Thus nonlinear substitution creates no new unchecked cross-family condition.

All vectors have squared norm96. Dense Z8 layers, transversals when applicable, Reed–Solomon six-support layers, roots and poles are unchanged apart from moving to the explicitly selected row section. Dense projections and all remaining layer counts are rechecked. The section has an explicit identity minor and an isometric identification with R^d. Distinct support sizes separate layers, and bijections preserve distinctness within a support.

## Production replay

Run `python ../verify.py` from this directory, or `python verify.py` from the parent package. The parent verifier checks the complete recursive manifest and replays both inherited Z8/support proofs, earlier section proofs, adapted-sign proofs, both new nonlinear proofs and all semantic corruption controls. No solver or search is used.

- `data/nr_partition.json` and `data/nr_pair4.json`: explicit nonlinear source and partition.
- `data/nonlinear_patterns.json` and `data/equal_nonlinear_patterns.json`: family, translation and group witnesses.
- `data/final_zero_rows.txt` and `data/equal_results.txt`: explicit row masks and sparse counts.
- `data/verified.json`: complete final claims, pattern histograms and family contributions.
- `source/replay.py`: independent integer-matrix pair checks, source reconstruction, pattern counts and physical maps.
- `source/independent.py`: separate scalar-popcount pair proof, exact span/partition checks, independent pattern counts and physical maps.
- `source/guards.py` and `source/negative_controls.py`: structural production guards and candidate-specific semantic rejection controls.
- `AUDIT.md`: archived independent research audit. Its research paths and pre-integration packaging notes refer to the original candidate; production integration supplies the pinned witnesses and controls.

Both nonlinear scripts independently reconstruct the source and compute counts. They share only redundant structural guards, not distance/count arithmetic. Exact input supports and Z8 generator are read from the parent `data/` directory. The full production replay freshly rechecks all 338,559,231 parent support pairs in each of two representations. Sign-family replacement changes none of those supports.
