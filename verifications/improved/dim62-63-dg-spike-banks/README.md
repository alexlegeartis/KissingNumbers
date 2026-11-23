# Delsarte–Goethals spikes and shorter shared guards in 62–63 dimensions

| Dimension | Previous verified lower bound | Added contacts | New verified lower bound |
|---:|---:|---:|---:|
| 62 | 172,134,412 | 600,576 | **172,734,988** |
| 63 | 230,159,620 | 2,031,616 | **232,191,236** |

These extend the [62D two-mark packing](../dim62-two-mark-diagonals/) and the [63D guarded equality section](../dim62-63-guarded-rephase/) with a new unequal-amplitude dense layer. The old dense, plane, pole, and sparse support coordinates and every old population remain; some old sparse sign cosets are rephased jointly. The 62D result also adds 213 sparse two-mark banks via one stronger parity guard rule.

## Dense donors in an omitted prefix

For each ordinary coordinate, assign a different one of the 80 existing cubic RM80 cosets. In that coset take the translate of a length-64 quadratic Delsarte–Goethals code, restricted to one designated omitted prefix. Within each cubic there are 1,024 chosen quadratic parts and all 128 affine parts. The checker constructs the parts by trace bilinear forms over GF(32) and extension to alternating forms on six variables. It checks all 15,903 distinct nonzero differences of quadratic parts against all affine words and proves binary minimum 24. Cubic translates preserve this within-site distance. Affine evaluation has rank four on the 62D marked sites (0,9,20,41) and rank two on the 63D marked sites (0,20), giving exactly eight and 32 choices per quadratic part, respectively.

- **62D:** Use hexagon prefix 5 (bits 1,0,1,0 in the ordered marks), which is at distance at most two from every old dense prefix. Retain the 60 ordinary signs, but raise the absolute value on the assigned site to 22/5, and put zero in both tail coordinates. This gives squared norm 78.36 and 60 × 8,192 = 491,520 points.
- **63D:** Use the omitted equality-section prefix (0,1), take 62 ordinary signs, raise the assigned absolute value to 4, and put zero on the tail. Squared norm is 77 and the population is 62 × 32,768 = 2,031,616.

At the same assigned site distinct donors have at least 24 ordinary sign differences. At distinct sites they lie in distinct cubic cosets at binary distance at least 16, with no distance lost on the common prefix. Against old dense points, delete at most two differences in 62D or exactly one at the prefix in 63D (the total old/candidate parent distance is at least 16). Direct dot bounds after weighting are checked by squaring against the product of the squared norms, not by rounding normalized coordinates. In particular the largest possible donor/old-dense raw dots are 35.4 in 62D and 35 in 63D.

## Signs of the old sparse supports

Every sparse mask is orthogonal modulo two to every row generating the original dense words, which the replay checks. Rephase zero-mark banks to odd total parity, preserving all single- and double-mark parities required for the earlier dense construction. The chosen omitted prefix makes a donor's ordinary restriction parity opposite to every single-mark sparse bank's required ordinary parity, as well as to all zero-mark banks after rephasing. Thus a mixed pair with one of those banks has at least one ordinary sign mismatch, discounting the raw dot by at least two. The largest such dots after this discount are 17.4 in 62D and 17 in 63D. In 62D the two-mark banks require no donor parity discount because they have only fourteen ordinary entries. Old planes and poles meet the exact angular limit by their fixed norms and the possible spike magnitudes.

All old overlap constraints remain in force. The proof uses the original guards with fresh opposite parity assignments from the compressed hexadecimal lists. Each bank's signed guard dependencies and total parity are checked simultaneously. The homogeneous span contains the original full support and every other nonzero row has weight eight. Its equal-size column classes can be successively bisected to give the extended Hamming dual, of rank five and minimum primal distance four. Fixing one or two omitted dummy signs retains respectively half or one quarter of the sign words because the full dual has no row of weight at most two. All old populations therefore remain unchanged.

## Short guards in 62D

Retain 213 more adjacent double-mark supports from the actual weight-16 parent shell, with the same two directed unit entries as in the prior 62D packing. There are 223 nonautomatic pair checks involving these supports; their worst raw dot is at most ten. Each has opposite parity on a common **ordinary** subset, now possibly of size six or seven: append one or two missing dummy coordinates separately in each bank to lift the ordinary subset to a weight-eight row. The two banks need not use the same source dummy names. Since all such dummy signs equal zero, the guards still force an ordinary mismatch and reduce the dot by two. The full signed-span replay certifies the assigned guards and all 512 signings per support. The extra gain is 213 × 512 = 109,056.

## Replay

Run `python3 verify.py` in the proof directory (Python 3 and a C++17 compiler). The complete replay first checks SHA-pinned predecessors, then the exact DG difference enumeration, omitted-prefix counts, every rephased parity bank and all 62D short guards and pair comparisons. `claims.json` contains both final totals. No floating-point comparisons or optimizer output are used as premises.
