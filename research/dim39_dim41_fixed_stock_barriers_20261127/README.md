# Lower63, 27 November 2026 UTC: exact fixed-stock gates in 39D and 41D

No new bound. At live base `167bdc80ec8736435f0adf6ae802798f83d2adb2`, our verified lower bounds remain **783,468 (39D)** and **1,460,066 (41D)**. All arguments below are scoped to stated geometric stocks.

## 39D: the four-mark axial route is impossible on the fixed binary stock

Source parent: Yang et al., [arXiv:2609.35051v1](https://arxiv.org/html/2609.35051v1), Appendix B Tables 11–12, encoded in pinned `verifications/improved/dim39-qr-three-coordinate-heads/parent39.json`. The ten even cosets of a rank-15 binary kernel give 327,680 points at ambient magnitudes 2. A rank-four mark leaves 2,048 donors per valid extra coset. With 438,272 complete signed octads, 2,908 roots, and 64 stars, eight extra donors *would* produce 785,308, but none of the proposed axial heads works.

Quotient all 39 coordinate flips by the binary kernel. Exact multi-source BFS on the 24-bit syndrome graph gives distance histogram 10, 390, 7,410, 91,390, 822,510, 4,292,719, 7,518,920, 4,004,073, 39,758, 36 (distances 0 to 9). Even parity is a well-defined character of the quotient, every source is even, and flips toggle parity. Therefore an even completion of any candidate's 35 body signs has a parent at full Hamming distance at most 8. With axial amplitude 17/5 on one of four marks, choose free marks to make that completion even. If the parent agrees on the axis, the projected distance is at most 8 (needs 9); if it disagrees, the projected distance is at most 7 (needs 8). The exact minimum squared scaled slack of the axial obstruction is 1,225.

More generally no head with four amplitudes 0≤a_i≤17/5 and at most one exceeding 2 can be inserted over this fixed full parent. Let m=min a_i. Complete with each chosen head sign, and flip a mark of amplitude m if parity is odd. Spending a marked mismatch raises the forced body-dot lower bound by 8, so the parent dot is at least
`L=76+sum_i min(2 a_i,8−2 a_i)−4 min(m,2)`.
The normalized half-angle cutoff is `T=sqrt(39(140+sum_i a_i^2))`. On the sorted pieces of the amplitude box, L−T is concave; the only vertices have sorted amplitudes in {0,2,17/5}, with at most one 17/5. The nine positive differences `25(L^2−T^2)` are 7900, 19600, 1225, 32100, 13165, 45400, 25905, 24300, 5365; L is positive. The blocker works for arbitrary marked choices. Two enhanced coordinates are *not* ruled out by this inequality, and preserving ten old cosets is an explicit assumption.

## 39D: even the published head architecture needs structural changes

From pinned (783{,}468) package use `certificate.json`, `parent39.json`, `base_certificate.json`, and SHA-pinned `supports_3424.txt.gz`. Exact 4096-word shortened kernel has a pair at Hamming 12. Signature of a sparse row records its physical three-coordinate head and maximum body dot against each complete dense bank. The 344,064 dense points and 436,432/2908/64 sparse bank sizes reproduce. Three missing octad signatures: multiplicities 176/192/1472; physical heads (2,−2,−2)/(2,2,−2)/(2,0,−2); maximal body dots 20/20/24 against each extra bank. Blocking extra-bank sets are {10,11}, {8,9}, and {8,9,10,11}, zero-based.

Keep merely original dense banks 0,4,5,6,7 fixed; remaining heads can move independently. For extra bank j write normalized four-vector ((a,t)=(h,12)/sqrt(144+||h||²)). Same-bank actual body distance 12 forces (||h||²≤48), so (t≥sqrt(3/4)), (t≤1), each (|a_k|≤1/2), and its norm is 1. For each fixed bank with head h_i and minimum body Hamming distance d_ij, the required linear inequality is ((a·h_i)+(144−8d_ij)t/12≤sqrt(144+||h_i||²)/2). Required sparse row with head z and max body dot b yields (a·z+bt/12≤sqrt8). Round right sides upward after multiplying by 12×10^6, and replace lower t bound by 866025/10^6. The resulting bounded rational polytope contains any admissible normalized head. Enumerate all full-rank 4-facet intersections using integer Cramer and check all walls. By convexity its squared norm is at most maximum at vertices, and every blocking pair has exact upper bound <1. Integer millionths, rounded upward:
`(176,10):971867, (176,11):978160, (192,8):948209, (192,9):978160, (1472,8):948477, (1472,9):978432, (1472,10):971867, (1472,11):978432`.
Thus fixed 5 dense banks alone forbid any missing row even if each extra head is varied independently and all other sparse constraints are dropped. Cannot conclude no new octad/support choices, new parent cosets, or different reroutings of the five fixed banks. Do not repeat previously known Q8 independent-sign and star swaps.

## 41D: 41D+42D common body, one real head line, and exact contact obstruction

The two pinned packages `verifications/improved/dim41-43-affine-sections/configuration_{41,42}.npz` have literally identical 40×40 even integral Gram matrix G. Body x has (x^T Gx=3r), physical Gram is (G/3), radii r in {10,12,14,16,18}. A complete family is (residue c mod3, radius r), deduplicated; its exact population comes from the existing verifiers. Family counts 144,145,29,8,2. All 328 residues are in `Gc=0 mod3`. Because G is even, a difference in residue δ has integer squared norm congruent to (δ^TGδ mod18). Existing configurations are independently checked for emptiness of every nonzero difference norm below exactly computed thresholds. Take maximal available cut for each unordered residue ±δ over both packages; this yields 2,647 certified difference lower bounds (d=30:324, 36:1120, 42:1080, 48:122, 54:1). For unseen classes use congruence and positivity; zero difference class has min positive at least 54. The pair body bound (U=3(r+s)-d) is on twice the body product in the integral Gram. For two copies of the same family include zero distance, U=6r.

Model: use both antipodal axial poles, any chosen complete bodies from the common inventory, arbitrary real scalar heads `h=√2 t`, with the safe self-bound r+h²≤18. For each pair require the exact scalar inequality from the inherited ledger U:
`U+12tt'≤3sqrt((r+2t²)(s+2t'²))`.
Pole constraint is h²≤r/3. For a fixed body the angle with body axis is arctan(h/√r) and is confined to [−π/6,π/6], hence at most two copies, possible only at opposite pole-tight endpoints; self-bound prohibits them for r≥14. This is a statement about the ledger-certified scalar model, not that the ledger is sharp for every pair.

A radius-18 core (population 1,084,590) must be at height zero. The two core classes cannot coexist. Without a core the coarse upper is 1,344,308. Its possible radius-16 neighbours (population 116,235) have U=54 exactly and so are forced to h=±√2; 2 possible neighbors for one core, 6 for the other. Every pair of such neighbors has U≥48 and so cannot share a sign. At most two radius-16 banks. By reflection it suffices to enumerate unordered pairs, assigning heights t=+1,−1 respectively. With 0/1 large neighbors similar single-sided bounds suffice and are strictly below the incumbent. For every minor family r=10/12/14 consider rational cells t∈[z,z+1]/256 with the exact radius caps and bound the right side of each core/large-bank inequality from above on the cell. This yields a rigorous overcount (not a grid search of chosen heights). At r≤12 a double family requires both fixed endpoints t=±sqrt(r/6); its two bank conditions reduce to `192r−U²≥0` and `(192r−U²)²≥96r U²` for applicable positive U.

Four top skeletons gross 1,460,561 = incumbent+495; the next best gross totals are strictly smaller than the incumbent. In each leading skeleton, one required 495-point family conflicts with a 9180-point family for all remaining heights. On the two cells, bound `U+12tt'` below at a rectangle corner and its RHS above using maximal squared heights; test exact integer squares. Index families by first occurrence in configuration 41 then configuration 42, merging (r,c) only. Top skeletons, then forbidden minor pairs:
- core 122, big (0,150): (251,146);
- core 151, big (152,153): (8,161);
- core 151, big (154,157): (8,162);
- core 151, big (155,156): (8,161) and (8,162).
Thus the ledger-certified model cannot exceed 1,460,066; the original construction establishes equality inside it. No claim over arbitrary real 41D packings, head planes, new residues, higher self-radius caps, or replacing only part of the complete body families.

Native follow-up within the authorized slot verified that each decisive U=42 is attained by the actual complete families, not just the conservative ledger. The tested cross differences have 495 vectors of squared norm 36 each. Joining the exact 495-point r12 shells to these differences yields respectively 4,407, 4,638, and 4,407 contacts with the r14 shell. Identity: with (v=x-y) of norm 36 and (x) of norm 36, (y) has norm 42 iff (2x^TGv=30). Exact output includes a literal pair per gate. A way around this obstruction would require skipping conflicting subbanks (the contact counts are substantial), move into a second head dimension while adjusting the physical body projection, or new residue classes; rerolling scalar heights is not enough within the certified stock.

## Replay

From the directory containing these scripts, in a clone of this repository:

```sh
pin=167bdc80ec8736435f0adf6ae802798f83d2adb2
mkdir -p baseline39
for name in parent39.json certificate.json base_certificate.json supports_3424.txt.gz; do
  git show "$pin:verifications/improved/dim39-qr-three-coordinate-heads/$name" > "baseline39/$name"
done
for n in 41 42; do
  git show "$pin:verifications/improved/dim41-43-affine-sections/configuration_$n.npz" > "configuration_$n.npz"
done
node verify39.mjs
python3 parse41.py > parsed41-42.json
node verify41_model.mjs
python3 verify_contacts.py
```

The two Python scripts require NumPy. The contact replay checks the three literal attainable dot-product witnesses, not the full shell-join multiplicities. The latter were counted by joining the exact residue shells; existence of one witness per disputed pair is enough for the obstruction. No lower-bound table changes.
