# 65D: three weighted donors per two marked body sites and one rail

**Verified lower bound: 424,312,658**, +135,166 over the predecessor 424,177,492. This bound uses the same 65D four-site packing existence premise as the baseline. It is not a statement of optimality or a universal literature priority. Historical comparison: Edel–Rains–Sloane's 64D construction, 331,737,984 (https://arxiv.org/abs/math/0207291), plus two perpendicular points, 331,737,986.

## Immutable predecessor and source identities

Parent repository snapshot: `934d082034d5bc6c64a5ec5d7b5e3f8576c29962`. Use its verified `verifications/improved/dim65-68-dg-diagonal-triples/` 65D row. Published complete-atlas head-selector Git blob: `fe8ba58f9864d1e57ab07e80082f642e12e8a181`. Its 65D triple profiles are exactly those in its pinned sibling `dim65-dg-three-site-donors/verify_gates.py`. Git blobs: actual 21-group site and cubic assignment with 46 accepted hyperplane IDs `b59d97220036ed6900db59c9837349212e4ff880`; literal 65D profile verifier `9e6532e37ebdf15dbbbf1aa929e76d3dbcb60ad7`; 80 mother cubics `5d374d5e49ab72f28d580008868330a47dcfa390`. The predecessor's `verify.py` separately establishes the binary RM80 carrier, 40,414 signed sparse banks with even total parity, all 46 atlas fibres, and the DG minimum distance 24. The paired group uses the same exact DG generator, reusing the same 1024 quadratic parts and all 128 affine words inside a listed cubic; it inserts no unpublished binary-code premise. The old 65D profile numerator norms reproduced with denominator 1000 are 95,991,025, 95,993,300 and 95,993,300.

## Construction and angular certificates

Keep all old sparse signed bank vectors and both unit tail poles, and keep the old profiles of the twenty triple groups other than (cubic label 12; sites 21,36,42). Revert that group to unweighted owners initially, so the 335,544,320-word mother dense carrier is the common starting stock. Pair (21,42) uses previously unused cubic index 8 and rail sign −1; pair (7,36) uses previously unused cubic index 13 and rail sign +1. Index 7 was the unique unused ordinary coordinate. Cubic *indices* 8 and 13 index the literal palette, not a binary coefficient value. The 20 triples and these two pairs are disjoint and partition the 64 truth-table positions. Every one of the 22 groups is assigned a distinct mother cubic. Every group retires its 131,072 designated DG owners and replaces each by 3 normalized weighted vectors. The old triples remain exactly their old signed profile families. For each new pair the three multipliers and physical tail coordinate, divided by 1000, are:
```
(-5679, 1250, -σ420)
( 1250,-5679, -σ420)
( 2780, 2780,  σ4290)
```
Multipliers multiply the corresponding owner's body signs; the physical tail does not depend on its signs. Other 62 coordinates retain their ±1 owner's signs. AB have squared norm 95,989,941/10^6, C has 95,860,900/10^6. This module uses both signs of the high rail because disjoint C profiles with identical rails would collide. The opposite rails ensure compatibility; this is not an unproved all-pair replication.

For a profile p let Δ(p)=Σ_{marks}|p_j−1| and N(p)=64−r+Σ_{marks}p_j²+h_p². An unweighted retained dense word is a different mother word (designated owners removed), so raw dot ≤32+Δ. Two disjoint groups are distinct cubics, with raw bound 32+Δ(p)+Δ(q)+h_ph_q. Different owners within one group have unweighted dot≤16 by DG minimum 24, giving 16+Σ ε_j(p_j q_j−1)+h_ph_q (all 2^r possible signs checked). Same owner uses exact 64−r+Σp_jq_j+h_ph_q. With a 16-body sparse sign bank old parity is even for both the mother sign and sparse coset. When the one flipped site of a profile is present, one disagreement at magnitude ≥1 lowers the magnitude bound by 2. With C no parity is needed, since the largest magnitude sum of 16 body coordinates is 14+2×2.78=19.56. Two rail poles are legal since |h_p|≤√N(p)/2. The two old ordinary poles at the now marked site 7 are removed. See below for the retained 30 of the predecessor's 92 signed atlas fibres.

Use the same Bao–Ji four-site existence theorem (https://arxiv.org/html/1401.2022#S4, Theorem 4.9) from the parent: 10,741 supports on 65 sites, any triple in at most one support. A four-site vector has 4 signed unit entries. The only gate not certified by the universal magnitude bound is a C profile compared to a four-site support containing both of its marked body positions **and** the tail coordinate. Choose the labeling of the independently selected four-site packing so that both triples are uncovered. In any such packing with `b` blocks, the leave contains `C(65,3)−4b` triples. The leave link at `v` has a vertex for each of the other 64 coordinates and an edge `ij` whenever `vij` is uncovered. For any pair `vi`, the 63 triples through it are covered in pairs by blocks, so every link vertex has positive odd degree. Some link has at most `floor(3(C(65,3)−4b)/65)` edges, hence at most 33 when `b=10,741`. In a 64-vertex graph with no isolates and `m` edges, a maximal matching of size `s` leaves `64−2s` nonadjacent unmatched vertices, each incident to a distinct nonmatching edge. Thus `m ≥ s+(64−2s) = 64−s`. In our link this yields at least 31 disjoint edges. Relabel the link center as the tail and any two matching edges as the disjoint marked pairs `(21,42)` and `(7,36)`, relabeling only this four-site layer. Both potentially failing triples are absent, so all 10,741 supports and all signs are kept.

This bound uses existence of the Bao–Ji packing, not a new literal list of 10,741 supports. The same constructive relabeling rule applies to the published grouped 10,726-support design, whose sparse leave link has at most 35 edges. Using it yields the separate literal-input count 424,312,418.

The companion Node checker tests every required cross-product using exact BigInt: positive raw dot d is admitted iff 4d²≤N(p)N(q). This is an incremental proof against the pinned established carrier; the initial code verification itself is in its parent package, not asserted to be recomputed by this small checker. Optional argument: the immutable predecessor worktree's root path, to verify actual source identities and claim.

## Keep 30 predecessor atlas fibres

The 65D predecessor has one fibre at each of the two heads ±2√2 along the tail axis on each selected 32-site affine hyperplane. A fibre means the exact 2,048 signed body supports and its existing head, unchanged. Every 32-entry body coordinate is ±√2. For a fixed old hyperplane, both old signed nearest-head families have degree-≥3 selector fixed by the literal atlas; both the remaining atlas phase and the projected quadratic part of any new DG owner are quadratic in the hyperplane's five affine variables. Hence all binary sign differences lie in the 65,536-word coset obtained by xor of the literal selector, the incoming projected mother cubic, and RM(2,5). This exact overapproximation includes both head-specific phase cosets without needing to know each quadratic offset.

Check all these coset words for each of the two new cubics and each of the 46 inherited hyperplanes. If the Hamming weight is d, let ε_j be 0 for a mark absent from the hyperplane, +1 for agreement and −1 for disagreement. The raw dot in units √2 is
`32−2d + Σ ε_j(p_j−1) + 2σ_head h_p`.
A squared integer comparison with numerator denominator 1000 accepts positive raw v iff v² ≤9 N_num (the other vector has norm squared 72). For d≥8 apply the uniform cap 16+Δ(p)+2|h_p|; the literal check evaluates all lower weights. The following 30 pairs of zero-based hyperplane ID and head sign pass **also against the opposite new group**, hence they can be retained unchanged:
```
4 -1
5 1
16 -1
19 -1
33 -1
34 -1
38 1
39 -1
44 -1
44 1
53 -1
53 1
54 1
60 1
70 -1
75 1
77 -1
79 -1
89 -1
89 1
91 1
94 -1
104 -1
105 1
107 -1
107 1
112 1
114 -1
114 1
119 -1
```
All other old pairs remain compatible by inheritance. The integer verifier tests angle geometry and recomputes these 30 fibres from all raw quadratic differences, and optionally validates source Git blobs.

## Count

Old base without all weighted donor deltas/atlas/ordinary poles, including both tail poles and the 10,741 signed four-site supports:
`335,544,320 + 40,414×2,048 + 10,741×16 + 2 = 418,484,050`.
The old bound adds `21×2×131,072 + 46×2×2,048 + 2`.
The new bound adds `22×2×131,072 + 30×2,048`.
Equivalently retain 20 old groups, replace one old net +262,144 group by two net +262,144 groups, pay the 62 missing atlas fibres (126,976 points), and two ordinary poles. No four-site points are deleted. The exact change is +135,166.

## Replay

Run `node verify.mjs` for the incremental angle, atlas, and count checks. Optionally pass the predecessor repository root as a single argument to check the published input blob identities. The full predecessor verification and the Bao–Ji existence premise are separate, as described above.
