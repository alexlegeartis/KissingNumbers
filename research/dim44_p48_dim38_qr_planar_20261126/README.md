# Lower57: P48 exact ceiling and planar QR in 38D

Nominal run 22:30–23:00 UTC, 26 Nov 2026. Main 26998bee3ff7338c043a885913d6e9a9cfacf4a0: 44D 4,940,961; 38D 652,076. External named comparisons: Conway–Sloane / Dorofeev–Sun–Wang 4,153,600 at 44D (https://arxiv.org/abs/2607.20359v4); Lin 605,322 at 38D (https://arxiv.org/abs/2611.03682). Yang et al. gave the QR39 parent (https://arxiv.org/html/2609.35051v1); Lin gave the unchanged octad code.

## 44D all-calibration obstruction
Incumbent already uses norm6/8 and a normal torus, so neither is new. The stronger proof uses the shortest-distance bound with absence of shorter source-shell differences: in a difference class with shortest normal norm q, use projected squared distance ≥6−q, or ≥8−q if its 6-shell bin is empty, or ≥10−q if both 6/8 are empty. The 324 scaled √3 D4 quotient classes have minimum representative norm≤3 by exact enumeration, so all complete fibres with body radius squared≤6 arise from shells6/8. The source 328 affine population identities are exact from the pinned 44D moments.py proof. The actual frame has 404 positive fibres; 3 empty short bins add 507 unordered compatibility edges, still exact cap 4,940,961.

For every possible source calibration, the 328 nonnegative affine count forms imply a 9-facet polytope with 12 integral vertices. Bounds c_i≥0 and sum c_i≤276 immediately eliminate the other inequalities. For a fixed clique, required absent bins are zeros of nonnegative linear forms and define a face. Counts are affine on that face, so a vertex of the full polytope attains at least the same weight; drop vertices of zero weight. Exact weighted clique exhaustion at all 12 vertices proves the same 4,940,961 ceiling. It is not an upper bound on actual endpoint-compatible or arbitrary 44D codes. Node stored-source readback PASS: `p48-d4-ceiling.mjs.`


## 38D planar aperture
Mark source 39-bit positions [0,1,31]. Replace signs on these sites by orthogonal projection of their ±2 triple onto the sum-zero plane, leaving 36 unchanged ±2 coordinates. The six nonconstant triples give a regular hexagon of squared radius32/3. Its equalities match the two tight original Hamming10 cases. The central coset provides the two constant prefix fibres as well: give 000 the head of singleton-mark31, and 111 its antipode. Weight10 differences from the center coset to all others avoid code coordinate31; independent exact 1891 checks prove the 62 fibres share a shortened [36,12,d12] kernel, 4096 each. Zero-head placement alone cannot improve on the 60 bright fibres in any of the 9139 markings by the exact center-neighborhood unions; sending the constant prefix on opposite hexagonal vertices is the useful move.

Every retained 38D signed octad/root/star survives uniformly: in the norm8 sparse scale max octad raw dot with dense body/head is 16, 14+r, or 12+√2 r for 0/1/2 plane coordinates, each below √(288+2r²) for r≤√(32/3). Root/star estimates are below15. The exact numerical tie is 62×4096+398124=652076; an additional whole fibre gives 656172. Suite: `qr38-hex-certificate.mjs`; exact discovery `qr38-hex-search.cpp` and `qr38-hex-input.txt.`

The discovery algorithms compute all 2^24 Hamming quotient distances by BFS on 36 syndrome generators. Arbitrary head in the safe disk gives outer linear intervals in z=3 head coordinates: z_j≥50−6d_i or z_j≤6d_i−50, from replacing variable norm by cap; sum z=0 and ||z||²≤96. This is an exact necessary condition for a new whole fibre keeping all 62. Separately check 61 finite half-integer head signatures and joint exchanges of 1 or 2 old slots. Independently check every positive witness directly by 4096 symbolic codeword differences. 

## Exact replay outcomes

The 44D clique searches finish at all twelve vertices and replay the limit 4,940,961. The independent 38D verifier reconstructs the 652,076-point tie. In the original safe disk the exhaustive syndrome scan finds no whole-fibre append, 324 candidates that each require an old fibre to be removed, and no candidate requiring exactly two removals. With verified mutual head constraints on the finite menu, none of 460 outgoing pairs gives a 1→2 or 2→3 gain; no pair was truncated.

The separate continuous append verifier increases the head radius to 93/20 with the same 62 full incumbent fibres. It checks the complete Hamming quotient of 16,777,216 syndromes. Of 880 coarse candidates, 18 satisfy the outer inequalities, all with one of two opposite six-hue patterns. In scaled plane coordinates z with sum zero, those patterns require z_i≤b−2 for all i or the sign-reversed constraints, where b=(sqrt(1392(1296+||z||²))/2−696)/12. Then b≥2 and the radius bound gives b≤25/12, so every |z_i|≤1/6. But the defining squared norm is 96+48b+(12/29)b²≥192. This excludes all 18. It does not cover a joint continuous exchange of multiple fibres or arbitrary 38D codes.

Independent commands (Node 24 and C++17):

- `node p48-d4-ceiling.mjs`
- `node qr38-hex-certificate.mjs`
- `g++ -std=c++17 -O2 qr38-hex-search.cpp -o qr38 && ./qr38 < qr38-hex-input.txt`
- `node qr38-hex-wide.mjs`

All four exact independent replays passed. No record table changed.
