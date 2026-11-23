# 7,906 phase words and four larger QR–Pless kissing configurations

The attached length-15 ternary phase code has 7,906 distinct words and minimum Hamming distance 5. The independent verifier checks all 31,248,465 pairs, with exact popcounts and no search dependency. This improves the 7,812-word phase used by four live QR48–Pless constructions; the separate 7,929-word coding-theory bound reported by Brouwer's table from Oscar Sanchez is not claimed beaten, and no public explicit 7,929-word witness was found.

Replace the entire old phase code only for zero-axial-occupancy sparse outer words. All old configurations with positive axial occupancy are unchanged. These are replacements, not unions of old and new phase codes. Keep all dense and small-shell layers from the pinned live configuration.

| Dimension | Live predecessor | Number of zero-axial outer words | Increase | Certified total |
|---:|---:|---:|---:|---:|
| 86 | 8,008,552,282 | 9,446 | 887,924 | 8,009,440,206 |
| 87 | 8,160,988,416 | 9,440 | 887,360 | 8,161,875,776 |
| 88 | 13,358,363,692 | 37,392 | 3,514,848 | 13,361,878,540 |
| 90 | 21,944,111,386 | 130,944 | 12,308,736 | 21,956,420,122 |

Source: the live 86/87 [rotated-plane layer](https://github.com/alexlegeartis/KissingNumbers/tree/main/verifications/improved/dim86-87-90-rotated-plane-trades) inherits the zero-occupancy counts from its [remarked-square predecessor](https://github.com/alexlegeartis/KissingNumbers/blob/main/verifications/improved/dim86-87-remarked-square/claims.json). The [88D live certificate](https://github.com/alexlegeartis/KissingNumbers/tree/main/verifications/improved/dim88-full-pless-chord-marks) explicitly enumerates the stated occupancy. The 90D live rotated-plane layer retains the full [pinned 45-block weight15 Pless list](https://github.com/alexlegeartis/KissingNumbers/blob/main/verifications/improved/dim90-92-hexagonal-phase-layers/hexagon_sparse15_phase7812_20261117/PROOF.md). No other row is claimed.

## Mixed-shell audit (all phase choices covered)

The changed sparse points have fifteen ordinary nonzero blocks, zero axial coordinates, and squared norm 45. In each occupied block the unchanged sign of the outer Pless symbol selects either the triangle of three A2 roots of squared norm3 or its negative. Against a unit triangle vector (either dense sign), each local absolute dot is at most3/2, so the dense mixed dot is at most22.5.

If two sparse outer words differ, each differing block contributes squared distance≥3 regardless of phases: zero/nonzero contributes3 and opposite root orientations contribute either3 or12. The outer word separation≥15 is inherited and unchanged. For equal outer word, phase separation≥5 contributes 9 at each differing phase, a total≥45. This also covers cross comparisons with the unchanged positive-axial-occupancy sparse points: their outer words necessarily differ (and any axial differences contribute at least3).

Against an old square point of squared norm4, at every occupied plane its unit entry is at least15° from every signed A2 root. The mixed dot is at most4√3 cos(15°), squared equal to 24+12√3<45. For the old rotated five-plane points of squared norm5, the angular separation is at least30° in each plane; the mixed dot is at most 5√3 cos(30°)=15/2, the exact threshold. For any old three-plane point of squared norm3, the dot is ≤3√3<√135/2. Each rank-two-coordinate D root of squared norm2 sees at most two ordinary blocks of squared norm3, giving squared mixed dot≤12<45/2; the same bound is homogeneous for the scaled pure-head D tails. The retained small-shell points occupy at most five ordinary blocks, so cannot coincide with any replacement point. Dense vectors occupy all 39/42/45 ordinary blocks; different outer words or different phases give distinct replacement points because the six signed root vectors and zero are distinct.

Run python3 verify.py with NumPy and assertions enabled from this directory. It SHA-pins all inherited source claims, rebuilds each ternary word, checks every pair and recomputes the exact layer counts. The pinned phase file SHA-256 is 220a8d800eaaf391898b424e5e27f54fe957c4f0548d6ce980e9b488637ab369.

## Discovery

Recovered ten exact isometries of the original phase set from single-coordinate and pairwise trit counts, then used explicit sweeps of short outside words and debt-limited blocker trades (a new compatible set may displace fewer old codewords than its size). Finally combined different trades with automorphism-conjugate bipartite matching. The maximum independent set in the union of two codes reduces to maximum matching between their exclusive parts, since both exclusive parts are internally independent and shared words are compatible with all. Release verification is independent of this search. 
