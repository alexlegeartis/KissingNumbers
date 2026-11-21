# Axial Pless layers in 86 and 87 dimensions

| Dimension | Antipodal dense QR48 section | Sparse Pless layer | Root tail | Lower bound | Increase |
|---|---:|---:|---:|---:|---:|
| 86 | 7,346,640,384 | 661,565,934 | 112 | **8,008,206,430** | +586,086,262 |
| 87 | 7,346,640,384 | 813,992,418 | 144 | **8,160,632,946** | +738,512,776 |

Run `python verifications/improved/dim86-87-axial-pless/verify.py` from the repository root with NumPy installed and assertions enabled. The replay uses the SHA-256-pinned QR48 generator from the preceding [chord construction](../dim86-87-ternary-chord-sections/) and the Git-blob-pinned [Pless data](../dim90-92-hexagonal-phase-layers/hexagon_sparse15_phase7812_20261117/). The two published parent ternary codes' minimum distance 15 is inherited; see the pinned base packages and [Nebe's discussion of the Pless code](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/DIM48.pdf).

## Construction

Use 39 mutually orthogonal triangle planes, eight chord axes, and in 87D a ninth axis. Fix QR48 coordinate 0 to zero and delete its constant triangle. Restrict coordinates 1, 2, 24, 28, 34, 39, 41 and 47 to symbols 1 and 2, translating each chord to its midpoint. The projection onto these nine marked columns has rank nine, so the word count is 2^8 × 3^15; take both antipodal images. This is a re-marked section, not the exact earlier dense point set. All points have squared norm 45. Distances between different words are 3d, with d a multiple of three between 15 and 45; hence their dot products lie between −45/2 and 45/2. A full triangle is disjoint from its negative, so the two images are disjoint. Leave the ninth axis zero for the dense layer.

For the sparse layer use weight-15 Pless words, with coordinate 0 forced zero in 86D. A nonzero outer trit is a sign: on each full plane place the corresponding signed choice of one of three A2 roots of squared norm 3, and on each axial slot place ±√3. Use the phase-code coordinates in increasing order of occupied full-plane indices. Every sparse point has squared norm 45. Against either dense sign, each local dot product has absolute value at most 3/2. Distinct outer words differ in at least 15 locations, each contributing at least 3 to squared distance. For the same outer word, different phase words differ in at least five locations, each contributing 9 to squared distance. The phase sizes by axial occupancy 0–9 are 7,812, 6,561, 2,187, 729, 729, 243, 81, 27, 9, 3. The first code is checked all-pairs; three sizes come from prefixes of a rank-six parity-check matrix and the rest from shortening a cyclic ternary Golay code.

The finite outer set is constructed from the pinned 130,944 weight-15 words by exact membership checks, cyclic shifts, and a signed column permutation checked on the full generator. The replay reconstructs 285,384 distinct words in 86D and 415,104 in 87D, then takes the exact weighted phase census. Completeness of the weight-15 enumeration is not needed.

Finally place all D8 or D9 roots on the axes, scaled to squared norm 45. A root sees at most two sparse entries, for an absolute dot product at most 2√(135/2) < 45/2. The dense tail has squared norm 6, so Cauchy–Schwarz gives √270 < 45/2. Occupancy of the full planes distinguishes the three layers. The earlier sparse and root layers are replaced, not added again.

- `verify.py`: replays the dimensions, actual marked QR48 projection, dense Gram bounds, and new layers.
- `axial_pless.py`: exact finite orbit, phase and root checks.
- `claims.json`: counts and comparison with the preceding construction.
