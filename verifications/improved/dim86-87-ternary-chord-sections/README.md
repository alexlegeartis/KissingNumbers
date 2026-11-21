# Ternary chord sections: 7,422,120,168 in 86D and 7,422,120,170 in 87D

| Dimension | Previous repository lower bound | Verified lower bound | Gain |
|---:|---:|---:|---:|
| 86 | 6,228,511,388 | **7,422,120,168** | **1,193,608,780 (+19.16%)** |
| 87 | 7,266,268,830 | **7,422,120,170** | **155,851,340 (+2.14%)** |

These are spherical kissing codes, not claims of optimality or unrestricted worldwide priority. The 86D result is a new construction, not an extension of the previous 86D code. The 87D result embeds it and adds two orthogonal poles.

## Dense construction

Start with the [published extremal ternary quadratic-residue code](https://arxiv.org/abs/1109.1681) with parameters [48,24,15]. The exact generator is the [previously replayed QR48 matrix](../dim68-91-lattice-code-continuation/ternary_qr48_20261117/certificate.json). Fix its first symbol to zero and restrict its next eight symbols to the two nonzero letters. The projection to these first nine coordinates has rank nine, so there are exactly `2^8 × 3^15` admissible words.

Replace each letter by one of three vertices of a unit equilateral triangle. Delete the wholly fixed block, which never contributes to pair distances. In each of the eight two-letter blocks subtract the midpoint of its two available vertices. The translated vertices are the two opposite points of a line, each with squared norm 3/4. The remaining 39 blocks are unrestricted unit triangle blocks. This is an affine section of the ternary simplex embedding: it preserves the squared distance `3d` for each pair at Hamming distance `d`, but lowers the common squared norm to `39+8×3/4=45` and the ambient dimension to `2×39+8=86`.

Self-orthogonality over GF(3) makes every Hamming distance divisible by three. Distinct words have Hamming distance at least 15; the fixed coordinate makes it at most 47, hence at most 45. The normalized body products before dividing by 45 are therefore `45−3d/2 ∈ [−45/2,45/2]`. We may include the negatives of all points. No negative is already present: every full triangle block uses an alphabet disjoint from its own negative. This supplies **7,346,640,384** points before auxiliary layers.

## Sparse and pure-head layers

Take the existing [certified 42-block Pless word list and nonlinear phase code](../dim90-92-hexagonal-phase-layers/hexagon_sparse15_phase7812_20261117/PROOF.md). Filter that word list to zero in its local positions 7, 9, 35 (zero-based) and delete those columns, leaving 9,662 weight-15 words in 39 blocks. Use those blocks as the full triangle blocks of the new dense construction, and put zero in all eight chord coordinates. For each outer word use the 7,812 phases of length 15. In an occupied block put one of the three A2 roots that are pairwise at angle 120°, with the sign determined by the outer letter. Each root has squared norm 3, and its absolute dot product with any unit triangle vertex is at most 3/2. Hence the sparse/dense cross bound is `15×3/2=45/2`, including both signs of the dense code. Each differently marked outer position costs squared distance at least 3, and each differently phased position of the same outer word costs 9; the outer and phase minimum distances 15 and 5 give the needed squared distance 45. This supplies **75,479,544** points, distinct from the dense layer because of their zero full blocks.

The eight chord coordinates have combined squared norm 6 in each dense vector and zero in every sparse vector. Place the **240 E8 roots**, scaled to squared norm 45, solely in this eight-dimensional space. Their product with any dense vector is at most `sqrt(45×6)<45/2` by Cauchy–Schwarz. They are orthogonal to the sparse points and disjoint from both other layers. In 87D add the two antipodes on a new coordinate axis.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim86-87-ternary-chord-sections/verify.py
```

The script pins all inherited inputs by SHA-256, independently checks the ranks over GF(3), verifies every possible dense distance and every E8 pair, and checks all 46,672,291 outer and 30,509,766 phase pairs. The standard QR48 minimum-distance theorem remains the explicit inherited input; this verifier does not enumerate the 3^24 parent. The sparse-layer verification does not need to assume a minimum distance for its original longer outer code because it checks the selected 39-block list directly.
