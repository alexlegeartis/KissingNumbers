# 60D: a common-class 48-frame and a 32-color tail

This gives **58,665,998** kissing points, an increase of **1,814**
over the preceding live bound of 58,664,184. It replaces the later
low/heavy auxiliary cubes over the [certified lifted
P48p body](../dim50-60-p48-global-owner-refill/); it is not a union with
the previous norm-12 auxiliary layer.

The literal P48p vectors in `witness.json` form an orthogonal 48-frame
of squared norm 12. Their congruence modulo 2 is the key property:
partition sixteen axes into eight pairs, take both half-sums and
half-differences, and obtain sixteen integral orthogonal norm-six
directions. Each unused frame axis is orthogonal to these controls.

The controls carry 4,096 unit-norm sign-cube points. For every even
sign word of length twelve take the control indexed by the bitwise XOR
of its flipped positions, with both signs on that control. Scale the
normalized control by 1/2 and the twelve signs by 1/4. The words of a
fixed XOR color have binary distance at least 4; across colors the
controls are orthogonal and the binary distance is at least 2.
The 512 retained existing tetrahedral sparse rays use sixteen other
head lines and, as in the inherited proof, cost their 32 equatorial
rays. We retain none of the former heavy cube modes.

The other 32 frame axes receive the 1,536 literal tail vectors Z of
squared norm 10. Each assigned axis w gives two new points
`(±w/(4√2), Z/4)`. Their pairwise test is
`|w·w'| + 2 Z·Z' ≤ 16` for distinct tails; opposite head signs at
a common tail are compatible. The orthogonal controls have no
horizontal cross product with this layer. The other cross products
are checked with the inherited coordinates or squared-radical gates.

For each of the 420 old antipodal tail classes the script checks
all *signed* candidate owners, including both signs in the earlier
representative arrays, and counts any higher-than-1/2 interaction
with the new paired layer. The 236 deletions are conservative:
the scan includes a superset of actual owners, so it cannot
undercharge already retired candidates. The retained count is
58,658,586 + 4,096 + (512−32) + 3,072 − 236 = **58,665,998**.

Run `python3 verify.py` with NumPy; the inherited files are pinned by
SHA-256. `python3 audit_angles.py` independently checks the literal
new sign words, all new cross geometry, and the signed deletion
threshold. Its matrix scan works with bounded exact integers
represented in binary64: every intermediate sum is strictly less
than 2^52, and the main verifier separately uses int64.
Neither script uses the optimizer that selected the frame or tail palette.
No upper bound is claimed.
