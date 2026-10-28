# Larger compressed-code bounds in 61D and 62D

The 28 October 2026 continuation gives two substantially larger constructions:

| Dimension | Previous repository bound | New bound | Gain |
|---|---:|---:|---:|
| 61 | 72,376,856 | **87,056,170** | **14,679,314 (20.28%)** |
| 62 | 100,437,548 | **120,620,044** | **20,182,496 (20.09%)** |

The comparison is with commit `452d929` (the same 61--62D values as `aa1b7c3`).
These replace the preceding constructions in their dimensions. They do not
append to the older bridge or section layers. Other dimensions retain their
existing certificates. The manuscript remains its earlier snapshot; no new
worldwide priority survey or unrestricted optimum is asserted.

## Common binary and sparse ingredients

Let C be the explicit binary [64,28,16] code rebuilt by
`../dim49plus/binary_62_code.py`. The verifier freshly enumerates all 2^28
words, checks minimum weight 16, and collects all 30,828 weight-16 supports.
It checks self-orthogonality directly on every pair of generator rows.
Consequently every word of C has even restriction to every support in C.
This parity fact enables the new 61D layer.

Remove the four coordinates F = (0,1,12,52), leaving a 60-coordinate head in
the original coordinate order. Exactly **9,673** of the weight-16 supports
avoid F. Their pairwise intersections are at most eight, by the minimum
distance of C. A separate packed-integer pass checks all **46,778,628** pairs.
An exhaustive search over all binomial(64,4) = **635,376** choices of F
reproduces the maximum 9,673; 1,008 choices attain it. This is an optimum
only for this finite four-zero shortening problem.

On each retained support use an extended Hamming [16,11,4] sign code:
its five parity checks are the all-one word and the four coordinate-bit
functions on F_2^4. Its kernel has 2,048 words. Use amplitudes +/-1/4, so
every point has norm one. Within one support, distance four bounds the inner
product by (16-8)/16 = 1/2. Between supports, an intersection of at most
eight gives the same bound, independently of the sign-code coset. This
layer contributes **19,810,304** points.

Both constructions also use every 16 sign pattern on each affine
four-support {a,b,c,a xor b xor c} in F_2^6, restricted to coordinates
0,...,d-1 of the *new* ambient space. The amplitudes are +/-1/2. Every
triple belongs to a unique affine block, so distinct supports intersect in
at most two. There are 8,555 blocks in 61D and 9,145 in 62D, giving 136,880
and 146,320 points. All coordinate poles +/-e_i are included.
These layers are mutually compatible: the weight-16/weight-four cross
bound is 4/(4*2) = 1/2; all coordinates of a sparse point have magnitude
at most 1/2. Different layers are disjoint by their numbers of nonzero
coordinates.

## 61D: odd parity repairs the dense/sparse angle

Shorten C to zero on coordinates 12 and 52. This loses two independent
dimensions, leaving 2^26 words of distance at least 16. Write their two
other removed bits as i = c_0 and j = c_1, and their head signs as s.
Use the existing exact latitude map:

    j = 0:  (s/8,             (-1)^i/4)
    j = 1:  (sqrt(17)*s/32,   (-1)^i/16).

Both have norm one. The universal dense/dense proof in
`../dim49plus/latitude_compression.py` checks every ordered label pair and
every permitted head distance with exact rational/radical comparisons.
It applies after any permutation of the binary coordinates.

The absolute-value cross bound with a weight-16 sparse point on the second
latitude is 16*sqrt(17)/128 > 1/2. This was the obstacle to using that layer.
Translate the extended Hamming sign code by a single bit: **every sparse
sign word now has odd parity**. For every selected support S, every dense
head restriction has even parity because S belongs to the self-orthogonal
code C and avoids all four removed coordinates. Thus the dense and sparse
signs disagree at least once. The cross bounds become

    j = 0:  14/32 <= 1/2,
    j = 1:  14*sqrt(17)/128 < 1/2.

The verifier checks the parity condition on every shortened generator row
and every support. It independently expands all 2,048 shifted sign words
as integer coordinates and checks their entire Gram matrix. A negative
control uses an even-parity coset containing the restriction of an actual
second-latitude dense word; the resulting sparse point has inner product
16*sqrt(17)/128 > 1/2. The odd coset is essential to this argument.

For the four-support layer there is at most one tail coordinate. Its
dense/sparse bound is at most (3/8+1/4)/2 = 5/16 on the first latitude.
On the second, the two cases give sqrt(17)/16 and
(3*sqrt(17)/32+1/16)/2, both below 1/2. The coordinate poles are compatible
as well.

    2^26 + 9,673*2,048 + 8,555*16 + 122 = 87,056,170.

## 62D: a hexagon replaces four binary coordinates

The projection of C on F has rank four, so each of the 16 labels has
exactly 2^24 words. Set bit c_52 to zero and retain the following six
labels (bits c_0,c_1,c_12, with c_0 the low bit), in cyclic order:

    0, 1, 3, 7, 6, 4.

Keep each head sign with amplitude 1/8 and replace the four removed
coordinates by a two-coordinate hexagonal tail, also scaled by 1/8:

    (2,0), (1,sqrt(3)), (-1,sqrt(3)),
    (-2,0), (-1,-sqrt(3)), (1,-sqrt(3)).

Every tail has squared norm 4 before scaling, so the full squared norm is
(60+4)/64 = 1. For labels at Hamming distance h, the head distance is at
least 16-h. Its inner product is at most (28+2h)/64. The scaled tail Gram
matrix has entries 4 on the diagonal, 2 for cyclic neighbours, -2 for
second neighbours and -4 for opposite vertices. These are respectively
at most 4-2h for h = 0,1,2,3. Hence every dense/dense inner product is at
most 32/64 = 1/2. The verifier checks all 36 ordered label cases. A
deliberate misordering fails, providing a negative control.

This retains **6*2^24 = 100,663,296** dense points. A weight-16 head
support has cross inner product at most 16/(8*4) = 1/2, so the ordinary
even extended Hamming signs work here. A four-support can meet either,
both, or neither tail coordinate. Checking those cases gives a bound of
at most max(1/4,5/16,(3+sqrt(3))/16) < 1/2. The poles are compatible.

    6*2^24 + 9,673*2,048 + 9,145*16 + 124 = 120,620,044.

## Replay and verification scope

From the package directory:

```bash
python verify_significant.py --cpp
python verify.py
```

`compressed_layers.py` rebuilds all proof inputs and writes the deterministic
`compressed_layers.json` receipt. No new binary certificate is required.
`verify_significant.py` runs it in a temporary copy, checks the full receipt
against the installed one, and checks the registered totals. With `--cpp`,
the existing `../dim61_96/enumerate_binary.cpp` independently enumerates the
same 268,435,456 words by scalar Gray-code updates. Its full histogram and
the SHA-256 of all 9,673 selected supports must agree with the NumPy replay.
The default package `verify.py` includes the new proof and retains all older
construction replays. The package SHA-256 manifest covers these files.

All angle comparisons are exact. Integer coordinate dot products stay
within int16 bounds; radical comparisons use rational squares. There is
no floating-point angle tolerance or trusted search score. The finite
support and sign-code checks, together with the layer inequalities above,
prove the full populations without materializing hundreds of millions of
real-coordinate points.

The earlier PackingStar-inspired fill/correct searches motivated the
continuation. These two larger gains come from parity information and
coordinate compression; they do not depend on PackingStar's implementation
or its reported 26--31D records. Exploratory dual-shell owner-sharing and
affine-coset searches are superseded by the simpler construction above.
