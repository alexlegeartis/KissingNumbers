# Affine-hyperplane heads and signed shelves in 63–68 dimensions

These are explicit kissing configurations with the following counts.
The comparisons are with the preceding verified repository constructions,
not assertions that the new numbers are optimal.

| Dimension | Before | New lower bound | Gain |
| ---: | ---: | ---: | ---: |
| 63 | 216,972,548 | 217,044,228 | 71,680 |
| 64 | 418,347,904 | 418,478,976 | 131,072 |
| 65 | 418,741,122 | 418,874,240 | 133,118 |
| 66 | 426,343,300 | 429,011,844 | 2,668,544 |
| 67 | 442,334,086 | 448,201,606 | 5,867,520 |
| 68 | 472,598,888 | 479,699,848 | 7,100,960 |

## The two new 68D ingredients

Start with the dense 80-coset Reed–Muller sign code, the signed affine
2-flats, the poles and the sparse supports of the
[previous construction](../dim63-67-signed-overlap/). The first new
ingredient adds 4,981 sparse supports, each of weight 14 in the original
64 coordinates and weight 2 in the new tail. On each combined 16-site
support, use a support-specific coset of an extended Hamming sign code
with 2,048 words. Use only support intersections of size at most ten,
and on a size-nine or size-ten intersection enforce opposite sign
parity on a shared eight-site subset. This forces at least one sign
disagreement. The exact verifier reconstructs all 30,667 such edges and
completes every partial dual check space using balanced fibres.

The second ingredient uses 94 affine hyperplanes of the original
six-bit coordinate space. Fix a hyperplane `H`, identified with
`F_2^5` by deleting the lowest nonzero normal coordinate. Put
`sqrt(2)` times a 32-bit sign word on `H`, zero elsewhere, then attach
one of 24 heads: eight axial points `±2sqrt(2)e_i` and 16 cube corners
`sqrt(2)(±1,±1,±1,±1)`. These are the vertices of a regular 24-cell
with squared norm eight and maximum distinct dot product four. The 24-cell
is especially natural here: for a four-dimensional head of squared norm
eight, bounding the sum of its two largest absolute coordinates by
\`2 sqrt(2)\` forces equality. Indeed, with magnitudes \`a>=b>=c>=d\`,
\`(a+b)^2-(a^2+b^2+c^2+d^2)=2ab-c^2-d^2>=0\`, with equality only
on an axis or when all magnitudes agree.

For each hyperplane and head, the body signs are one coset of the same
size-2,048 binary code, spanned by the affine functions and the trace
quadratics `Tr(a x^3)` over `F_32`. The exact replay checks this is a
[32,11,12] code in RM(2,5). The 24 heads on a given hyperplane use
24 distinct quadratic phases modulo this code. Thus any two body words
for different heads on the same hyperplane differ in at least eight
places. On distinct hyperplanes the intersection is empty or has size
16. Distinct heads then give dot at most `2*16+4=36`; equal heads
need one disagreement, forced by opposite fixed parity on an affine
3-flat or on the full intersection. The new vectors have norm squared
72. The contribution in 68D is `94*24*2048=4,620,288`. In 66D and 67D
retain just the four or six axial heads, respectively, and restrict
the sparse stock to the tails present in that dimension.

For old/new pairs, the replay directly verifies the minimum distance
four on the hyperplanes against the restrictions of the dense sign
code, which gives dot at most `24 sqrt(2)`. For an old sparse support,
the margin in units of `sqrt(2)` is 12. For each relevant intersection,
the verifier uses the largest possible contribution of the new head
on that sparse tail. Where one disagreement is needed, it inserts a
shared weight-eight check of opposite parity. Where two are needed,
the intersection is the whole old 16-site support: check equal total
parity as well, so the number of disagreements is even. These
additional checks are inserted without changing any earlier syndrome
or the 2,048 signs per support. The affine 2-flats and poles obey the
mixed bounds directly.

## Completing the sparse shelves

After the hyperplane heads are installed, add 64 further weight-16
zero-tail supports and 128 weight-12 body supports with all four tail
positions nonzero. The last 192 supports add 393,216 points in 68D.
Their literal proof is in `joint_witness.txt`. The entire sparse
stock now has 31,512 guarded overlap edges. All 192 are compatible
with the 94 hyperplanes by the raw intersection inequalities alone.
The zero-tail additions also pass to 64–67D. In 63D, use the
coordinate-equality section of the enlarged equator.

## Separate 65D polar extension

The equator is as in the earlier certificate. Use three 65,536-word
cells on each side at height `sqrt(17)`, with a new choice of
representatives and sign shift, and eight 128-word cells on each side
at height 6. The two old axial poles are removed. The small polar
verifier checks all binary minimum-distance gates and gives a net gain
of 2,046 beyond the enlarged 65D equator with the previous polar stock. This polar construction
is not assumed to coexist with the higher-dimensional heads.

## Reproduce

Run `python3 verify.py` with a C++20 compiler on x86-64 with BMI2.
The older gzip witness is read from its adjoining directory, and all
literal witness SHA-256 digests are pinned in the Python driver. No
network access, floating-point tolerances or unpublished lattice
enumerations are needed. The historical external comparison is
[Edel–Rains–Sloane, On Kissing Numbers in Dimensions 32 to 128](https://arxiv.org/abs/math/0207291)
and its elementary dimension extensions; the before/after column
above is the separate comparison with the immediate repository head.
