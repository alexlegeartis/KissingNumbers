# Complete exact proof for the 197,586-point certificate

The input is only `source/certificate25.json`. The checker imports no
construction/search code or stored shell. Its only external dependency is
NumPy, used with bounded signed-integer products. Acceptance uses no
floating-point geometry or tolerance.

## 1. Independently reconstructed shell

Construct all binary multiples of

    g(x) = 1 + x^2 + x^4 + x^5 + x^6 + x^10 + x^11

with multiplier degree at most 11, and extend each by a parity coordinate.
The checker verifies 4,096 distinct words, weight distribution

    0:1, 8:759, 12:2576, 16:759, 24:1,

and pairwise orthogonality of the 12 basis words over F_2. The construction
is linear, so distinct words have Hamming distance at least 8, and every
codeword is orthogonal to every codeword.

The independently generated shell has three families:

1. Two coordinates of absolute value 4: 1,104 vectors.
2. Absolute value 2 on an octad, with an even number of minuses:
   759 * 128 = 97,152 vectors.
3. A codeword sign vector with one coordinate multiplied by -3:
   4,096 * 24 = 98,304 vectors.

All 196,560 vectors are checked distinct and of squared norm 32. Their
mutual products are at most 16 by this complete family argument:

- Family 1 with family 1: distinct vectors have product at most 16.
- Family 1 with family 2: at most 4*(2+2) = 16.
- Family 1 with family 3: at most 4*(3+1) = 16.
- Family 2 on distinct octads: code distance at least 8 limits the
  intersection to four positions, giving at most 4*4 = 16. On the same
  octad, distinct even-minus sign patterns differ at least twice, giving
  at most 32 - 8*2 = 16.
- Family 2 with family 3: if the triple is outside the octad, the bound is
  16. Otherwise code orthogonality and the even-minus rule give odd parity
  of coordinate-product signs after the triple reversal. At least one
  product is negative, so the bound is 2*(3+7-2) = 16.
- Family 3 with family 3: let d be the distance of the two original
  codewords and z_i their original sign products. For the same exceptional
  position i, the product is 24-2d+8*z_i, at most 16 when distinct since
  d >= 8. For different positions i,j, it is 24-2d-4*z_i-4*z_j. This is
  exactly 16 if d=0, and at most 24-16+8 = 16 if d >= 8.

These checked finite code facts and elementary inequalities prove all
shell/shell bounds without an assumed lattice-minimum theorem.

## 2. Circles and retained-shell pairs

For each recorded pair u,u+v, the checker verifies that both endpoints are
in the reconstructed shell, u.v = -24 and v.v = 48. Put c=(2u+v)/2.
Then c.c = 20 and c.v = 0. A circle point is

    (c + cos(theta)*v/2, 2*sqrt(3)*sin(theta)).

Its squared norm is 32. For each initial phase, a rational half-angle
p/q rotates its cosine and sine by

    cos(delta) = (q^2-p^2)/(q^2+p^2),
    sin(delta) = 2pq/(q^2+p^2).

All reconstructed coordinates lie in Q(sqrt(2),sqrt(3)). The checker also
checks each squared norm directly in this field. Both endpoints of every
circle are removed. Exact signs check -1 <= cos(theta) <= 1, so the foot
is on their closed segment. Its product with any retained shell vector is
at most 16 by convexity and Section 1; the last coordinate contributes zero.

For two different circles of one class, the common circular plane is
orthogonal to both centres. The checker directly checks
(2u+v).(2u'+v) <= 16, so the centre product is at most 4. Both circular
components have radius sqrt(12). Cauchy-Schwarz bounds their product by
12, and the full product by 16 for any phases. These comparisons cover
exactly 1,658,754 exceptional pairs in this witness.

## 3. All remaining exceptional pairs

Represent each coordinate by four integer coefficients for
1,sqrt(2),sqrt(3),sqrt(6), divided by its positive denominator.
Integer square roots at precision 10^80 give outward radical endpoints.
Integer floor and ceiling operations enclose each coordinate in [L/S,U/S]
with S=10^8. The largest enclosure width is at most two integer units.

Write x=(L_x+e_x)/S with 0 <= e_x <= E_x, coordinatewise. An upper bound
for S^2*x.y is

    L_x.L_y + E_x*||L_y||_1 + E_y*||L_x||_1 + 25*E_x*E_y.

The checker bounds every signed-int64 multiplication, its partial sums,
and all correction additions before evaluation. All unresolved comparisons
are evaluated using arbitrary-precision integer field arithmetic.
Together with Section 2's analytic same-class pairs, an explicit coverage
identity accounts for all 3,025*3,024/2 = 4,573,800 exceptional pairs.
There are 2,915,046 interval-screened pairs and 3,292 exact fallbacks.

To decide the sign of a+b*sqrt(3), check the coefficient signs first.
Opposite signs reduce to comparison of a^2 and 3*b^2. For
x+y*sqrt(2), with x,y in Q(sqrt(3)), equal signs settle the result.
Opposite signs reduce to the sign of x^2-2*y^2, again in Q(sqrt(3)).
The irrationalities ensure that the tested nontrivial zero cases cannot
occur. Squaring is only used after the relevant signs have been checked.

## 4. Rational general points and retained shell

Each of the 247 general points has explicit integer numerators and a
positive denominator. Its norm is checked exactly. The same outward
coordinate enclosures give a rigorous upper bound on its product with
every retained integer shell vector. Every inconclusive case would be
checked by an arbitrary-precision integer dot product. All 247*194,561 =
48,056,567 comparisons pass, with 12,492 exact fallbacks.

## 5. Count, uniqueness and kissing normalization

The union of all deleted shell endpoints and additional deletion lists
contains exactly 1,999 vectors. Every listed deletion is checked to be a
shell vector; duplicates in deletion lists are counted only once. Hence

    196560 - 1999 + 3*926 + 247 = 197586.

The common squared norm is 32 and every distinct-index pair has product
at most 16, so coincident points are impossible. Every squared distance
is at least 32. Scaling by 1/sqrt(8) gives centres at distance 2 from the
origin and mutual distances at least 2, proving the unit-ball kissing
lower bound. Scaling by 1/sqrt(32) gives the equivalent unit-sphere code.

The numerical search and provenance hashes are not assumptions in this
proof. The self-tests provide regression and corruption checks, not a
replacement for these exact arguments.
