# Upper bounds for one-axis lifts into 26D and 27D

A one-axis lift starts with a unit-sphere code C in dimension d. It
retains the unselected points on the equator. For each selected parent
x, it replaces (x, 0) by the two points

```text
(cos(a_x) x, sin(a_x))  and  (cos(a_x) x, −sin(a_x)).
```

Each parent can have its own angle a_x in [0, π/2]. The two children of
one parent have inner product cos(2a_x), so their compatibility
requires a_x ≥ π/6. At most one parent can use a_x = π/2, since that
pair consists of the two poles.

For two nonpolar parents x and y, compatibility of their positive
children implies

```text
x · y ≤ (1/2 − sin(a_x) sin(a_y)) / (cos(a_x) cos(a_y)) ≤ 1/3.
```

For the second inequality use
sin(a) sin(b) + (1/3) cos(a) cos(b) ≥ 1/2
on [π/6, π/2]². As a function of either variable with the other
fixed, the left side is concave. Its minimum therefore occurs at a
corner, where its values are 1/2, 1/2, 1/2 and 1.

## Into dimension 27

The nonpolar parents form a spherical code in R²⁶ with maximum product
1/3. Use the following rational Delsarte polynomial.

```text
f(t) = (t − 1/3)(t + 9/14)²(t + 5/18)²(t − 1/38)².
```

It is nonpositive for −1 ≤ t ≤ 1/3. Its normalized Gegenbauer
coefficients in dimension 26 are nonnegative, with positive constant
coefficient. The addition formula gives
Σ(x,y in X) f(x·y) ≥ f₀ |X|² for the set X of nonpolar parents.
The off-diagonal sign bounds the same sum above by |X| f(1). The exact ratio f(1)/f₀ is
5,363,432,606 / 606,523 < 8,843. Thus at most 8,842 nonpolar parents
can be used, plus at most one polar parent. Starting at 201,192 gives
at most **210,035** points in dimension 27. The bound applies to
every parent of that cardinality and to independently chosen symmetric
latitudes.

## Into dimension 26

Start with the repository's 197,586-point critical-cluster witness in
dimension 25. The Leech-shell parents must have mutual products in
{−1, −1/2, −1/4, 0, 1/4}, since the only larger distinct-point Leech
product is 1/2. In normalized Gegenbauer polynomials for dimension 24
use f₀ = 1 and the following other nonzero coefficients.

| Degree | Coefficient |
| ---: | ---: |
| 1 | 928/77 |
| 2 | 5083/77 |
| 3 | 25826056/129591 |
| 4 | 27600/77 |
| 9 | 27747200/129591 |

The polynomial vanishes at every allowed product. Its ratio f(1)/f₀
is 9360/11, so at most 850 nonpolar parents come from the Leech shell.

Each of the witness's 926 circles has a centre of squared norm 20
and squared radius 12 in the squared-norm-32 convention. Any three
points on one such circle have a pair with angular separation at most
120 degrees. That pair has normalized product at least
(20 − 12/2)/32 = 7/16 > 1/3, so at most two nonpolar parents can
come from each circle. Count all 247 general points as potentially
eligible and allow one polar parent. This yields at most
197,586 + 850 + 2·926 + 247 + 1 = **200,536** in dimension 26.
[Source geometry and counts](../../verifications/improved/dim25-critical-cluster/)

## Replay

Run `python verify_symmetric_lifts.py`. The checker verifies the two
polynomial certificates with the standard library's rational arithmetic.
This checker covers the symmetric two-child construction above.

## Arbitrary latitudes along the parent rays

The following bounds allow different positive and negative latitudes,
three- and four-point replacements, and motion of singly represented
parents. Every child of x has coordinates (cos(a) x, sin(a)) with
−π/2 ≤ a ≤ π/2. The projection of each child onto the original space
therefore follows its original parent's ray.

Each parent's children are separated by at least π/3 in latitude.
A parent with no child at absolute latitude at least π/4 has at most
two children. If it has two, one is positive and one negative, with both
absolute latitudes between π/12 and π/4. Call such parents balanced.
For any two balanced parents the positive children impose

```text
x · y ≤ g(a,b) = (1/2 − sin(a) sin(b)) / (cos(a) cos(b)).
```

For fixed b, the derivative of g in a has the sign of
(1/2) sin(a) − sin(b), which is increasing on the interval.
Thus the maximum on the square is at a corner. Corner evaluation
gives g ≤ 2√3 − 3 < 15/32.

Now consider parents with a child of absolute latitude at least π/4.
Assign each such parent to the sign of one such child. Within each
sign, the nonpolar children imply pairwise nonpositive products for
their parent directions, since sin(a) sin(b) ≥ 1/2. The orthoplex bound
gives at most 2d such parents per sign, and the pole allows one more
per sign. Three children contribute a gain of two. Four children
require both poles and can occur for at most one parent. The total
contribution of these parents is at most 2(4d + 2) + 1.

### Counts for the supplied parent codes

In the 25D critical-cluster code, the above 850-point bound still
applies to the balanced parents from the Leech shell. Exact checks of
the field coordinates find a pair of product above 15/32 in each of
the 926 circles, so each circle contributes at most two balanced
parents. The remaining 247 points can all contribute.

For each of the two 26D certificates, the verifier checks that the
integer equator consists of Leech points. On each eight-point
midpoint fibre, the residual vectors form four orthogonal axes with
both signs. Two different residual axes give full normalized product
1/2 because they share the same squared-norm-64 projection. Hence
at most two balanced parents come from each such fibre. The remaining
columns include all integer and rational points outside the stated layers.

| Parent witness | Balanced Leech cap | Balanced cap off the equator | Gain from the larger latitudes | Lift upper bound |
| --- | ---: | ---: | ---: | ---: |
| 197,586-point critical-cluster code in 25D | 850 | 2·926 + 247 | 205 | **200,740** |
| Takhanov’s 201,192-point code in 26D | 850 | 2·1,024 + 888 | 213 | **205,191** |
| Alternate 201,192-point code in 26D | 850 | 2·1,088 + 376 | 213 | **204,807** |

Run `python verify_radial_replacements.py` from this directory with
NumPy. Supply `--published26 /tmp/d26_coordinates.npz` to include the
pinned Takhanov parent. The replay reconstructs the 25D circles and
uses exact comparisons in Q(√2, √3). The two 26D fibre decompositions
use integer arithmetic. Parent configurations and other geometric
premises are documented in the source certificates linked by the
[research README](README.md).
