# Nested affine face trades in 66–68 dimensions

In 68D the four-site supports form an exact Steiner quadruple system:
every triple of coordinates occurs in exactly one support.

A four-site vector has entries ±1 on its support and zero elsewhere.
Take all sixteen signs on each selected four-site support. Two such
supports are compatible exactly when their intersection has at most two
coordinates. The complete hyperplane atlas includes all 10,416 affine
2-flats in the 64 binary body coordinates as four-site supports.

A face trade removes an all-body affine 2-flat and sends distinct faces
(three of its four points) to distinct tail coordinates. The old
all-body flat is the unique one containing any of its faces. Faces
using distinct tails can share a body pair; faces using the same tail
must not reuse a pair. A face must also avoid using the horizontal pair
of an old two-tail support when its tail is among those two.

| Dimension | Old bound | Four-site supports | Gain | New bound |
|---:|---:|---:|---:|---:|
| 66 | 429,274,500 | 11,088 | 10,240 | 429,284,740 |
| 67 | 448,596,358 | 11,792 | 20,480 | 448,616,838 |
| 68 | 481,013,656 | 12,529 | 30,720 | 481,044,376 |

## A reusable eight-point template

Label the points (s,b) in F2² × F2. Each affine function
f(s)=l·s+c gives the four-point graph of f. For tail j in F2²,
omit the point with label s=j+T(l), where T(x,y)=(x+y,y). For every
nonzero h, h·T(h)=1. Thus if two distinct labels s and s' are prescribed
with two prescribed b values, their two compatible graphs have
gradients differing by the nonzero vector h perpendicular to s+s'.
Exactly one of the two omitted points lies on those labels. For each
tail, the faces cover every edge joining different labels exactly once,
leaving just the four edges of difference (0,1).

Place this template in F8² using the 64 affine lines with nonhorizontal
directions, where (0,y) is horizontal. For slopes s≠0 choose the
residual edge direction (s,s²). The graph of Frobenius x↦x² is an
additive 3-space, and each of its seven nonzero directions is one such
residual direction. Apply the same template to its eight cosets. The
eighth original line direction and the residual graph direction are
chosen as (a,0) and (a,7), with a²=7. Now use the same eight-point template on the eight horizontal
rows (constant x). The uncovered edges are exactly the three complete
matchings with differences 1 in the horizontal coordinate, (a,0),
and (a,7). Attach a complementary pair of tail tags to each matching,
and retain the all-tail four-set. This gives three tiers of one
eight-point rule, with no literal search witness. For 66D and 67D,
keep exactly the supports avoiding the removed tail coordinates.

In 68D the exact 12,529 supports exhaust all 50,116 coordinate triples.
This is optimal **within the four-site ±1 alphabet**, not an upper bound
on the kissing number. To extend the cap to partial signs, fix a global
sign assignment and retain points matching it on all four entries.
Retained supports are distinct and share no three coordinates, as their
shared signs agree. A random assignment retains each point with
probability 1/16. Hence any signed selection in the branch has size at
most 16 times 12,529, attained here.

## Exact replay

Run `python3 verify.py` on a BMI2 x86-64 machine with a C++20 compiler.
It first reruns the [complete hyperplane atlas](../dim65-68-complete-hyperplane-atlas/).
The composition check enumerates the complete resulting four-site stock
and checks uniqueness of every three-site subset. The other old types
are unchanged: against any new support the absolute dot is at most 4
for the norm-squared-16 and norm-squared-64 vectors, at most 2 for a
coordinate pole, and at most 5√2 for a norm-squared-72 atlas vector.
The rematched two-tail and pure-tail supports have mixed dot at most
4√2 against the old atlas. The deletion and projection leave every
other layer of the corresponding dimension untouched.
No floating-point calculations or search libraries are used in the replay.
