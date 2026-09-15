# -*- coding: utf-8 -*-
"""The axis of the dimension-29 configuration is at its exact ceiling: 10 points, no more.
(python axbound29.py -- exact, sympy only, no data files.)

Once the layer occupies all 32 deep holes, an axis point (0, 2a) with |a| = 1 in the D5 frame has
to clear two families at once:

    the layer  2 sqrt(5/2) <a,w> <= 2 for every w = s/sqrt5, s in {+-1}^5   <=>   ||a||_1 <= sqrt2
    the caps   (4/sqrt3) <a,z> <= 2 for every cap direction z = r/sqrt2     <=>   |a_i| + |a_j| <= sqrt(3/2)

The first is the whole point: the 32 sign patterns are exactly the extreme points of the
l_infinity ball, so "cosine at most 2/sqrt10 with all of them" is literally an l_1 bound.  Both
together force every admissible a to be within 22.7 degrees of a coordinate axis, which leaves
room for one axis point per +-E_i and no more.

Step 1.  Write m = max |a_i| and let the second largest be t.  Then
             1 = sum a_i^2 <= m^2 + t (||a||_1 - m) <= m^2 + (sqrt(3/2) - m)(sqrt2 - m),
         because the tail is bounded coordinatewise by t and t <= sqrt(3/2) - m is the cap
         constraint on the two largest coordinates.  So
             2 m^2 - (sqrt(3/2) + sqrt2) m + (sqrt3 - 1) >= 0,
         and since m >= 1/sqrt5 lies strictly above the smaller root, m >= m0, the larger root.

Step 2.  m0 > sqrt3/2, hence 2 m0^2 - 1 > 1/2.  Two admissible points whose largest coordinate is
         the same index with the same sign have
             <a,b> >= m0^2 - sqrt((1-m0^2)(1-m0^2)) = 2 m0^2 - 1 > 1/2,
         which the axis constraint <a,b> <= 1/2 forbids.  So at most one point per (index, sign),
         i.e. at most 2*5 = 10; and m0 > 1/sqrt2 makes the largest coordinate unique, so the ten
         blobs are disjoint and every admissible point lies in exactly one.

Step 3.  +-E_i are admissible (||E_i||_1 = 1 <= sqrt2 and |E_i| + 0 = 1 <= sqrt(3/2)), pairwise at
         inner product 0 or -1 <= 1/2, so the bound is attained.  10 is the maximum.

For comparison, the published 40-point axis (Ma et al. 2025, arXiv:2511.13391) keeps only 8 of its
points against the same layer
(research/collab2531/dim29/price29.py), so replacing it with +-E_i is worth 2 points on top of the
layer's 128: the gain over the published 209496 is 128 + 10 - 40 = 98.
"""
from __future__ import print_function
import itertools as it
import sympy as sp

L1 = sp.sqrt(2)                      # the l_1 bound from the 32 deep holes
PAIR = sp.sqrt(sp.Rational(3, 2))    # the two-coordinate bound from the 40 cap roots

# ---- the two regions really are what the docstring says
m = sp.Symbol('m', nonnegative=True)
H = sp.sqrt(sp.Rational(5, 2))
a = sp.symbols('a0:5', real=True)
# layer:  2 H <a, s/sqrt5> <= 2  for all sign patterns s  <=>  max_s <a,s> <= sqrt5/H = sqrt2
assert sp.simplify(sp.sqrt(5) / H - L1) == 0, 'the deep-hole constraint is exactly ||a||_1 <= sqrt2'
# caps:  (2/sqrt3)*2 <a, r/sqrt2> <= 2  <=>  |a_i| + |a_j| <= sqrt2 * sqrt3/2 = sqrt(3/2)
assert sp.simplify(sp.sqrt(2) * sp.sqrt(3) / 2 - PAIR) == 0, 'the cap constraint is |a_i|+|a_j| <= sqrt(3/2)'
print('  region: ||a||_1 <= sqrt2 = %.9f   and   |a_i| + |a_j| <= sqrt(3/2) = %.9f'
      % (float(L1), float(PAIR)))

# ---- step 1: the quadratic and its larger root
q = sp.expand(2 * m ** 2 - (PAIR + L1) * m + (sp.sqrt(3) - 1))
assert sp.simplify(q - (m ** 2 + (PAIR - m) * (L1 - m) - 1)) == 0, 'the quadratic IS  m^2 + (sqrt(3/2)-m)(sqrt2-m) - 1'
roots = sorted(sp.solve(sp.Eq(q, 0), m), key=lambda r: sp.N(r, 40))
m0 = sp.simplify(roots[1])
assert sp.N(roots[0], 40) < sp.N(1 / sp.sqrt(5), 40), 'the smaller root is below the floor 1/sqrt5'
print('  step 1: 2m^2 - (sqrt(3/2)+sqrt2) m + (sqrt3 - 1) >= 0 with m >= 1/sqrt5 = %.6f'
      % float(1 / sp.sqrt(5)))
print('          smaller root %.9f < 1/sqrt5, so max |a_i| >= m0 = %.9f  (%.3f degrees off axis)'
      % (float(roots[0]), float(m0), float(sp.deg(sp.acos(m0)))))

# ---- step 2: one point per blob
assert sp.N(m0, 40) > sp.N(sp.sqrt(3) / 2, 40), 'm0 > sqrt3/2, so 2 m0^2 - 1 > 1/2'
assert sp.N(m0, 40) > sp.N(1 / sp.sqrt(2), 40), 'm0 > 1/sqrt2, so the largest coordinate is unique'
lo = sp.simplify(2 * m0 ** 2 - 1)
assert sp.N(lo, 40) > sp.Rational(1, 2)
print('  step 2: two admissible points in one +-E_i blob are at inner product >= 2 m0^2 - 1 = '
      '%.9f > 1/2, which the axis constraint forbids  ->  at most 1 per blob, at most 10 in all'
      % float(lo))

# ---- step 3: +-E_i attain it
E = [sp.Matrix([1 if k == i else 0 for k in range(5)]) for i in range(5)]
AX = [s * e for e in E for s in (1, -1)]
for v in AX:
    assert sum(abs(c) for c in v) <= L1, 'axis point outside the l_1 region'
    assert max(abs(v[i]) + abs(v[j]) for i, j in it.combinations(range(5), 2)) <= PAIR
for p, r in it.combinations(AX, 2):
    assert sp.simplify(p.dot(r)) <= sp.Rational(1, 2), 'two axis points at inner product > 1/2'
print('  step 3: the 10 points +-E_i are admissible and pairwise at inner product 0 or -1')
print('EXACT: the axis of this configuration has at most 10 points, and +-E_i attain it.')
