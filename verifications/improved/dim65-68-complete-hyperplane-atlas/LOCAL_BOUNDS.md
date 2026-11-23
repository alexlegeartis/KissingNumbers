# Exact caps for two local fibre models

These caps concern particular layers of the construction. They are not
upper bounds on kissing numbers in any dimension.

## One body sign coset of RM(2,5)

On one affine hyperplane and at one fixed head, the desired Hamming
minimum is at least nine. Inside a single translate of RM(2,5), all
differences have weights in `{0,8,12,16,20,24,32}`, so the allowed
distances between distinct words are `{12,16,20,24,32}`.

Put `k_j(d)=K_j(d)/binom(32,j)`, with the usual binary Krawtchouk
polynomials. The polynomial

```
1 + (5517/340) k_1 + 124 k_2 + (83669/170) k_3 + 899 k_4 + (35061/68) k_9
```

is nonpositive at each allowed distance and has value exactly 2048 at
zero, as can be checked in `local_bounds.py`. Its Krawtchouk
coefficients are nonnegative with constant coefficient one, so summing
over ordered pairs bounds the code size by 2048, even if it is nonlinear.
The trace-quadratic code in this certificate attains this cap.
Separately, exhaustive 21-bit quotient BFS (`cover32.cpp`) gives
covering radius eight for this particular `[32,11,12]` trace code.
Thus no single additional body word, even outside that Reed–Muller
coset, can augment an intact trace-code fibre at the same head.

## One support with arbitrary real heads, body weight 12–15

Use a fixed support of (w) body coordinates with signs of magnitude
one, and heads of squared norm (16-w) in at most four real
coordinates. Pairwise dot products must not exceed eight. The whole
fibre has at most 2048 vectors, even with arbitrary real head positions.
Let (d) denote Hamming distance on the body.

- (w=15): d ≤ 2 is impossible even for opposite heads. The
  radius-one balls of the body words in the 15-cube are disjoint;
  their volume is 16.
- (w=14): d ≤ 1 is impossible; (d=2) forces opposite heads.
  A word can have at most seven distance-two neighbours, since the
  respective two-site flips must be pairwise disjoint. A cube point
  is in at most two radius-one balls because three pairwise opposite
  heads are impossible. Two adjacent balls of distance-two centres
  have exactly two common points. Thus the union of balls has at
  least (15N-7N=8N) points in the 14-cube.
- (w=13): (d=0) is impossible, (d=1) forces opposite heads,
  and at (d=2) unit-normalized heads have dot at most (-1/3).
  An occupied cube point has at most one occupied neighbour, whose
  head must be its opposite; an empty point has at most four occupied
  neighbours by the simplex bound on their mutually negative heads.
  Counting neighbour incidences gives `13N ≤ 4(2^13 − N) + N`.
- (w=12): a body word carries at most two heads, then opposite.
  On the neighbours of a body word all normalized heads have
  nonpositive dot, so an unoccupied word has at most eight such
  heads by the orthoplex bound in four dimensions. A singly
  occupied word has at most four neighbouring heads: each has
  dot at most (-1/2) with the central head, and their sum has
  squared norm at most their number. A doubly occupied word has
  none. If (n_1,n_2) count body words with one or two heads, then
  (N=n_1+2n_2) and `12N ≤ 8(2^12 − n_1 − n_2) + 4n_1`.

The bounds show why enlarging head geometry alone does not increase the
capacity of these particular fibres. They do not exclude joint support
changes or several different amplitude patterns.
