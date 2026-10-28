# Fano-grid and cube-edge latitudes in 38D

The construction gives **611,522** points, **1,128** more than the retained
610,394-point certificate. The latitude population increases from 5,408 to
6,656; deleting 120 original cap points pays for the larger layer.

| Component | Points |
|---|---:|
| Retained Leech caps | 589,560 |
| Lin equator | 13,710 |
| Rotated flat axes and coordinate axes | 1,596 |
| Fano-grid latitudes | 6,272 |
| Cube-edge latitudes | 384 |
| Total | **611,522** |

The cap and equator construction, attribution, lattice minimum proof and
normalizations are retained in [the original proof](../source/proof.txt).
The original certificate is unchanged. This continuation changes its latitude
layer and removes exactly the six caps of each of twenty body lines.

## The Fano grid

Split the fourteen tail coordinates into two sets of seven. Each carries the
seven lines of a Fano plane listed in `recipe.json`. For every pair of lines
`(i,j)`, use their union as a six-coordinate tail support. Assign the cell
the body-coordinate pair `(2c,2c+1)`, where `c=(i+j) mod 7`.

On the resulting eight-coordinate support, include all 128 sign patterns
of even parity, with every nonzero entry equal to 4, normalized by
sqrt(128). There are 49 supports and 49*128 = 6,272 points.

Each support has two body coordinates and six tail coordinates, so the body
has norm squared 32 and the tail 96. Distinct grid supports intersect in at
most four positions: cells sharing a row or column have different colors
and disjoint body pairs, while cells in different rows and columns share
only two tail positions. Their inner products are therefore at most 64.
On the same support, distinct even sign patterns differ in at least two
positions and also have inner product at most 64.

The certificate encodes this as 1,568 four-point families
`(+/-u,+/-2Z)/sqrt(128)`, where Z has six entries +/-2. Body-root parity
matches the tail-sign parity. Only fourteen body lines are used.

## The cube edges

The remaining latitude points use body coordinates 14,15,16, which are
disjoint from every grid body. In each tail half, distinguish three
positions s0,s1,s2 and let O be its other four positions.

For each cyclic triple (i,j,k)=(0,1,2),(1,2,0),(2,0,1), and each of the
64 sign assignments on all six tail positions except si, put

    u = 4 * product(y[t] for t in O) * (y[sj] e[14+j] - y[sk] e[14+k])
    point = (u, 4y) / sqrt(128).

This gives 192 points in each half, or 384 in total. They use six additional
body lines. Their compatibility is also transparent:

- With grid points the bodies are orthogonal and the tails overlap in
  at most three positions, giving product at most 48 before normalization.
- Between different tail halves the tail product is zero and the body
  product is at most 32.
- With the same omitted position, one flipped ordinary sign reverses u,
  and one flipped distinguished sign changes one body component. Either
  compensates the otherwise excessive tail product. Two or more flips
  already give a tail product at most 32.
- With different omitted positions, agreeing on all five shared signs
  gives tail product 80 and body product -16. Any disagreement gives
  tail product at most 48 and body product at most 16.

Every distinct pair thus has inner product at most 64/128 = 1/2.

## Retained layers

For all latitude points, the exact checks establish `|Z.w| <= 8` for
every original cap head w, and the same bound for every retained rotated
axis. All six caps belonging to every latitude body line are removed.
Every remaining cap has a different Leech line and hence body product of
absolute value at most 16. Its product with a latitude point is at most

    (16 + 16 sqrt(2)) / (32 sqrt(6)) <= 1/2.

The two verifiers prove this inequality with different integer comparisons.
Equatorial compatibility follows from Cauchy-Schwarz because the normalized
latitude body has length 1/2. Coordinate-axis inequalities reduce to
coordinate squares. The four layers have different body norm fractions,
which proves that they are disjoint.

## A precise ceiling for this support model

These 63 tail supports comprise the 49 Fano-grid unions and the fourteen
six-subsets of the two seven-coordinate halves. With this exact support
set, entries +/-4 in the tail and body norm squared 32, **6,656 latitude
points is an upper bound**, even allowing arbitrary real body vectors.

For a fixed oriented six-sign head, two compatible bodies must be
antipodal; hence there are at most two. The 49 grid supports contribute
at most 49*64*2 = 6,272.

In one seven-coordinate half, a signed six-support is an edge of the
seven-dimensional cube: its missing sign is free. At any cube vertex,
bodies of incident selected points have pairwise products at most -16.
For m such norm-32 bodies, their sum has norm squared at most
32m-16m(m-1), so m<=3. There are 128 vertices and each point is counted at
two vertices. Thus each half has at most 192 points. The construction
attains both bounds.

This ceiling concerns the stated latitude support model. It does not
bound kissing numbers, other head supports, other norm splits, or the
best possible cap deletion count. Those remain directions for larger gains.

## Reproduction

Run from the package root:

```bash
python verify.py
python refill/negative_controls.py
```

`construct.py --check` rebuilds `certificate.json` byte for byte from
`recipe.json`. `verify.py` checks compressed family products.
`verify_coordinates.py` independently reconstructs the Golay code and
Leech shell, replays all original layer facts, expands all 6,656 latitude
points and checks all **22,147,840** distinct latitude pairs.

`verification.json` and `coordinate_verification.json` archive passing CI
receipts for this exact certificate.

`negative_controls.py` requires both proof representations to reject a
duplicated family, a wrong cube-body sign, an undeleted body owner and a
false count. The package wrapper performs isolated replays, so receipts
never modify the archived source. Python and NumPy suffice; no optimizer,
network retrieval, cached shell, or floating-point tolerance is used.
