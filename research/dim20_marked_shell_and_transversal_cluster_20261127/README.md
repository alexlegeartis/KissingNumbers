# Dimension 20: marked-shell trade costs and an exact transversal cluster

No new lower bound is claimed. The comparator is the 19,448-point
[Cohn–Li construction](https://arxiv.org/abs/2411.04916v3): 760 signed
roots, 130 odd-signed Golay octads (128 points each), and 2,048 dense
flat signs. The marked trade below is larger than a uniform support
replacement: it permits arbitrary subsets of signs on every admitted
support.

## A half-dense cut cannot pay for a signed marked-octad refill

Fix any coordinate i. Remove the 1,024 dense points whose ith sign is
positive, and remove the 3,328 old octad points with ith coordinate +1.
Retain all other points. Consider arbitrary new signed weight-eight
points with coordinate +1 at i; the seven other coordinates have
amplitude 1. Their norm squared is 8.

The possible seven-site supports are exactly the 260 shapes checked
below. A new support cannot overlap an unmarked old octad in five or
more places: any choices of signs on this overlap extend to the old odd
parity, and the dot product is then greater than 4. A seven-site
support cannot overlap the seven leaves of a negatively marked old
octad in exactly six places: there is one unused old sign to adjust,
and the marked contribution is −1, giving a possible dot of 5. An
identical seven-site support is allowed subject to sign parity, which
we need not impose for the upper bound. The checker reconstructs the
shortened extended Golay code and computes these admissible shapes for
all 20 marks.

An exact four-subset cover bounds even fully flexible sign choices.
There are 256 bundles of five admitted supports with a common
four-subset and pairwise intersection exactly that four-subset. Each
support belongs to four or sixteen bundles (240 and 20 shapes,
respectively). Keep the 240 bundles containing exactly one of the
twenty degree-sixteen shapes. Every support remains in at least four
bundles. These finite facts can be replayed with integer set
operations only.

For one bundle, fix the signs on its four common residual coordinates.
Two points from different support shapes cannot both appear: their
inner product is 1 from the mark and 4 from the shared signs. From one
shape the three remaining signs give a binary code with distance at
least 2, of size at most four. There are 16 common sign patterns, so
one bundle admits at most 64 points. Counting bundles with weight
1/4 proves that the entire signed refill has size at most
`240·64/4 = 3,840`. Consequently this whole trade gives **at most
18,936**, which is 512 below the incumbent. No numerical optimizer
or full-packet assumption is used.

There is a separate cheap obstruction for sparse new directions with
the entire octad layer retained. Every four-set of retained coordinates
lies in a retained octad: the five octads through that four-set have
disjoint complementary tetrads partitioning the other twenty
coordinates, and the four dropped coordinates cannot meet them all.
If a new point has support at most seven, parity in an old octad is
always adjustable outside the new support, so the sum of its four
largest magnitudes is ≤4. With all magnitudes ≤2, this and squared norm
8 force the root magnitudes (2,2). Indeed, writing a for the fourth
largest, convexity bounds the squared norm by
`4+(2−2a)^2+5a²` for 0≤a≤2/3, or
`(4−3a)^2+6a²` for 2/3≤a≤1; equality 8 is possible only at a=0
with magnitudes 2,2. Thus every genuinely new such point has a unique
heavy signed coordinate of magnitude >2. A used signed head requires
at least 32 root deletions (19 other positions, at least 13 with
both signs because they are zero); two heads share at most one deleted
root. At most 20 new points can share a signed head, because their
tails in dimension 19 have strictly negative pairwise dots. Hence
with at most 25 used signed heads, a sparse trade cannot gain even
before any dense conflicts are charged. This does not bound dense or
many-head trades.

## A sharp cluster with one matching class per block

In the [five-tetrad triality framework](../dim20_leech_tetrad_triality_20261126/),
fix one of the three four-axis root alphabets in each of the five
tetrads. Use all 16 signings of each transversal axis word (four
occupied tetrads, one omission). All requirements reduce to at most
two coincident occupied axes between distinct words: with signs
allowed independently, every coincident axis contributes 2, and
distinct axes in a fixed alphabet are orthogonal.

In a fixed omission class every projection onto three of the occupied
tetrads is injective. Two different omission classes have exactly
three common positions and their projected banks are disjoint.
For their sizes n_i,n_j this gives n_i+n_j≤4³=64 for each pair.
Summing yields total ≤160 axis words, or **2,560 signed
transversals** in this cluster.

The upper bound is attained. Use F₄ with its standard polynomial
basis; take the five columns (1,t,t²), t∈F₄, and (0,0,1) in F₄³.
Any three columns are independent. Put u=(0,1,0). For omission i
keep a four-coordinate tuple when its image under the four-column map
lies on either of two specified affine lines parallel to c_i.
Take the following basepoints for the two lines in each class, using
the integers 0,1,2,3 for 0,1,a,a+1:

| i | first basepoint | second basepoint |
|---|---|---|
| 0 | (0,0,0) | (0,1,0) |
| 1 | (0,0,2) | (0,1,2) |
| 2 | (0,0,3) | (0,1,3) |
| 3 | (0,0,1) | (0,1,1) |
| 4 | (0,2,0) | (0,3,0) |

The ten affine lines are disjoint. Within one omission, the two
parallel lines lie in one affine plane directed by c_i,u. Its
direction does not contain any other c_j, so changing one coordinate
cannot take an admitted tuple to another. Between omissions,
coincidence of all three common axes would give a point in the
intersection of the corresponding affine lines. Each line has
16 preimages in its omission (rank three); hence each bank has 32 and
the total is 160. This construction uses the actual overlap-three
mechanism not permitted by the old distinct pair-receipt design.
It is a standalone cluster, not a claim that several clusters would
be compatible with the record's transversals or dense stock.

## Replay and next step

Run `python verify_one_sign.py` (NumPy; exact integer operations) and
`python verify_cluster.py` (standard library, exact GF(4) arithmetic).
The most positive direction now is to fit several clusters with
different matching alphabets while jointly paying the dense loss.
A new transversal support missing one tetrad and not itself a Golay
octad conflicts with exactly eight dense signs in the full 2,048 bank:
the projected dense code on that eight-set is onto. Its dual is
the shortened Golay code plus the coset of a selection of one
coordinate per tetrad. That extra coset is odd in all five blocks
and cannot be contained in such a support; the only possible nonzero
shortened word contained in the support would be the support itself.
Thus a positive construction must share these blocker cylinders,
not count individual new transversals against an unchanged dense bank.
