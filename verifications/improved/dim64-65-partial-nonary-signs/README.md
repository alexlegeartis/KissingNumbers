# Partial signed nonary layers in dimensions 64 and 65

The earlier 40-point nonary fibers were selected using an unsigned intersection
cap of four. On some shared five-coordinate sets, the existing sign code
misses several sign patterns. These missing patterns allow partial fibers
whose support overlaps may exceed the unsigned limit. No point in the
preceding constructions is removed or moved.

| Dimension | Preceding bound | New bound | Gain |
|---:|---:|---:|---:|
| 64 | 335,789,656 | **335,795,354** | +5,698 |
| 65 | 335,857,882 | **335,863,580** | +5,698 |

In the inherited squared norm 64 convention, each new signed point has nine
entries of magnitude 8/3. It is encoded by a 64-bit support and by a submask
indicating negative entries. If two supports overlap in k coordinates and
their signs differ on e of those, the inner product is (64/9)(k-2e).
The verifier checks k-2e ≤ 4 against every point of every inherited nonary
fiber where k≥5, and between all new points (on the same support it checks
the equivalent sign distance bound). Inherited nonary signs are generated
directly from the punctured Gray image of the pentacode in the preceding
[nonary proof](../dim61-66-nonary-affine-caps/).

All 771 new supports avoid the affine Steiner quadruple blocks because their
36 pair-XORs are distinct, and every support intersects each of the 32,801
inherited weight-16 supports at most six times. Consequently the products
against the inherited four- and sixteen-coordinate layers are at most
(8/3)·4·3 and (8/3)·2·6, both equal to 32. Dense ±1 and axis ±8
products are at most 24 and 64/3. This proves the 64D augmentation. In 65D
the very same points are equatorial. The parent proof establishes that
all lifted bodies are dense or GF(4)-line owners from the 64D parent,
contracted by √3/2, and pole products vanish. Distinct new points are checked
by distinct mask pairs; the original nonary supports are explicitly excluded,
and the other old equatorial support sizes differ.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` in this directory with NumPy
installed. It pins all earlier inputs and does not need the discovery pool or
solver. The inherited configurations' own validity and counts are the
preceding certificates' claims, not replaced by this replay. The claims are
ordinary kissing-number lower bounds, not optimality statements.
