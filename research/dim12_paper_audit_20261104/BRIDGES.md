# An alternative saturated bridge coupling

This is a construction beyond the fixed bridge choice in
[arXiv:2606.18984](https://arxiv.org/abs/2606.18984). It does not refute the
paper's uniqueness statement: K6 has one factorization up to relabeling, but
six different factorizations on a fixed labeled vertex set.

Let F0,...,F5 be those six labeled factorizations. Each is a set of five
perfect matchings of K6. Any two different factorizations share exactly one
matching. `alternative_bridges.py` verifies these finite facts by exhaustive
enumeration from the fifteen perfect matchings.

For a left edge {u,v}, put Muv = Fu intersect Fv, its unique shared matching
on the right. For each of the three edges {k,l} in Muv, include the support
{u,v,6+k,6+l} with all sixteen sign choices, each coordinate of magnitude 1/2.
This gives 15*3*16=720 vectors. Each row neighborhood is a perfect matching,
and each column neighborhood is also a perfect matching.

Two different supports cannot intersect in three coordinates. If their left
edges agree, their right edges are disjoint. The same holds with left and
right exchanged. Otherwise the two left edges and the two right edges each
intersect in at most one vertex. Thus all resulting vectors have pairwise
inner product at most 1/2. Within one support, changing a sign also gives at
most 1/2. Norms are one.

Any of the paper's 60-point blocks remain compatible: a 48-system vector has
all coordinates bounded by 1/2, while an axis vector has one coordinate of
magnitude one. Their inner products with a bridge are at most 1/2. Consequently
two such blocks and the alternative bridges give a valid 840-point arrangement.

The 15-by-15 support-incidence graph is connected. The paper's coupling instead
has five K3,3 components. More intrinsically, the completed arrangement with
the paper's explicit Clebsch blocks has 351136 contact triangles, compared with
353696 for the exported start. Both have 34560 contact edges. Exact integer
arithmetic in Q(sqrt(2)) establishes these counts, so the point configurations
are not isometric.

There are exactly two combinatorial types under independent relabeling of the
left and right coordinates, within the saturated family of 45 distinct 2+2
supports carrying every sign pattern. Here is a proof of completeness for that
limited family.

At any left vertex u, the matchings Muv for its five incident edges must partition
all right edges: each right-edge column is itself a matching of left vertices.
They therefore form a factorization Fu. For any left edge uv, Muv belongs to
both Fu and Fv. If two vertices have the same factorization F and a third has
a different factorization G, the two edges to that third vertex must both use
the unique matching in F intersect G. That repeats a matching at the third
vertex, which is impossible. Hence either all six Fu agree, or all six differ.

If all agree, the left edges are partitioned into five matchings, giving the
paper's product coupling. There are 6*6*5!=4320 labeled product couplings.
If they all differ, their assignment to the six left vertices is a bijection,
giving the alternative construction in 6!=720 ways. `couplings.py` enumerates
all 5040 distinct matrices and pins their hashes. This classification does not
cover partial sign patterns or other support sizes.

The alternative is an exact 840-point seed, not a proof of 842. The campaign's
147 alternative-bridge 842-point trials did not reach the kissing threshold.
