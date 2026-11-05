# Exact ceiling for signed four-support equators in 12 coordinates

Let E be any set of vectors with exactly four entries in {−2,2}, all other entries zero, and pairwise inner products at most 8. Then |E|≤816. This is an upper bound only for this specified family.

For each sign pattern σ∈{−1,1}^12, collect supports of vectors in E whose four signs agree with σ. Their supports are distinct four-subsets and no two can share a triple: such a pair would have dot product 12. Thus the supports form a 3-(12,4,1) packing.

Any such packing has at most 51 blocks. Let L be its uncovered triples. Every pair {i,j} belongs to 10 triples in all and every block containing the pair covers two of them, so each pair has even degree in L. For each point i, |L_i|=55−3r_i is congruent to 1 modulo 3. The link graph of i in L has all degrees even. It cannot have just one edge; if nonempty it has at least three edges, and congruence forces at least four. Therefore sum_i |L_i|≥48, hence |L|≥16. As |L|=220−4b, b≤51. Moreover r_i=(55−|L_i|)/3≤17.

Each signed four-support vector agrees with exactly 2^8 global sign patterns. Summing the 51-block bound over 2^12 sign patterns yields 2^8|E|≤2^12·51, so |E|≤816. The 816-point equator of the saved source attains this. Restricting to sign patterns with a fixed signed coordinate and using r_i≤17 similarly gives at most 136 members of E containing that signed coordinate. For a fixed signed pair, every remaining signed coordinate occurs in at most two members (the fourth coordinate is forced, with at most two signs); counting those 20 signed coordinates twice per vector gives at most 20 members.

Since every coordinate axis ±4e_i is compatible with every vector in E, a 12D configuration consisting only of E and all 24 axes has at most 840 points; it cannot improve the known 841-point 12D construction. This says nothing about 12D constructions with additional vector shapes or changed axes.
