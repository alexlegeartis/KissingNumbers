# 16D quaternionic packet matchings: ties, exact blocker prices, and the full-packet ceiling

This is a supplemental structural route under the same fixed complete 1,248-point axial/two-block scaffold as the [main 16D argument](README.md). It does not improve **4,320**. Credit for the three uniform Latin-cube layers belongs to [Rustem Takhanov](https://arxiv.org/html/2611.00927v1#S2). These constructions are not claimed to be a new set of isometry classes.

Split four quaternionic blocks into two pairs I=01 and J=23. In either pair, identify the four projective axes of each \(C\)-level \(t\) with \(\mathbb F_2^2\). For each \(h\in\mathbb F_2^2\), let \(L_{t,h}\) be all 16 signed two-letter rows in level \(t\) in both slots, with XOR of their two axis labels equal to \(h\). The distinct packets form twelve full orthoplexes (norm squared two); define \(R_{t,h}\) independently on the other pair. All cross products between distinct packets on one side are at most one: if the two levels are equal, at most one axis agrees; otherwise both per-slot products are at most \(1/2\). Within a packet the distinct two-letter rows are mutually nonpositive.

For any perfect matching \(\pi\) of the twelve left and right packet labels, use all Cartesian products \(L_i\times R_{\pi(i)}\), with overall weight 1/2 on each quaternionic block. Within a packet Cartesian product, a repeated half can contribute two and the other half is nonpositive. Between packets the product sum is at most two. Thus there are \(12\cdot 16^2=3,072\) valid new points, or 4,320 with the full scaffold. Different matchings need not leave a full uniform-level Latin layer intact.

## Exact prices for individual signed replacements

With respect to one half of a candidate, call it pure if both C-levels coincide, mixed otherwise. The dot distribution over each 16-point half-packet is independent of axes and actual sign choices:

* For a pure half, the unique exact packet containing it has values \(2,0,-2\) with multiplicities \(1,14,1\); the other eleven have \(1,0,-1\) with multiplicities \(4,8,4\).
* For a mixed half, the four packets of each of its two involved levels (eight in all) have \(3/2,1/2,-1/2,-3/2\) with multiplicities \(1,7,7,1\). One of four axis pairings matches the involved level axis; it produces one of each signed inner product. The other twelve phase choices produce six of each \(±1/2\). The four packets of the third level have \(1,0,-1\) with multiplicities \(4,8,4\).

Summing strict pairs with combined dot above two gives:
* Pure/pure candidate not already present: **8**; the two exact half-packets run in different matched pairs and each contributes four blockers.
* Pure/mixed: **36** regardless of matching. If the exact pure packet meets an involved mixed level, its eight positives are blockers and there are seven other involved packets each contributing four. If it meets the unused level, it contributes four and there are eight other involved packets each contributing four.
* Mixed/mixed: **\(64-7x\)**, where \(x\) is the number of matched pairs joining involved levels on both sides. Here \(4\le x\le8\); each of the \(x\) such pairs supplies one blocker and each of the \(16-2x\) involved/unused pairs supplies four. Thus the intermediate prices **15,22,29** are possible (and precise), as well as 8/36. For any desired \(x\) in this range one can choose a matching joining exactly \(x\) of the eight designated left labels to the eight designated right labels.

This also proves literal maximality against single additions for *every* complete matching. The cheaper intermediate options may be useful only in a nontrivial simultaneous trade.

## No gain from any collection of complete orthoplex Cartesian packets

Allow arbitrary 16-point norm-squared-two orthoplexes in \(C\times C\) for the two halves, not necessarily the XOR atoms or even pure levels. Take entire products across matched half-packets. For two distinct half-orthoplexes in \(\mathbb R^8\), the greatest cross product is at least \(2/\sqrt8\) by orthonormal-basis completeness and antipodes (the sum of squared normalized projections of any given unit row on an orthonormal basis is one). All possible products are half-integers, so this greatest product is at least one. Because any left and right crossing pair of rows can be simultaneously realized in the full Cartesian products, the two half-maxima must sum to at most two. Thus each maximum is exactly one and the half-packets in a common side form a mutually compatible union with distinct rows.

Every arbitrary kissing code in the equal-block \(C\times C\) half-alphabet has at most \(24\cdot8=192\) rows: fix one exact signed C-letter on the first block; its second letters must be pairwise nonpositive in \(\mathbb R^4\), hence at most eight. Thus at most twelve distinct 16-vertex packets can be used; a thirteenth full packet is impossible, regardless of mixed-frame design. Even without a packet condition, use of **only pure halves on one fixed pair of positions** caps the full 16D four-block stock at 3,072: that half ranges over just \(3\cdot8^2=192\) exact signed pairs, and for a fixed half the other half is a nonpositive code of size at most 16 in \(\mathbb R^8\).

So a genuinely different gain here needs thin/signed partial packets, or to modify the fixed scaffold. This closes only the full-packet and pure-half scopes.
