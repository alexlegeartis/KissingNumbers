# Three-site donors with a 46-hyperplane atlas

| Dimension | Preceding bound | New lower bound | Gain |
|---:|---:|---:|---:|
| 65 | 424,169,298 | **424,177,492** | 8,194 |
| 66 | 432,812,452 | **434,139,558** | 1,327,106 |
| 67 | 451,803,398 | **453,138,696** | 1,335,298 |
| 68 | 481,928,984 | **483,338,010** | 1,409,026 |

The 65-dimensional four-site stock still uses the Bao–Ji existence theorem. This release does not assert a worldwide priority comparison beyond the named older bounds.

## Dense donors

Use the 21 disjoint triples and distinct cubic labels in `triple_assignment.txt`. In each listed cubic coset remove its 131,072 Delsarte–Goethals owners and replace each by three weighted sign vectors. The 61 unmarked body signs remain unchanged and the designated triple multiplies its three signs by a row's first three entries. In 65D use the profiles of the [previous three-site construction](../dim65-dg-three-site-donors/). In 66–68D use the following numerators over 100,000; copy the fourth entry to each of the first two tail axes and put zero in other tail coordinates:

| First | Second | Third | Tail coordinate (twice) |
|---:|---:|---:|---:|
| −558229 | 105952 | 154047 | 40547 |
| 160809 | −557420 | 100000 | 40547 |
| 103127 | 157232 | −557870 | −40547 |

Normalize each donor after constructing it. Same-owner distances use the literal profiles. For different owners in the same triple, the unweighted Delsarte–Goethals distance is at least 24; across triples and between donor owners and the rest of the dense stock, the parent distance is at least 16. The included checker verifies the exact angular bounds for all cases. Restore the two ordinary poles on the one unused truth-table coordinate (index 7), whose donor magnitudes are one, and retain the tail poles.

## Tagged sparse banks

The old signed stock is in the [retained-atlas certificate](../dim65-68-dg-retained-atlas/). Each zero-tail support has 16 ordinary sites and even parity; a donor that flips a sign on a site in the support forces a mismatch because the parent word has even parity on that support. If a weight-16 support uses four tail coordinates or overlaps at most one nonzero tail entry of our donors, the exact cap holds without a parity comparison. For a 14-body bank using exactly the two donor tail axes, the only unsafe magnitude patterns hit the negatively weighted site and at least one other marked site. Both tail entries have the same sign. Every one of the 4,767 body masks of this type is orthogonal to the RM(2,6) basis, so the donor sign parity when the large mark is present is the assigned cubic word's restriction parity plus one. Opposite total parity forces a mismatch of magnitude at least 40547/100000. In the previous signing, 51 banks would fail to guarantee this mismatch. Flip their total parity, along with the 68 shared overlap-guard bits in `rephase.json`. The native checker verifies that all sparse–sparse guards remain balanced and consistent, and the exact tag check then finds no violations. No sparse banks are deleted.

## Atlas and replay

Keep all 44 old atlas hyperplanes and add rows 27 and 115 in the original zero-based indexing. All 46 hyperplanes keep their full available head stock. The replay constructs the rephased sparse duals first, preserving every original overlap guard, and then chooses balanced parity guards against the 46 heads from those rephased duals. It also checks the head-to-head angular gates. The old C++ proof is run separately and SHA-pinned, so the reconstructed donor families and historical counts are independently checked before applying the new delta. `cross.cpp` checks 241,172,480 body coset words, including the complete parent cubic palette, and confirms that all 4,344 weight-four words avoid the marked positions of the appropriate cubic (unassigned cubics have no marks). All reachable weights are even.

Run `python3 -B verify.py` with a C++20 compiler on a BMI2-capable x86-64 machine, with predecessor directories from this repository present. The exact scripts need no optimization or network. The supporting `verify_diagonal.py` can be run alone to check the new rational angular inequalities. The final dimensions differ only in their inherited stocks and number of allowed atlas heads per hyperplane.
