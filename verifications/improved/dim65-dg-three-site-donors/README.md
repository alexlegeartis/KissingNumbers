# Three-site Delsarte–Goethals replacements in 65 dimensions

This construction gives **424,169,298 contacts**, an exact increase of **1,310,720** over the previous repository bound (422,858,578). It replaces selected dense words by triples instead of pairs. As in the existing construction, the 65-dimensional four-site layer uses the Bao–Ji existence theorem.

## Construction

Start with the full 335,544,320-word dense carrier of the [pinned retained-atlas result](../dim65-68-dg-retained-atlas/). Keep its 40,414 rephased zero-tail support banks, its 10,741 four-site supports with all 16 signs, and the two tail poles. Keep all 44 previously retained hyperplanes, with both head types per hyperplane (180,224 contacts), and omit all ordinary poles. The older pair replacements are not used.

Partition 63 body coordinates into 21 disjoint triples. Assign each to a distinct cubic coset using `triple_assignment.txt` (truth-table ordering). In each assigned coset, designate the full 131,072-word Delsarte–Goethals subfamily as owners. Replace each owner with three vectors. At its designated triple, multiply its signs by the first three entries of one row below; keep the other 61 body signs unchanged and use the fourth entry as the signed tail:

| First | Second | Third | Tail |
|---:|---:|---:|---:|
| −1129/200 | 5/4 | 5/4 | 0 |
| 7/5 | −11/2 | 32/25 | 107/100 |
| 7/5 | 32/25 | −11/2 | −107/100 |

Normalize after applying the profiles. Their squared norms are 3839641/40000 and 959933/10000 (twice). We remove 2,752,512 owners and insert 8,257,536 donors. The resulting count is 418,484,050 + 180,224 + 5,505,024.

## Certificate

The parent code distance is at least 16, and within a Delsarte–Goethals family it is at least 24. Use those distance bounds for different owners. The designated triples are disjoint across cosets; the two signed tail products are taken explicitly. Same-owner pairs are computed directly. Every owner has even parity on every retained weight-16 support. The donor flips exactly one body sign, creating a forced disagreement on any retained support that contains that site. Any four-site vector or tail pole is bounded directly using its support size and the new profile magnitudes. For the retained 44 hyperplanes, the 60,555,264 cubic-vs-quadratic words under the assigned triples all have even Hamming weight at least four. In each of the 392 weight-four words, none of its disagreements is a marked coordinate. Flipping the one negative multiplier can therefore only help at distance four. The other cases use distance at least six. The tail contribution is included explicitly. All gates are checked with exact arithmetic.

Run `python3 -B verify.py` with a C++20 compiler on a BMI2-capable x86-64 machine. The main verifier checks the hashes and reruns the inherited certificate, including 2,330,703,675 sparse support-pair comparisons and 31,512 parity guards. It then checks the new gates and count. `verify_gates.py` is a separate integer implementation of the new gates; it does not replace the inherited replay. The pinned predecessor files are available in commit `06943cf8d8237fac6e1f7ddb968b8e494618f00c`.
