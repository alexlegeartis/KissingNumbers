# Simultaneous support rephasing in 62D and 63D

| Dimension | Verified lower bound | Gain over preceding bound |
|---:|---:|---:|
| 62 | 170,279,948 | 3,358,720 |
| 63 | 230,159,620 | 3,878,912 |

The [previous construction](../dim62-single-mark-square/) uses two types of signed banks, sixteen-site ordinary supports and fifteen-site supports with one fixed length-two head. We keep all its dense words, planes, and poles. We keep every previous support and the number of signs it contributes, but change the sign cosets on some of those supports. This allows new supports without paying a full bank for each conflicting overlap.

A new support has size sixteen on the original coordinate set. A support hitting one mark is orthogonal modulo two to every dense sign word, and its body coset uses the same compulsory total-parity condition as before. All other new supports avoid the marks. Any potentially conflicting overlap between a new and an old bank has an eight-site guard contained in the shared ordinary support. Set the new bank's guard parity to zero and the old bank's guard parity to one. The verified overlap exceeds the magnitude threshold by at most two ordinary hits, so the forced sign disagreement reduces the dot product by the required eight. All new–new and old–old pairs satisfy their ordinary magnitude thresholds outright. The dense, plane and pole bounds continue to hold.

Changing a coset does not reduce its size. For each bank the incidence patterns of all of its guards form an affine subspace of bits with equal multiplicities across sixteen sites, using one dummy site outside every guard on a fifteen-site bank. This assigns distinct four-bit column labels to the sites, making every guard an affine extended-Hamming parity check. The replay checks the joint coset parities explicitly. Thus each sixteen-site bank contains 2,048 signs and each fifteen-site bank 1,024, all with internal Hamming distance at least four. This is a simultaneous certificate, not a sum of independently priced single edits.

Run `python3 verify.py` with the same dependencies as the parent (Python and a C++17 compiler). It reconstructs the old supports, decompresses the two fixed witnesses, checks every potentially conflicting overlap and verifies the parity labeling and counts. These are lower bounds, not optimality claims.
