# 37D: one complete QR mother center cannot be profitably replaced by the natural five rails

The old 544,666-point code is pinned in main `5fe68ce190db9760cbf36ee0775a59d775d48730`. Source and quotient instructions are in the companion flat-body graded QR note. Remove all four 4,096-point grades of zero-based mother center 7 (16,384 rays). Keep the other ten centers and all nondense layers. Try to refill with arbitrary unions of complete translates of H0, the rank-12 projected zero-grade QR kernel, on all 36 ordinary sites with amplitudes ±2. Use one constant scalar head from {0,+1,+4,−1,−4} in each new bank.

One exact quotient BFS from the ten remaining grade-zero centers gives the distances in each of the four grades by XOR translates as in the main flat-body proof. The exact allowed distance lists, in grade order (+4,+1,−4,−1), are:
```
 0 :  9  9  9  9
+1 :  9 10  8  9
+4 : 10  9  6  8
−1 :  8  9  9 10
−4 :  6  8 10  9
```
These are simply the first integer d satisfying `4(144−8d+zm)^2 <= (144+z²)(144+m²)` whenever raw dot is positive; the raw product decreases in d. The candidates' conflicts with existing nondense norm-32 vectors and poles are ruled out by the published ±1/±4 heads and direct bounds for head 0 (raw at most 32 with any octad, 24 with any radical, 16 with any root; threshold sqrt(144*32)/2). Distances within a complete new bank are ≥12 from the unchanged H0. The census finds 4 candidates at head 0, 39 at each nonzero head, exactly 160.

Candidate-to-candidate exact minimum Hamming distance is calculated over all 4,096 H0 differences, then the two raw head amplitudes are included in the same signed dot gate. The maximum clique in this compatibility graph is **4**, verified by the simple exhaustive include-first/exclude-first recursion: an include branch intersects the remaining candidates with neighbors, while the exclude branch leaves them; a branch with remaining size plus chosen size no larger than the incumbent can be skipped. The loop exits early only if five are found, which the complete run does not do. There is no sampling, solver tolerance or inferred clique upper. E.g. a literal 4-clique consists of `(16767921,−1),(16713146,+1),(16443056,+4),(16120108,−4)` (scalar head and decimal mask representative as in the four-bank tie proof). Therefore this replacement class ties the old score and cannot improve it; the result says nothing about non-full banks, other magnitudes, different codes or multiple deleted centers. The individual arbitrary real-head obstruction for the *intact* old dense code is a separate stronger statement with a different scope.