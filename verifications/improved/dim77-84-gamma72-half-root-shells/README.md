# Exact kissing counts for half-root Γ72 attachments and orthogonal sidecars

For 77–80 dimensions place four vertical coordinates on D4. Attach a tail z chosen from the first one, two, three or four simple roots of that D4, with combined vertical squared norm
```
3 ||z||² / 4 + ||x + z/2||²,  x in D4.
```
The corresponding horizontal vectors are represented by eight literal triples in the integral Leech block model of Γ72.

For dimensions 81–96 add up to four further orthogonal D4 tails, each with its own horizontal phase map. Prefixes of the simple roots give the intermediate dimensions. Their A-label planes are chosen disjointly from the earlier A span.

The verifier independently reconstructs all relevant vertical patterns and their multiplicities. For the half-root block, a vertical squared norm at most 8 forces |z_i|<=3 and |x_i|<=4, and the sidecars are ordinary D4. It convolves their exact histograms, adding horizontal syndromes by XOR. Then it verifies the clearance and counts exact horizontal vectors at the required square length. In particular, for horizontal norm five the two A labels are independent, so the raw Leech-block norm ten consists of one norm-four and one norm-six block. Each compatible pairing gives four signed choices. At the other needed norms the verifier uses the exact pinned parent shell rules.

Replay with `python verify.py ../dim73-82-gamma72-depth6-lattice`; NumPy and SymPy are needed. Both frozen count tables are compared to independent replay. A certified count in a row is not necessarily larger than other constructions in that dimension.

## Published comparison

The [comparison rows](comparisons.json) compare to repository commit [`691f6183fa1c`](https://github.com/alexlegeartis/KissingNumbers/tree/691f6183fa1c24a114b055ac183f8c5ddf9147a8) and to the named derived baselines in its table. Above 84D these recorded examples do not improve that live comparison. No optimality or worldwide priority is asserted.

| Dim. | Verified lower bound | Gain over preceding repository bound |
|---:|---:|---:|
| 77 | 6,222,013,132 | 275,906 |
| 78 | 6,223,119,444 | 827,718 |
| 79 | 6,224,780,712 | 1,656,324 |
| 80 | 6,228,104,952 | 3,314,328 |
| 81 | 6,228,378,586 | 3,314,392 |
| 82 | 6,228,925,854 | 3,314,392 |
| 83 | 6,229,747,428 | 2,292,250 |
| 84 | 6,231,391,584 | 132,368 |
