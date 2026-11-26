# Exact non-flatness gates for extending the 40D QR split

This is a scoped upper bound on **augmentation** of the fixed 589,823 QR sign vectors and 40 daughters in the [published 1,139,159-point construction](../../verifications/improved/dim40-qr-central-split/). It does not require any of the Smart supports or the roots. It is not a new kissing-number bound or a statement about changing the retained QR vectors. The allowed code threshold is 1/2.

Put θ=2 arctan(9/100). Every unit vector within angular distance θ of any of the 2^40 flat signs is strictly blocked by the retained bank. Thus an additional vector must obey Σ|x_i|<√40·9919/10081. On the 28 deepest QR quotient classes the stronger gate is θ=2 arctan(3/16). In particular neither flat signs nor small perturbations of the deepest Hamming holes can extend this bank.

## Exact syndrome geometry

Project the zero-prefix shortened QR48 code to the 40 body coordinates. It is linear of dimension 16, and the published 589,824 original dense signs are nine fibres of it. The 24-dimensional syndrome map and nine fibre centres are regenerated from the same QR mother code, not accepted as extra data. Reachability in the XOR graph of the forty coordinate syndromes gives 1,414, 28 and 0 uncovered quotient elements at Hamming radii seven, eight and nine, respectively.

A radius-seven covered flat has a blocker at dot at least 13/20 unless the only usable original row was the deleted all-plus. For a flat at distances 1 to 7 from all-plus, take a daughter hot coordinate at a disagreement. Its dot is at least (27+a)/√(40D)>649/1000, where a=√41−2 and D=2(a+1)^2. At distance zero the daughter dot is √(41/80), again larger. With cos θ=9919/10081 and sin θ=1800/10081, a blocker with flat dot at least 649/1000 still has dot strictly above 1/2 throughout the angular neighborhood. The replay checks the relevant squared rational and surd comparisons.

For each of the remaining 1,414 quotient elements, collect 72 to 96 distinct error masks of weight at most nine from noncentral retained fibres. If the sample size is m, the sum of weights T and per-coordinate counts c_i, the average of these blockers has radial dot ≥11/20 and squared tangential norm exactly Σ(40c_i−T)^2/(16000m²)<1/25. These identities follow by averaging their normalized sign rows componentwise. Every flat in the same quotient class has the same profile by a common sign change, and all sampled points avoid the deleted row. For an arbitrary vector x within θ of the flat, the centroid gives
`x·mean > ((11/20)·9919 − (1/5)·1800)/10081 > 1/2`.
Therefore at least one actual retained point blocks x.

At each of the deepest 28 quotient elements, all distance-nine errors are in noncentral fibres. There are 397 to 424 of them; including them all makes the radial dot 11/20 and tangent norm <1/30. The extreme squared tangent ratio is 2,605,400/(16,000·405²). At θ=2 arctan(3/16) the strict centroid bound is 2653/5300>1/2.

## Replay

Run `node verify_qr40_gate.js`. The standalone code regenerates the syndrome columns from translated quadratic residues; computes the nine radius expansions by packed word permutations; and checks each near-face and deepest centroid using integer Hamming masks. It uses several MBytes of memory and a bounded finite search. No Smart code data or random sample is required.
