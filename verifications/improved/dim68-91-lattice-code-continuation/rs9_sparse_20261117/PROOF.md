# Extended Reed–Solomon support-ten augmentation

Use n coordinates, n=87,...,91. Partition the first81 coordinates into nine groups of nine, indexed by F9. A last group consists of k=min(9,n−81) coordinates. In compressed QR48 constructions, coordinates0,...,83 are the42 triangle blocks and84,...,n−1 are head coordinates; thus ALL head coordinates lie in the last group. For90 and91 use the45 triangle blocks and, in91, one unused zero coordinate.

For each polynomial f over F9 of degree≤5 whose leading coefficient a5 belongs to the selected k-element subset, take one coordinate labeled f(t) in each finite group t∈F9, and the coordinate labeled a5 in the last group. There are exactly k*9^5 distinct supports of size10. If two leading coefficients differ, the difference polynomial has at most5 finite roots and infinity does not agree. If they agree, degree of the difference is≤4, so at most4 finite positions plus infinity agree. Hence distinct supports intersect in at most5 positions.

For each support use all64 signs of the binary linear[10,6,≥3] code with parity-check columns the nonzero4-bit vectors except{1,2,4,8,15}. They are distinct and nonzero, so no nonzero word has weight1 or2. Exhaustive1024-syndrome checking verifies dimension6. Their sum is zero, so the all-ones word belongs and the signs form an antipodal set. Assign each nonzero real entry magnitude sqrt(45/10). All points have squared norm45. Same-support inner products are≤(10−6)*45/10=18; different-support inner products are≤5*45/10=22.5.

Add all2n(n−1) scaled D_n roots with two entries ±sqrt(45/2). Their mutual inner products are≤22.5, and their inner products against support-ten points are≤2*sqrt(45/2)*sqrt(45/10)=45/sqrt5<22.5.

Dense QR48 body coordinates have absolute value≤1; each support includes at mostone head coordinate, whose absolute value is≤sqrt(3/2). Thus dense/support-ten inner products are≤(9+sqrt(3/2))*sqrt(45/10)<22.5. Squaring reduces this to sqrt(3/2)<5/3. Dense/root inner products are≤sqrt135<22.5. Dense points have at least42 nonzero coordinates, so cannot duplicate support-ten or root points.

Dense configurations are the previously independently audited9,12,16-head constructions in87,88,89D, and shortened45-triangle code in90,91D. Their counts are respectively2*M*3^18 for M=9,12,16,27,27. This new layer REPLACES every previous sparse layer; it is not a union with them. Independent audit of this new extension is pending.
