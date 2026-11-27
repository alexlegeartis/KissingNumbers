# 48D addendum: a parity-screen relaxation for the norm-eight shell still cannot gain

The full-period proof forbids every norm-six vector r in the period K. A weaker sufficient idea is to allow such vectors but retain only norm-eight candidates y with (y,r) even for every residual r. This is safe: if y,z have norm eight and y−z=r of norm six, then (y,r)=3 is odd. It does not require the full set of ambient vectors to have squared minimum eight. Here is an exact obstruction to this entire parity-screen strategy; no claim about arbitrary shell-sensitive conflict exclusions.

Let W be the period residue space in V=L/2L and R the span of the norm-six residues in W. All parity-screened norm-eight points in an affine W-coset, if any, lie in an affine coset of W'=W∩R^⊥. Every norm-six residue of W' lies in R and is orthogonal to W', so the residual roots in W' pair evenly. Distinct nonantipodal norm-six vectors from this residual root stock have absolute normalized dot either 0 or 1/3. For their projective lines, the even Gegenbauer polynomial
  1+(14053/839)G_2+(4228308/839)G_6
is zero at both dots, with positive coefficients in dimension 48. Consequently there are at most floor(4243200/839)=5,057 lines, or r<=10,114 oriented roots, in the period W'.

Let M=[V:W'], U=(W')^⊥, h=F6 and g=F8. As previously, min h>=−12,800 on nonzero characters, while sum_{u∈U}h(u)=M r, so r>=12800(4096/M−1). If M≤2048 this contradicts r<=10,114; hence M≥4096. Also
  sum_{u∈U,u≠0}(h(u)+12800)=12800(M−4096)+M r.
In Fourier inversion replace g by 304800 and error bounded by 154(h+12800). The bound for any affine W' coset is therefore its ordinary periodic dyadic upper estimate, at most 9,828,000 for M≥4096, plus 154r<=1,557,556. The screened set has no more than **11,385,556** points. Thus merely imposing additional parity checks to erase the odd positions of every residual short root cannot rescue this large shell. The truly different route would exploit exact dot 3 rather than every odd dot, or nonbinary phases.

The polynomial roots can be checked directly with normalized degree-2 and degree-6 Gegenbauer polynomials for \(\lambda=23\). The ratio is \(4243200/839\); \(154\cdot10114=1557556\).
