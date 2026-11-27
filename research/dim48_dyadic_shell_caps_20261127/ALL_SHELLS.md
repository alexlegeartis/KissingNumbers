# Addendum: no ambient-minimum affine dyadic shell beats 52,416,000 in 48D

Use the definitions, theta character proof, and 8, 10, and 12 results in the companion `README.md`. This addendum treats every squared norm m≥14. Thus the complete affine mechanism is excluded: for every extremal even unimodular L of dimension 48, every 2L⊂K⊂L of minimum at least m, and every affine coset a+K, the norm-m vectors of that coset number at most 52,416,000. This does not bound shell-sensitive affine selections with min(K)<m, multi-coset nonlinear ambient mechanisms beyond the ranges stated in the README, or arbitrary kissing configurations. All norms mean squared norm. No new kissing lower bound.

## 1. More coefficients in the same cusp table

Here again s is the minimum norm in a nonzero class t, and a=m12(t) in the s=12 row. The identical finite rational systems (no new modular assumption) supply two more columns:

| s | F14(t) | F16(t) |
|--:|--:|--:|
| 6 | −31471524864 | −1453140567648 |
| 8 | −1271878656 | 64577837472 |
| 10 | 973644800 | −22098568800 |
| 12 | −892032000−3956736a | 19203856800+95158272a |
| 14 | 902734848 | −19443574368 |
| 16 | −892032000 | 19203856800 |

The untwisted values are N14=15,173,208,925,056,000 and N16=327,259,384,199,748,000. For s=12, a≤96: if x and y give distinct lines of norm twelve in the same residue, (x±y)/2 are nonzero, so 6±(x,y)/2≥6 and (x,y)=0. At most 48 mutually orthogonal lines occur.

## 2. Spectral density bounds

For minimum fourteen use H14=F6+(185/1758)F8+(25/3516)F10+(1/9376)F12. Its degree is 28668906630000/293 and its minimum on nonzero characters is at least −70800/293. Check each table row, and for s=12 both endpoints a=2 and a=96. The ratio formula in the README gives inverse density at least 23890755584/59; in particular an affine dyadic kernel has index at least 2^29.

For minimum sixteen use H16=F6+(128/601)F8+(1750/48681)F10+(8/1803)F12+(11/97362)F14. Its degree is 2345624805888000/601 and its minimum on nonzero characters is at least −102400/1803, by the same endpoint table. This gives inverse density at least 2^36 even for nonlinear 2L-periodic sets. Sharpness in P48n follows from the Dursthoff ζ8 in §7.5 of his [thesis](https://www.math.rwth-aachen.de/~David.Dursthoff/DursthoffThesis.pdf): use (1+ζ8)(1+ζ8²)L, the √2 dilation of the norm-8 construction. The index is 2^36, and 2=(1−ζ8)(1+ζ8)(1+ζ8²) ensures 2L containment.

## 3. Reciprocal affine coset bound

Let W=K/2L, U=W^⊥, M=|U|, P=|V|=2^48. For m=14 or 16 the full Fourier table has the exact form
  F_m(t)=β+δ ε_t+sum_{j=6,8,10,12} c_j m_j(t)   (t≠0).
The coefficients are:
  m14: β=5351424, δ=−897383424, c=(−16187129856,−189923328,35454976,−3956736).
  m16: β=−119858784, δ=19323715584, c=(−716848496640,22686990336,−1327497216,95158272).
This follows by comparing the six rows: m6,m8,m10 are two exactly at their own types and zero otherwise, m12 is a at type12 and zero elsewhere. In terms of J(t)=sum c_j F_j(t) and Q_W(a)=|W|^−1 sum_{x∈a+W}ε_x, affine Fourier inversion, its reciprocal, and the plus-type Gauss identity give exactly
  n_m(a)=(N_m−β−δ)/M + β 1_{a∈W} + (δ 2^24/M) Q_W(a)
          + P^−1 sum_{x∈a+W}J(x).
For example, the reciprocity follows by summing F_j over a+W, which yields |W| sum_{u∈U}m_j(u)e_u(a). This also proves the identity without any assumptions on sign or frequencies of short types in the affine complement.

For m=14 and t≠0, J(t)/P≤303/500 (its maximum is at s=14). The s=12 row is constant as its a-slope vanishes. The exceptional x=0 contribution together with β is negative: β+J(0)/P<0. If M=2^29, b restricted to W must have positive rank. Otherwise q on W would be linear, and W∩{q=0} is a subspace of index ≤2 in W whose nonzero preimage has minimum 16 (from minimum 14 and squared norms divisible by four). Its index in L would be at most 2^30, contradicting the spectral min-16 threshold. Positive rank bounds |Q_W(a)|≤1/2 for all affine translates. Hence |(2^24/M)Q_W(a)|≤1/64. For M≥2^30 the same bound holds by the trivial |Q_W(a)|≤1. Using these estimates gives, uniformly,
  n_14(a) ≤ (N14−β−δ)/2^29 + |δ|/64 + (303/500)2^19 < 43,000,000.

For m=16 and t≠0, J(t)/P≤179428077/2097152, attained at type16, and the s=12 slope is −24. The zero-class correction β+J(0)/P is negative. Here M≥2^36, |Q_W(a)|≤1, hence the same formula gives
  n_16(a) ≤ (N16−β−δ)/2^36 + |δ|/4096 + (179428077/2097152)·4096 = 9,830,400.
In the origin coset the zero term replaces one of those copies of the positive upper bound by β+J(0)/P, lowering it by exactly 2400 to **9,828,000**, sharp by the (1+ζ8)(1+ζ8²)L construction and its common exact eigenspace (or the scaled norm-8 count in the main proof).

## 4. All remaining shells

Every nonzero class of L/2L has minimum ≤16 by the two-cusp table. Therefore any 2L-periodic subset with ambient minimum ≥18 has at most one residue, and in the affine setting necessarily K=2L. The exact counts are read directly at the zero cusp. For norm 18 only the odd types contribute, the counts are 36,848, 34,300, 34,928 for s=6,10,14. Norm 20 has even counts 385,024 (s=8), 393,216+228a (s=12), 393,216 (s=16), giving at most 415,104. Norm 22 has odd counts 3,357,774, 3,537,800, 3,527,064. Norm 24 has even nonzero coset counts 26,300,824, 26,099,200−896a, and 26,099,200. Its zero coset has exactly 52,416,000 as 2L doubles the original minimal vectors. The zero coset does not contribute to lower norms. A nonempty 2L-periodic set cannot have minimum >24 because differences within a single residue include the shortest points of 2L. Finally, the entire norm-6 shell already has exactly 52,416,000, and norms below six are empty (apart from the zero vector, which gives no direction).

Independent elementary check for 18,20,22: a coset of 2L is antipodal, and distinct off-antipodal raw dots have |t|≤m−12 and t≡m mod4. These give at most three absolute dot values: {2,6}, {0,4,8}, {2,6,10}. Homogeneous projective interpolation with degree at most six gives at most binom(53,6) lines, hence fewer than 52,416,000 points. The cusp counts are much smaller and exactly checked.

Run the dependency-free `verify-extension.mjs` to replay every numeric identity in the addendum. Together with the companion theorem this rules out the full *ambient-minimum affine dyadic whole-shell* strategy. The exact remaining place for constructive progress is a shell-specific exception: a shorter lattice difference in K is acceptable if no two selected shell vectors realize it, or if phases do not have a full affine dyadic period.
