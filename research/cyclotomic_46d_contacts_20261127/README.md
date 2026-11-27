# Lower83: true cyclotomic contact shells in 46D

No verified gain; live 46D bound remains 13,387,797 (P48n rank-40/E8-complement body with 6D tags, repository verifications/improved/dim41-46-p48-heads). The older published P48 comparison is 12,309,600, historically Conway–Sloane, with counts in [Dorofeev–Sun–Wang](https://arxiv.org/html/2607.20359v4). This is not a reoptimization of the 47D coarse distance menu. The object here is the **actual** 46D orthogonal projection, with no head coordinate.

## Geometric model

Use the order-47 fixed plane in P48q from [Nebe](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/DIM48.pdf), and its raw ternary-coordinate realization from the independently checked Lower81 census. In lattice units the two normal vectors are u=(1/2,…,1/2,−5/2)/√3 and v=(1/2,…,1/2,7/2)/√3. Their Gram is (6,1;1,8). The rank-46 body projection deletes infinity and subtracts the mean of the 47 finite coordinates. Using doubled raw integral coordinates x (actual lattice x/(2√3)), body point angles are rational dot tests with integer Gram
  G(x,y) = 47∑_(i finite) x_i y_i − (∑ x_i)(∑ y_i).
In particular the test is G≤0 or 4G²≤G(x,x)G(y,y).

The two largest outer radius groups of the norm-eight shell are A=(s,t)=(0,4) of population 7,945,350 and B=(2,−3) of population 6,777,870. Their actual projected squared radii are a=280/47 and b=278/47, both below six; because projection differences within one fixed profile are in the rank-46 lattice of minimum six, each group is a valid kissing configuration by itself. Write −A,−B for their negations. Any two of these four groups have a combined cardinality above the incumbent (even B∪−B gives 13,555,740), so a genuine two-bank trade, not one more small group, is the meaningful scale.

## Exact disagreement graphs, not lower bounds based on shortest differences

For any two groups of projected squared radii a_i,a_j, the angular failure test on their actual pair is
  ||x−y||² < a_i+a_j−√(a_i a_j) + Q(Δs,Δt)/47,
where Q(s,t)=8s²−2st+6t² and x,y on the left are unprojected vectors in lattice units. To get lower bounds on the attainable norm of an actual difference, subtract the indicated lattice normal; the body is unchanged, and either the parent lattice minimum applies or the projected points would coincide. The latter cannot occur because the four residues 8s−t mod47 are distinct. The following exhaustive obstruction shells result:

| Full join | Difference profile | Normal subtracted | Squared norm lost | Exact forbidden difference norm |
|---|---|---|---:|---:|
| A,−A | (0,8) | v | 8 | 14 |
| A,B | (−2,7) | v | 6 | 12 |
| A,−B | (2,1) | 0 | 0 | 6 |
| B,−B | (4,−6) | u−v | 8 | 14 |

Each threshold lies strictly between the stated even norm and the next larger even norm. The lattice minimum after translation excludes every smaller even norm in the row. The other two pair types follow by global negation. In particular any partially selected subfamilies in these joins have a direct exact compatibility check by *one* translated shell per pair type.

All four types really occur, including the mixed A,B case: the companion standard-library script checks two literal actual lattice points for each of the four rows, along with their exact kissing conflicts. To verify membership it builds the canonical ternary QR48 code via the 47 Paley circulant rows I+χ(j−i), infinity −1, checks self-duality and rank 24, then tests each literal witness. For doubled integer x, membership in the even 3-frame sublattice means x/2 has even coordinate sum and residue in the code. For odd x, test (x−(1,…,1,−5))/2 instead. The min-six even-neighbor identification is in Nebe/Lower81. The finite projected Gram triples for the four bad pairs are (1716,3360,3360), (2088,3360,3336), (1860,3360,3336), (1740,3336,3336); each first entry is the cross term, and each violates 4G²≤G_xx G_yy exactly. This excludes every complete two-bank join for these two giant profiles, not all large partial joins. It is independent of the coarse normal lower bound used in the preceding exploratory programme.

One particularly transparent real contact rule is the A,B row. Here the unscaled infinity difference is 4.5. Squared norm twelve requires squared raw finite difference 15.75; any such diff is a half-integer in every finite coordinate. Its low norm forces exactly **two** ordinary sites of magnitude 1.5 and all other finite magnitudes 0.5. Conversely these magnitudes give the unique offending shell. Thus for flat A, if B has three high ±1.5 sites (the usual B type), a true conflict must align at least one of those high sites with a nonzero A sign; every opposite alignment at a high site is forbidden for this lightweight difference. Retaining subsets by true high-site interactions, rather than rejecting a body residue wholesale, is an explicit discrete route. A wholesale deletion of an A bank for a B bank loses points.

## Actual contacts to the 9,372,270-point center

Every A donor has raw infinity 2, raw finite sum 10 and finite squared length 20. Integers of magnitude at least three among the finite coordinates would leave at most eleven unit magnitude weight contributors and violate minimum ternary codeweight 15 (the infinity residue is nonzero), so ordinary magnitudes are 0,1,2 and at most two sites can have magnitude two. If there is any such site, norm and finite sum force the spread max(x_i)−min(x_i)≥3. Explicitly with one magnitude two there are sixteen unit sites, and with two there are twelve: positive high sites require negative unit sites unless an opposite high is already present, and the analogous assertion holds if all high sites are negative. The central root √3(e_i−e_j) between maximizer and minimizer is an actual lattice minimum vector whose dot with the candidate is its spread, at least 3. But the mixed threshold with this center is √(6a)/2<3. Hence **every nonflat member of the A donor family is genuinely blocked by the A46 rail subset of the center**, not merely by a sufficient residue criterion.

The only remaining A candidates for retaining the full center are flat: 15 finite +1, five finite −1, 27 zeros and raw infinity 2. They are not automatically safe. Let K6 be the actual central minimum shell, and A6 the norm-six profile (0,4) (1,081 vectors by the prior exact census). A donor x of norm eight in A is center-compatible **if and only if**
  (x−K6)∩(A6∪A)=∅.
Indeed a central dot is an integer, at most four by the lattice minimum, and fails exactly at 3 or 4; subtraction then produces squared norm eight or six, respectively, without changing the profile. In particular any neighbour at squared distance six within A gives a literal center blocker from their difference. The verifier includes a flat A donor and its same-profile neighbour, yielding an actual central blocker; this demonstrates the second mechanism independently of the rails. No claim is made that every flat A donor is blocked, or that either of the large fibres has been optimally thinned.

The next quantitative target would be an explicit large independent pair of subsets using these actual single-shell contacts, or an exact packing/matching obstruction at the threshold: A∪B permits deletion of at most 1,335,422 points for a strict record; B∪−B only 167,942. The useful constructive move is therefore a high-retention block partition of the true bipartite graph, not another enumeration of whole-body residue angles. All support tables and literal validation in this note are integer arithmetic; computation was hosted only.

### Explicit completion of the 1,081-rail cut

A6 is one full orbit of the affine coordinate group i↦a i+b (mod 47), where a is a nonzero square and infinity is fixed. Such coordinate permutations preserve the QR code and the fixed normal plane. The updated replay supplies a literal A6 word (12 ordinary positive unit coordinates, 2 negative, infinity 2), reconstructs all 1,081 distinct orbit points, and independently checks their norms, membership and label. Completeness follows from the old exact shell census. It also verifies that the literal flat A donor suffers exactly ten true K6 blockers from subtraction of A6, separate from its same-A blocker. Thus the first half of the actual central-admissibility condition can be coded by exactly 1,081 simple signed linear inequalities: its finite raw inner product with every member of this explicit orbit must be ≤7 (the norm difference six event is at eight). This is precise enough for future optimization without any unknown lattice cataloguing. Not a proof every flat donor fails that cut.
