# Mixed-shell 71-dimensional kissing construction

Let L be an extremal even unimodular 72-dimensional lattice of minimum8 (e.g. Nebe's Gamma72). Fix x in L with squared norm8. For (s,j) in {(8,3),(8,4),(10,4),(10,5)}, take every y in L of squared norm s and inner product (x,y)=j. Replace y by its orthogonal projection p=y−jx/8, normalized to unit length. All points lie in x-perpendicular, dimension71.

The residual squared norm is r=s−j²/8>0. For two distinct original points y,z in slices(s,j),(t,k), minimum8 yields (y,z)≤(s+t−8)/2. Their projected inner product is therefore at most B=(s+t−8)/2−jk/8. The supplied standard-library verifier checks all10 pair types exactly: B≤sqrt(r_s r_t)/2. This also proves no projected normalized points collide. These inequalities include pairs within a slice (distinct original points only).

The slice populations respectively are13,959,168;127,800;2,595,532,800;13,959,168. Their sum is **2,623,578,936**. These follow from shell11design moments and the theta coefficients6,218,175,600 and15,281,788,354,560. The verifier checks every even moment throughdegree10. The norm10 values are also Proposition3 of Scott Duke Kominers, https://arxiv.org/html/2606.08734 . The bijection y↦x−y identifies the first and fourth populations and negates their projections; the second and third slices are each centrally symmetric under the same map.

Source premises remain explicit: Nebe's extremal lattice exists and has minimum8; the spherical11design theorem for its nonempty shells; the extremal theta series. No individual enumeration of the norm10 shell is required. The result is a construction lower bound, not an exhaustive worldwide-priority claim. Live repository comparison observed17November2026:2,603,658,750; gain19,920,186 (+0.765084%). Independent exact arithmetic and geometry checks accompany this construction.


For shell8, the only nonnegative inner products are0,1,2,3,4,8, with8 occurring only at x. For shell10 they are0 through5: larger values contradict minimum8 for y−x. Antipodal symmetry and the six even moments give a unique solution for the six unknown multiplicities (the norm8 endpoint8 is already known). Thus the moment checks determine counts, rather than merely test a proposed distribution.
