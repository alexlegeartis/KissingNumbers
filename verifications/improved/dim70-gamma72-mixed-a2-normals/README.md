# A mixed-normal Γ72 A2 section in dimension 70

The three compatible projected fibres have **610,937,568**, **610,937,568**,
and **139,276,672** vectors. This gives **1,361,151,808** kissing contacts,
**9,408** above the preceding bound of 1,361,142,400.

Both normals are norm-eight vectors of Nebe's Γ72 lattice and their inner
product is four. The first lies in one Leech block; the second uses all
three. This is not another choice of two roots within one Leech block.

## Exact corner count

Use Leech triples (a,b,c) with inner product half the sum of the block
inner products. The pinned parent data give its Leech Gram F, the two
polarisation maps A,B, and the integer map T from the Γ72 catalogue
basis. We check T diag(F,F,F) Tᵀ=2G, that A+B=I and AB=2I, and that the
binary membership conditions are

    aB = bB = cB   (mod 2)
    (a+b+c)A = 0   (mod 2).

The rows of T satisfy all the conditions. The complementary maps both
have binary rank twelve, so the congruence lattice has index 2^36, which
is also the absolute determinant of T. Thus the equations are sufficient
as well as necessary. The Γ72 Gram G is the table from the
[Nebe–Sloane catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html).
Its minimum and unimodularity are inherited from [Nebe's construction](https://arxiv.org/abs/1008.2862).

In the raw Leech coordinates the first normal is (0,0,2r), with r of
squared norm four. We generate 127,800 distinct Γ72 minimum vectors
having product four with that normal. The search uses the three cases
c=r, c=r+u for norm-four u orthogonal to r, and c=2u for a norm-four u with
r·u=2. An integer LDL enumeration supplies the short Leech vectors and
one of its polarisation cosets. Every generated triple is rechecked with
integer arithmetic for membership, norm, distinctness and its product
with the first normal. This is the complete link even without relying
on exhaustive generation: the degree-eight design moments give exactly
127,800 vectors with this product. Counting the vectors in this complete
link whose product with the second normal is also four gives **x=1,187**.

## Counts and kissing geometry

For norm-eight normals with mutual product four, the three-shell harmonic
profile equations have 141 antipodally paired unknown bins and rank 140.
The exact arithmetic check regenerates the support restrictions, short
translations, low-order design moments and cross-shell harmonic theta
equations up to degree 38, verifies a rank-140 modular minor, and verifies
a one-parameter affine family satisfies every equality. The parameter x
is the (4,4) corner. The weighted theta premise is explained in
[Nebe's account of Venkov's theory](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf);
the weight for a harmonic polynomial of degree k is 36+k, and the
weighted series is divisible by Δ^4 because the minimum is eight.

The three norm-eight slices with products (−1,−1), (−1,0), (−2,−1) have
populations 610,980,300−36x, 610,980,300−36x and 139,219,696+48x.
Project them to the orthogonal complement of the two normals and
normalize. The Gram of the normals is [[8,4],[4,8]]. Because distinct
parent vectors have product at most 4, the six exact pair-type checks
show all distinct projected points have product at most 1/2 and none
coincide. Normalization thus gives a 70-dimensional kissing construction.
This no longer requires enumerating the hundreds of millions of vectors
in each chosen fibre.

## Replay

Run `python verifications/improved/dim70-gamma72-mixed-a2-normals/verify.py`
with assertions enabled, a C++17 compiler, NumPy, SymPy, and Node.js.
The replay regenerates the link in a temporary directory and uses the
sibling proof package's SHA-pinned Leech structure. There is no randomized
step, external download, or floating-point geometric test in this replay.
