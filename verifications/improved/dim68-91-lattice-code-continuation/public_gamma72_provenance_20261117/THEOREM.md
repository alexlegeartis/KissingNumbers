# Gamma72: source identification and complete norm-8 fiber weights

## Source and conventions

The source is Gabriele Nebe and N. J. A. Sloane's [Catalogue entry for Gamma72](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html). The entry supplies a72-dimensional Gram matrix G and states minimum8 and kissing number6,218,175,600. The downloaded HTML's SHA256 is92daeaf0d93672e82ca5e5bc4479e5775daf59029aa4e983fb464e36ef002044. All5184 matrix entries were independently parsed and found identical to the supplied catalogue_gram72.json. source.json contains the matrix, a precise canonical serialization/hash, and the exact matrices below. This identifies the actual lattice, not merely its dimension or theta series.

Vectors are ROW vectors. Let F be the24-dimensional Gram matrix in source.json and write x²=xFxᵀ. The integer matrices A,B satisfy
A+B=I, AB=BA=2I, AFAᵀ=BFBᵀ=2F.
For
T = [[I,I,B],[0,A,A],[0,0,2I]],
the verifier checks the EXACT identity
2G = T diag(F,F,F) Tᵀ.
Thus Gamma72 consists precisely of triples (a,b,c)/sqrt(2) in the three orthogonal F-spaces satisfying
**aB≡cB (mod2), and b≡c+aA (mod2).**
Indeed a=u,b=u+vA,c=uB+vA+2w gives these congruences. Conversely set u=a and v=(b−a)B/2; the congruences make v and w=(c−aB−vA)/2 integral, proving sufficiency.

The source minimum8 is an explicit external premise. Since every (0,0,2x)/sqrt(2) occurs, it implies the F-lattice has minimum at least4. F is positive definite, even, and determinant1 (checked exactly). Its rank24 theta series is therefore the standard rootless even-unimodular series E4³−720Delta. Its populations at squared norms4,6,8 are respectively196560,16773120,398034000. This standard theta-series theorem is the other explicit external premise; the verifier reproduces these coefficients arithmetically.

## Leech parity classes: multiplicities proved, not assumed

Let L be the F-lattice and V=L/2L. The quadratic form is q(x)=x²/2 mod2. There are2^24 parity classes. If distinct non-antipodal x,y have the same parity, then both (x+y)/2 and (x−y)/2 are nonzero lattice vectors and have norm at least4.

For norm4 or norm6 this is impossible: their two norms sum respectively4 or6, less than8. Consequently such classes contain at most the two vectors ±x. A norm4 class and a norm8 vector cannot coincide, by the same argument with sum6. Norm6 classes have odd q, while norms4 and8 have even q, so these types are disjoint. The zero class has no nonzero vector below norm16. For two non-antipodal norm8 vectors of equal parity, the two midpoint norms are4±(x,y)/2, forcing (x,y)=0. Hence a norm8 parity class contains at most24 antipodal orthogonal pairs, or48 vectors.

The independently known populations force equality everywhere, because
1 +196560/2 +16773120/2 +398034000/48 =2^24.
Thus every parity class is exactly one of: zero;98280 type2 classes with two norm4 vectors;8386560 type3 classes with two norm6 vectors;8292375 type4 classes with48 norm8 vectors. These assertions concern the indicated lowest shells, not absence of larger vectors in a class.

## The two finite parity facts

Put U=ker(B mod2). The verifier checks U=im(A mod2), dimension12. It is totally isotropic since (xA)²/2=x² is even. Nondegeneracy of q follows from det(F)=1, so U is maximal isotropic. Every nonzero U-class is type4: a norm4 vector with aB even would equal (aB/2)A, contradicting the minimum4 of L. Thus U has4095 nonzero type4 classes.

The bundled196560 DISTINCT integer norm4 witnesses are checked directly against F. Their number equals the theta coefficient, so they comprise the complete norm4 shell; no correctness of floating enumeration pruning is required. Their parity keys give98280 classes, each twice. Applying B to these keys gives all4095 nonzero B-values, with EXACTLY24 type2 classes over every such value. This is the finite24-per-fiber fact required below.

For any nonzero B-value, its affine fiber is c+U with c not in U. On this fiber q(c+u)=q(c)+(c,u) mod2. Since U=U-perpendicular and c is not in U, this linear functional is nonzero. Therefore exactly2048 of the4096 classes have each q-value. The verifier additionally checks nontriviality of the functional for a representative of EVERY4095 nonzero B-fiber. The identity AF+FAᵀ=F also gives, whenever aB=cB,
q(c+aA)=q(c)+q(a).

## Completion theorem

Fix an integer c in L. The number of ordered pairs (a,b) for which (a,b,c)/sqrt(2) is in Gamma72 and has squared norm8 is as follows:

| c² | Exact completion count |
|---:|---:|
|0|9,828,000|
|4|12,800|
|6|192|
|8|f8(c mod2)|
|16|1 if c is in2L; otherwise0|
|all other values|0|

Here f8 is an explicitly defined finite parity table. Let K be the set of98280 type2 parity keys, let K_g={t in K:tB=g}, and let XOR denote addition modulo2. Then
f8(s)=4 sum over nonzero g of #{(t,u) in K_g²:s=t XOR u XOR g} +96*1[s in U and s≠0].
The verifier rebuilds all2^24 entries and checks the SHA256 against the frozen fiber8_table.bin. Thus its normalization and indexing are independently identified.

Proof by components: norm8 means a²+b²+c²=16, and nonzero component norms are at least4 and even.

* If c=0, the non-pure case has a²=b²=8 and equal nonzero U-parities. It contributes4095*48². The two pure possibilities (a,0,0) and (0,b,0), with a or b twice a norm4 vector, contribute2*196560. Total9,828,000.
* If c²=4, the possible component norms are (4,8),(6,6),(8,4). Each of the24 type2 classes in its B-fiber has two a-vectors. The resulting even b-parity must be type4: a zero or type2 b-parity would produce a Gamma72 vector of norm below8, contradicting the SOURCE minimum. Hence the two ordered (4,8)/(8,4) cases contribute48*96=4608. The2048 odd-q a-classes pair with odd-q b-classes, each with2*2 representatives, giving8192. Total12,800. Interchanging a,b preserves the congruences, as A+B=I and AB=0 modulo2.
* If c²=6, only (4,6) and (6,4) occur. There are48 ordered type2/type3 parity pairs, each with2*2 representatives. Total192.
* If c²=8, (4,4) gives the ordered-pair term4N22 above. The cases (0,8) and (8,0) occur exactly for a nonzero U-parity, each contributing48. This gives f8.
* If c²=16, a=b=0; the congruences force c≡0. Conversely every such c works. If c² is10 or12, the remaining positive component would have norm6 or4 and B-value0, impossible because a vector in LA has norm twice an even norm of L and minimum8. Norm14 leaves insufficient norm; larger norms are impossible. This exhausts every case.

As a check independent of slice enumeration, f8 is supported on type4 parities, and its sum over parity keys is9,828,000. Therefore the total norm8 population obtained from these weights is
9,828,000 +196560*12800 +16773120*192 +48*9,828,000 +196560
= **6,218,175,600**, agreeing with the published source total.

## Verification and scope

Run `python verify.py` (NumPy required). Optionally supply a separately downloaded source HTML path to recheck extraction of G. All matrix identities, witness norms/distinctness, parity fibers, quadratic nontriviality, completion-table hash and total arithmetic are exact integer checks. The norm4 witness count closes completeness using the independently derived theta coefficient, rather than trusting a floating enumerator. This theorem certifies the norm8 completion weights used in weighted slice counting; it does not by itself certify any floating slice traversal, any norm10 completion formula, or any additional kissing construction.
