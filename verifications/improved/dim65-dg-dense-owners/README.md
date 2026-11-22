# A larger dense owner code in dimension 65

**65D kissing lower bound: 335,929,116**, an improvement of **65,536** over 335,863,580. Start from the entire 335,795,354-point 64D equator including the [latest partial sign additions](../dim64-65-partial-nonary-signs/). Replace the original 65,536 dense owners in the two-latitude construction by 131,072 dense owners. Keep the 2,688 sparse line owners and the two pure poles.

## The larger dense owner set

In the inherited dense BCH radical, label the 64 positions by the successive powers of the primitive element of F64, followed by zero, in its actual code order. The verifier reconstructs the radical and checks that the constant, six linear, and fifteen square-free quadratic evaluation rows all lie in it. Thus the parent contains the whole binary RM(2,6) in this coordinate convention.

Use the Delsarte–Goethals construction; see [the binary background and alternating-form extension in Fernández-Alcober, Kwashira and Martínez-Moro](https://digibuo.uniovi.es/dspace/bitstream/10651/56772/1/Nonassociativeringsandnonlinearcodes.pdf). For each a,b in F32 take the symmetric binary form B(x,y)=Tr(a xy + b(x²y + xy²)) on F32. The alternating form on six bits has off-diagonal entries B_ij+B_ii B_jj on the first five coordinates and B_ii in the last column. Use its upper-triangle entries as the coefficients of x_i x_j in an ordinary binary quadratic word, then add any affine word. The replay checks 1,024 distinct quadratic parts and exhausts all 15,903 different nonzero symmetric-difference masks against every one of the 128 affine words. The minimum Hamming distance is exactly 24. This computation certifies the capacity and distance directly; the named construction is background, not an unverified coding-bound premise.

## Lift and exact count

Retain the known compatible 64D geometry as the equator, except that every selected dense or sparse owner u is removed and replaced by (√3 u/2,+4) and (√3 u/2,−4). Owners of the old dense set outside the new set go back into the equator; the old dense lifts are not counted. Distinct dense owners have raw dot at most 16 by the verified Hamming distance. Dense/sparse pairs also have raw dot ≤16 because the sparse owners have amplitude 4 on four coordinates. The 336 GF(4) affine lines meet pairwise in at most one coordinate; each carries eight even-parity signs, so distinct sparse owners also have raw dot ≤16.

All norms are squared 64. Same-tail lifted dots are at most (3/4)·16+16=28. Opposite-tail lifted dots are at most (3/4)·64−16=32, including the two lifts of one owner. Against any retained equatorial point the old compatible-body cap 32 is contracted; the two poles have absolute products at most 8·4=32. The layers are disjoint.

The count is 335,795,354 + 131,072 + 2,688 + 2 = 335,929,116. Run `python verifications/improved/dim65-dg-dense-owners/verify.py` with NumPy installed (used by the inherited source module). The verification pins its source files and uses only exact binary-code arithmetic. The previous 64D construction and sparse lift remain inherited.
