# A separate long-radius donor shell in dimension 93

**49,546,479,128 kissing directions**, an increase of **819,200** over 49,545,659,928 (about **0.00165%**). The addition keeps every preceding point. This is an explicit lower-bound improvement, not a claim of unrestricted optimality or of external publication priority. It uses 32 long heads, each carrying 25 exact cosets of size 1,024, instead of trying to find a 47th full short-head fibre.

The [preceding 46-fibre 93D design](https://github.com/alexlegeartis/KissingNumbers/tree/6e9b639b64e985af96b054e541481cb17c16fbaa/verifications/improved/dim93-z8-split-d5-heads) gives 49,392,123,904 dense words, 153,534,224 body-only directions and 1,800 head-incident D93 roots. The [published exact generator](https://github.com/alexlegeartis/KissingNumbers/blob/6e9b639b64e985af96b054e541481cb17c16fbaa/verifications/improved/dim93-z8-variable-radius-heads-continuation/z8_generator.json) is unchanged (Git blob b9b7faca47e39c04091a73098f7e82dd5c76b065; SHA256 903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73). We inherit homogeneous minimum 24 from [Duursma, Greferath, Litsyn and Schmidt (2001)](https://doi.org/10.1109/18.923742), not from the rank check. Their 96D binary nonlinear code, shortened by three binary coordinates, yields the named earlier 93D comparator 17,179,869,184 used by the project table.

## Why a separate long shell exists

Retain the same 22 physical four-coordinate blocks, body norm squared 88, and three cross-polytope frames. In one colour the row dot is `4 - 2w(a-b)`, where the homogeneous weight function is `w(0)=0`, `w(4)=4`, and 2 otherwise. For different colours, each row dot is ±2, so the body cap is 44. In the same colour, the parent minimum 24 bounds the body dot by 48 minus the actual two-row prefix dot. With one fixed prefix, two different body points have dot at most 40. This caps the squared head radius at 8 in the full-fibre proof, which is why looking only for complete fibres misses the new shell.

Put a new head `h=(14/5)(±1,±1,±1,±1,±1)`, with squared head norm 196/5 and total squared norm 636/5. A parent subset with conditioned homogeneous minimum 32 has body dot at most 24, so two points with the same new head have dot at most `24+39.2=63.2≤63.6`. The bigger radius raises the allowed mixed dot with the old short heads. The 32 new heads have exact compatible prefix labels; the checker tests every new/old mixed inequality by signed integer squaring.

There is an actual separation of radii in this proof method. If compatibility is proved with all forty old D5 heads using the conservative cross-fibre body envelope b≥44, and r denotes squared head norm, the two largest absolute coordinates have sum S satisfying

```
44 + 2S ≤ sqrt(24(88+r)),     r ≤ (5/4) S².
```

The second bound follows by maximizing `x₁²+4x₂²` over fixed sum `x₁+x₂=S` and `x₁≥x₂≥0`. Hence for r=5a² one needs `13a²−44a+22≥0`, with roots `(22±3 sqrt(22))/13`. Intermediate radii are forbidden to this interface proof, not necessarily to the exact pair geometry with stronger body estimates. The balanced long directions lie beyond the larger root, while the minimum-32 subgroup still allows r≤40.

## The coset code and phase square

Work with the pinned 12×24 triangular generator G over Z/8Z. For a fixed first-two-symbol prefix, the first two coefficients are fixed by its invertible pivots; the last ten range independently. Let U be rows 2 through 11, zero-based. The subgroup 4〈U〉 has 1,024 words whose body homogeneous minimum is 32. The checker generates its entire binary body code G₀, checks rank ten and minimum Hamming weight 8. The quotient by the subgroup within a fixed prefix is (Z/4Z)¹⁰ with little-endian coefficient digits.

Each literal quotient coset includes all 1,024 words. For a difference x in the body, let A be the positions where x is nonzero modulo 4 and B the positions where x=4 outside A. The minimum weight over the coset is exactly

```
2 |A| + 4 min_{q in G₀} |(q XOR B) minus A|.
```

Indeed adding four keeps weight 2 unchanged on A and toggles weight 0 versus 4 on its complement. The checker evaluates that formula without rounding.

The 25 quotient labels W and four offsets Z are literal in the replay. The 300 inner differences all meet weight 32. Group the 32 head signs by whether coordinates 2 and 3 are negative. For the group sign masks 0,12,4,8 on those coordinates use respective packed prefixes 8,16,24,32 in colour 1, with packed digit 8p₀+p₁. For each of the eight choices o of free signs in coordinates 0,1,4, translate W+Z by a common shift Sₒ. Each group has eight disjoint copies, each of 25 cosets.

If d head signs differ, the head squared distance is 4(14/5)²d and the body squared distance is four times the homogeneous distance. The required minimum body distances for d=0,1,2,3,4,5 are 32,24,17,9,1,0:
- One shared head: the 25 quotient cosets and their internal subgroup clear 32.
- Different heads on the same prefix: the translations prevent duplicate bodies; the full inherited minimum 24 then applies to the body.
- Different groups and d=1: the free sign labels coincide, so Sₒ cancels. The replay tests all 25² coset pairs on each of the four sign-square edges, verifying distance 24.
- Different groups and d≥2: their second prefix symbols agree so at most 4 homogeneous weight was deleted. The inherited body minimum 20 exceeds the required 17.

The 32 distinct heads and their 800 disjoint cosets supply exactly 819,200 points, not 32 times 8¹⁰.

## Remaining interfaces and exact replay

All predecessor body-only interface bounds depend only on the unchanged body frame alphabets and the floor 88 on dense squared norm. The new norm exceeds 88; the replay evaluates the exact published inequalities for the Golay, matching, Reed–Solomon, seven-row flat, four-row triality and D4 layers. For the head-incident D93 roots, two head entries have joint magnitude at most 28/5; one head plus one body entry has joint magnitude at most 14/5+2. Both are below the normalized cap. Neither a sparse point, a root, nor an old dense direction can duplicate a new point. The old dense Gram caps are explicitly rechecked.

Run the standalone `verify.mjs` with Node; in hosted V8 replace its last `console.log` with `text`. It reconstructs the literal pinned source and witness geometry, the 1,024-word subgroup, the quotient distance gates, the prefix algebra, and all new/old geometry. Use `python verify.py` for the repository replay, which first checks the retained predecessor and then uses Node.js for the new shell. The new checker requires no solver and prints `{"dimension":93,"old":49545659928,"retainedDense":49392123904,"retainedBodyOnly":153534224,"retainedRoots":1800,"heads":32,"cosets":800,"cosetSize":1024,"added":819200,"lowerBound":49546479128,"edgePairs":2500}`.

The transferable observation is that a larger near-corner head may require a much stronger *conditional* body code but can be compatible with the entire old short-head sector. Here the throughput is controlled by the outer Z4 quotient code and four prefix offsets rather than the long-stalled 47th complete short-head fibre.