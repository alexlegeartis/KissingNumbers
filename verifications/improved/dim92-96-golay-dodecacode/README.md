# A Golay–dodecacode augmentation in dimensions 92–96

## Exact counts

| Dimension | Published repository incumbent | New layer | Combined configuration |
|---:|---:|---:|---:|
| 92 | 29,797,911,144 | 324,009,984 | 30,121,921,128 |
| 93 | 42,951,379,608 | 983,564,288 | 43,934,943,896 |
| 94 | 73,016,281,804 | 983,564,288 | 73,999,846,092 |
| 95 | 109,523,504,196 | 983,564,288 | 110,507,068,484 |
| 96 | 206,160,530,240 | 2,708,471,808 | 208,869,002,048 |

These are kissing-number lower bounds using the pinned incumbent constructions and their inherited proofs. The appended layer is checked independently. This is not a claim of optimality or an exhaustive claim of worldwide priority.

[92D incumbent](https://github.com/alexlegeartis/KissingNumbers/tree/35370272757ec8d6a68a3fea61ca4b6b6e5f43d7/verifications/improved/dim92-four-row-partition-heads) [93D and 96D incumbents](https://github.com/alexlegeartis/KissingNumbers/tree/35370272757ec8d6a68a3fea61ca4b6b6e5f43d7/verifications/improved/dim93-96-e7-rs-root-augmentation) [94D and 95D incumbents](https://github.com/alexlegeartis/KissingNumbers/tree/35370272757ec8d6a68a3fea61ca4b6b6e5f43d7/verifications/improved/dim94-95-e6-e7-head-exchange).

## Outer code and geometric rule

Use 24 ordered four-coordinate blocks, with the original row numbering of the Z8 parent. For shorter dimensions keep the new points at zero in the 4 or 2 replaced blocks and put zero in the replacement head coordinates. A pattern is a length-24 word over four symbols, with exactly 12 nonzero entries. Its three nonzero symbols assign the three perfect matchings of the four coordinates. In a chosen matching, the four projective axes are the two signed diagonals on each of its two coordinate pairs. Scale each active axis to give exactly two coordinates ±2. Each point has 24 nonzero coordinates and squared norm 96.

For two patterns let t be the number of doubly occupied blocks and s the number of those with the same nonzero tag. Their Hamming distance is 24−t−s. The block dot product is at most 8 for equal tags and 4 for unequal tags. The latter follows because different matchings share at most one coordinate between their chosen pairs. Thus total dot is at most 4(t+s). If the pattern distance is at least 12 this is at most 48 for every choice of signs and axes.

For each pattern, use the 4,096 projective-axis words of the additive quaternary dodecacode on its active rows, ordered by ascending block number. Its minimum Hamming distance 6 means at most six coincident axes for distinct words, so their dot product is at most 48. For a fixed word, independently use the [12,8,3] binary code given by parity-check columns 1 through 12 as four-bit vectors. There are 256 sign words and any two differ in at least three positions, giving dot product at most 8(12−6)=48. Each pattern thus supplies 2^20 distinct vectors.

## Compatibility with the inherited code

Each new active block has two entries of absolute value 2. The three incumbent dense row alphabets contain either four entries of absolute value 1 or a single entry of absolute value 2. The dot product with any such row is at most 4 on an occupied block. At most twelve blocks are occupied, so the total is at most 48. In 92D the dense code replaces the first four rows with its own head and retains the other twenty. In 93–95D the body retains the last twenty-two rows. The new points vanish on the omitted rows as well as on all head coordinates. Thus the head construction changes nothing in the dot bound.

The inherited Reed–Solomon points have six entries of amplitude 4. Their dot product with the new layer is at most 6·4·2=48. The inherited roots have two entries of amplitude sqrt(48), so the dot bound is 4 sqrt(48)<48. The new points cannot equal either sparse type since their supports have size 24. They cannot equal any dense old point since their body has empty rows and each old body has all its rows nonzero. Tags are recovered from the two-position support in each active block, and the axes and signs are recovered from its values, proving injectivity.

## Outer witness search and encoding

A pair (a,b) of binary extended Golay words defines a quaternary word with the two binary entries interpreted as a label 0 through 3, using bit a as the low bit. The binary code uses the 12 repository generator rows modulo 2. If a,b,a+b are all nonzero, their weights are at least eight and the union weight is one half their sum. Thus weight 12 occurs precisely for three octads. When one vanishes the other two are a single dodecad. There are 220,248 such vertices before shortening. Two distinct patterns conflict exactly when the difference in one of those three binary directions is zero and the other two are an octad. This means the search graph has only the three types of exact-coordinate fibres. The selected lists below were found by growing and exchanging an independent set, seeded by an affine quotient construction. The verifier independently checks all actual quaternary pair distances and does not assume any search optimization is correct.

There is also a compact algebraic 96D witness of 1,884 patterns. Let H be the binary message subspace spanned by [1798,1057,2678,460,2714], where the 12 bits are generator coefficients. Its complement basis is [1,2,4,8,16,32,64]. All 31 nonzero words of H have physical weight at least 12. Let f give the seven complement coefficients. Select (a,b) with f(b)=M f(a), where the columns of M are [52,70,112,26,103,34,91]. Both M and I+M are invertible. If a difference has one binary direction zero then the other belongs to H, and all other differences have union weight at least 12. The weight-12 restriction gives 1,884 patterns and 1,975,517,184 points. This smaller construction is an alternative to the larger selected list.

## Replay

Run `node verify.js` and `python3 verify_independent.py` from this directory or provide their paths from the repository root. The latter is a separately written standard-library replay. The three base-36 lists in `patterns.json` encode (a+4096b), where (a,b) are binary Golay message indices. Pair distances, support weights, both inner codes, and the block geometry are checked exactly. The mathematical compatibility proof above certifies the interactions with the pinned incumbent layers, without expanding billions of vectors.
