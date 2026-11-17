# Fresh independent audit, 17 November 2026

94D kissing construction: **73,016,281,804** points of squared norm96 with all distinct inner products at most48. This is not a sphere-packing density claim.

## Construction checked

Input candidate.json SHA256 de956543b714631172ccde72301b6312031a1d8f03cd2a9bfa19a881d319e4c3 matches current producer candidate. There are68 distinct norm8 heads and distinct (color,prefix) assignments. Within a color new head Gram is <= original two-symbol crosspolytope Gram; all cross-color head products<=4. The six-dimensional map is (sqrt3*h0,h2,h3,h5,h6,h7), since h1=h4=-h0. It preserves Gram and all resulting coordinates have absolute value<=2.

The unchanged rank12 Z8 parent has8^12 distinct words, homogeneous minimum24, and exactly8^10 words per fixed two-symbol prefix. Replacing the first8 coordinates by the6D head preserves squared norm96 and only decreases same-color inner products. Different colors have22 four-coordinate body blocks with products<=2 each, total<=44; adding head product<=4 gives48. Same-prefix distinct words remain distinct and retain original inequality. Thus dense count68*8^10=73,014,444,032.

Sparse augmentation is1,835,008 RS signed points plus2,764 roots. Six coordinate groups of sizes14,16,16,16,16,16 partition94 coordinates. Evaluate degree<=3 polynomials overGF16 at six distinct field elements, restrict first value to14 symbols;14*16^3 supports. Each receives32 even-parity signs at amplitude4. Distinct supports meet at most3, same-support signs differ in>=2, giving all products<=48. Dense/RS products<=6*4*2=48. Roots sqrt48*(±e_i±e_j) lie within one coordinate group and meet each RS support at most1; root/RS and dense/root products<=4sqrt48<48. Distinct roots have product<=48. Sparse support sizes6/2 differ from dense body support at least22, so layers are disjoint.

Fresh check.py independently reconstructs alphabets and checks all2278 head pairs, GF16 irreducibility/no-zero-divisors and all65535 nonzero cubic root counts on the six evaluation points. It does not import producer checker or earlier independent head checker.

## Inherited proof replay and review

Fresh complete parent/augmentation driver replay exited0; inherited_replay.log retained. Reviewed parent_proof code: ring generator is free rank12 via triangular odd pivots. Even-residue words reduce to Z4 Lee minimum12, established by all2^24 words; nonzero odd residue has weight>=8. Residue weight>=12 immediately gives homogeneous weight>=24. For each of759 weight8 residues, shortened dual conditions exclude all5489 errors of outside Lee weight<=3, yielding total>=16+8. All arithmetic ranges are safe for stated dtypes. This fresh replay reuses the published independent parent implementation, not a second newly authored parent-minimum algorithm.

## Comparison

The preceding published 62-head bound was 66,573,830,860. The gain is 6,442,450,944 or 9.677152209474%. The preceding published package is commit 27327c5079d1fdeedd42da650956b3ccef35050e.

The repository comparison34,359,738,368 is a derived two-coordinate shortening of Duursma–Greferath–Litsyn–Schmidt2001 binary parent, not a tabulated exact94D record. Searches did not locate a stronger external94D bound, but this is not an exhaustive worldwide priority audit; Cohn table does not cover94D. Report improvement over live repository, and do not equate limited-source absence with global SOTA proof.

