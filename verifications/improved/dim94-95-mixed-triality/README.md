# Mixed triality continuation in dimensions 94 and 95

| Dimension | Matching-row predecessor | Appended | New lower bound |
|---:|---:|---:|---:|
| 94 | 74,139,180,748 | 1,937,408 | 74,141,118,156 |
| 95 | 117,088,854,084 | 1,937,408 | 117,090,791,492 |

Every predecessor point is preserved. Add 473 six-row patterns (26 with one axial row, 359 with two, 88 with three), each with 4,096 points. The 22 ordinary four-coordinate rows and zero head are identified with original row indexes 2 through 23. Each line of witness.txt gives the three disjoint masks for axial, even-parity-flat and odd-parity-flat rows. Each row frame has four projective axes of norm two and their negatives. Use 128 length-six quaternary axis words of Hamming distance at least three and all 32 even sign words. All new vectors have squared radius 24.

Pairs of new outer patterns have raw dot at most twice the sum of their common occupied rows and their common same-frame rows. The exact checker verifies this does not exceed 12. For the same outer pattern, at most three equal axes or at least two sign reversals give the same bound. Against the earlier ordinary, paired, and hexary BCH supports the checker uses the row-frame caps directly, and also uses the support-group bound for the earlier flat Reed–Solomon layer. The new vectors use both axial and flat rows, so are distinct from the earlier pure flat and matching rows. Against the preceding matching append, every shared row contributes at most 2√2; the checker compares this raw envelope to each actual product of radii. The older dense and shifted Golay/dodecacode vectors have squared radius 96 and per-shared-row dot at most four.

The physical amplitude-four Reed–Solomon coordinate groups have lengths 14,16,16,16,16,16. In dimension 95 the exceptional coordinate 88 is omitted and is zero in every new point. Importantly a four-coordinate row can straddle two of these groups; this proof checks actual intersection with each physical group, not a four-row approximation. Each group is met by at most one axial row and the groupwise maximum coordinate budget is at most six in aggregate, giving raw mixed dot ≤24 against amplitude-four Reed–Solomon. The root supports use either two coordinates in one group, or one coordinate and the excluded head axis, so they contribute at most 3√48<24. Each axial row has just one nonzero coordinate of magnitude two.

Run python3 verify.py with the standard library, assertions enabled. The proof uses the explicit witness and source profiles, not the candidate search. The 96D construction is not needed; only the small row-frame and inner-code formulas are reused.
