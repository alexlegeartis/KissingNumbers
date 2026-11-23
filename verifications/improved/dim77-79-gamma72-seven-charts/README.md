# Seven-frame Γ72 tail fibres in dimensions 77–79

| Dimension | Exact lower bound | Increase over preceding construction | Tilted fibres | Pure tail points |
|---:|---:|---:|---:|---:|
| 77 | 6,222,276,760 | 263,436 | 30 | 40 |
| 78 | 6,223,917,240 | 797,572 | 42 | 72 |
| 79 | 6,226,104,558 | 1,322,542 | 58 | 126 |

These are constructive lower bounds, not optimality claims. The comparison is the [rephased half-root data at e7daacf8](https://github.com/alexlegeartis/KissingNumbers/blob/e7daacf8ce180b8f01f722dda82af7edfcb1b454/verifications/improved/dim77-84-gamma72-half-root-rephased/phase_claims.json), pinned locally as `predecessor.json`. We keep the original Γ72 norm-eight shell and replace its previous nonzero-tail layers.

## Seven phase charts from an automorphism

In the exact three-Leech decomposition of Γ72, let (A_F,B_F=1-A_F) be its two integral similarity operators, with (A_FB_F=2). Reflection in the subspace ({(x,x,xA_F):x∈Λ}) sends ((u,v,z)) to ((u-p,v-p,z-pA_F)), where (p=(u+v+zB_F)/2). The Γ72 parity conditions make (p) integral. The verifier checks that the reflection is integral in the 72D lattice basis and preserves the published Gram. Its action together with block permutations and sign changes produces seven charts for the integral-block 36D half-phase subspace. The new charts contain phases outside the starting subspace, and their depths are transported by the exact automorphisms.

Start with one explicitly admissible pair of binary A labels (a,b). Use the phases of zero B syndrome for every ordered pair of elements of the two-dimensional space ({0,a,b,a+b}), in all seven charts. Their union has exactly 64 phase classes: zero, 21 of horizontal squared depth four and 42 of depth six. The depth-six tests are the original Leech norm-four and norm-six exclusion; the pinned theta replay gives exactly 136,704 shortest vectors in each such isotropic half-coset.

## Tail root selection

Embed the doubled D5, E6 or E7 roots in doubled E8 coordinates. A witness row selects one antipodal root direction, divided by two so its squared length is two, and assigns it a single depth-six phase for both signs. Put the entire 136,704-point shortest shell of that phase over each sign. Distinct rays have original tail inner products of absolute value at most one. For a nonorthogonal pair, the verifier proves their phase difference has squared depth six; for orthogonal pairs, depth at least four suffices. The pair separation is therefore at least squared eight including the vertical difference. The same-ray antipodes are already vertically separated by squared eight. Finally the full doubled D5, E6 or E7 root set gives the indicated pure vertical points, which are automatically compatible with all tilted fibres by Cauchy–Schwarz.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` from this package. To replay elsewhere, append the path to the pinned `dim73-82-gamma72-depth6-lattice` source package. Requires NumPy and SymPy.
