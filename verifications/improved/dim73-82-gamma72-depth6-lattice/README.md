# Depth-six Γ72 lattice kissing configurations in dimensions 73–82

The exact shortest shells of ten explicit lattices give the following lower bounds.
These improve the preceding repository bounds; they are not claims of optimality.

| Dimension | Explicit kissing points | Increase over previous repository bound |
|---:|---:|---:|
| 73 | 6,218,449,010 | 242,078 |
| 74 | 6,218,995,830 | 781,973 |
| 75 | 6,219,816,924 | 1,484,822 |
| 76 | 6,221,459,976 | 2,820,530 |
| 77 | 6,221,737,226 | 2,716,350 |
| 78 | 6,222,291,726 | 2,556,884 |
| 79 | 6,223,124,388 | 2,502,244 |
| 80 | 6,224,790,624 | 3,028,346 |
| 81 | 6,225,064,194 | 2,834,750 |
| 82 | 6,225,611,462 | 1,131,770 |

The preceding bounds are pinned at [commit 424a5de](https://github.com/alexlegeartis/KissingNumbers/tree/424a5de8a27da94f3c5789f30ddaaf3516494fc8). [Exact comparison data](comparisons.json) separately retain the named reference bounds.

## Construction and minimum

Start with [Nebe's even unimodular lattice Γ72 of squared minimum 8](https://arxiv.org/abs/1008.2862). [Nebe's catalogue Gram](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/neb72.html) is reproduced as a canonical integer array. The verifier reconstructs its three-Leech-block model and checks the exact Gram identity and its binary glue. The metric is half the direct sum of the three Leech metrics. One polarized half gives two differences ((a,b)) between the block labels and the other half gives their total (t). They specify the horizontal shift modulo Γ72. A literal integer Leech-block triple represents a half shift of Γ72 since doubling it meets all binary glue equations.

For vertical squared norm 2 we require horizontal depth at least 6, so no Leech triple of total norm below 12 can represent the specified shift. The independent nonzero labels (a,b) rule out one occupied block. The exact norm-four shell reports the possible second labels for each first label and excludes all 4+4 combinations. For a 4+6 combination, the norm-four label and total second label annihilate the two-dimensional first-label plane. The norm-six vector's two labels would therefore pair to zero, contradicting its nonzero quadratic value. For vertical squared norm 4 only a single Leech block can violate separation. The verifier excludes norm four by shell lookup and norm six by the mod-two quadratic form. For norm 6 a nonzero glue class suffices. The parent has its attained squared minimum 8, so the full lifts do as well.

The tails are (A_1,A_2,A_3,D_4,D_4+A_1,D_4+A_2,D_4+A_3,D_4+D_4,D_4+D_4+A_1,D_4+D_4+A_2), with minimum squared norm 2. The two (D_4) blocks use opposite permutations of a three-point first-label plane and can reuse the same kernel of second labels. The final (A_1,A_2) come from an independent first-label plane. This rules out the dangerous cross-block phase cancellations; all short vertical states and their phases are checked directly. The dimensions 73–78 specialize the [corresponding density certificates](https://github.com/alexlegeartis/SpherePackings/tree/a0ecc6aef01d50b4038e89e0ca2b18b9e5307eb6/research/gamma72_depth6_20261122).

## Exact shell counts

The parent shell has 6,218,175,600 points. For a norm-two vertical coefficient, exactly **136,704** horizontal contacts occur. In fact any isotropic Γ72 half coset of depth at least 6 has this number. To see this, let (u) be its doubled representative in Γ72 and let (f(z)=sum_{xinGamma_{72}}(-1)^{ucdot x}q^{|x|^2/2}). The character is trivial on (Gamma_0(2)) for isotropic (u). Poisson summation gives (f|W_2=2^{18}Theta_{Gamma_{72}+u/2}(2z)), using weight 36. The ten forms (A^{18-2j}B^j), (0le jle9), with (A=2E_2(2z)-E_2(z)), (B=E_4(z)-4E_4(2z)), span by the valence bound and the verified independence of their first ten coefficients. Both forms are Fricke-odd. The first four coefficients of (f) are (1,0,0,0), and the shifted series has six vanishing coefficients. The resulting invertible linear system forces the claimed number. Our shifts are isotropic because an integer block-coordinate triple has even total Leech norm.

In the other layers the horizontal norms are only 4, 2 or 0. Modulo twice the Leech lattice, each norm-four or norm-six class has precisely two representatives of that norm. A norm-eight class has at most 48 representatives: non-antipodal representatives must be orthogonal by applying the Leech minimum to their half-sums and half-differences. The exact Leech theta coefficients fill all remaining classes with 48 each. The code counts the single-block and two-block lifts, and checks the zero-horizontal cases by the glue equations.

Replay with `OPENBLAS_NUM_THREADS=1 python verify.py` (NumPy and SymPy). It reconstructs all projected norm-four labels from an integer frame, checks the horizontal and vertical minimum gates, solves the theta coefficient systems, and verifies [the exact counts by layer](results.json).

`python audit_cartesian.py` repeats the vertical enumeration in Cartesian root coordinates and checks the layer totals independently.
