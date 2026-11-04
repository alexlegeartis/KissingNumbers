# Classical root-system operator

`lattice.root_system` takes no inputs and parameters `family` (`A`, `D`, or
`E`) and positive integer `rank`. It generates roots from formulas and replays
the shared exact coordinate verifier. It reads no saved record coordinates.

`A_n` has `n(n+1)` roots, obtained from all `e_i-e_j` in the zero-sum subspace
of `R^(n+1)`. Orthogonal integer Helmert coordinates give an explicit
`n`-coordinate representation with a positive integral diagonal metric.
For `b_k=(1,...,1,-k,0,...)`, let `L=lcm(k(k+1))`. The returned coordinates are
`x dot b_k`, weights are `L/(k(k+1))`, and squared norm is `2L`. The weighted
Gram matrix is exactly `L` times the original root Gram matrix.

`D_n`, for `n>=2`, consists of all four signs on each pair of coordinate axes:
`2n(n-1)` points of squared norm 2. These two families accept any rank allowed
by explicit execution budgets; there is no threshold at dimension 8 or 24.

`E_8` uses integer roots of squared norm 8: the 112 signed vectors `2e_i+2e_j`
and the 128 all-sign vectors having even sign parity. `E_7` is the section
`x_6+x_7=0`; deleting `x_7` and giving `x_6` weight 2 is an exact isometry.
`E_6` additionally imposes `x_5+x_7=0`; its last three coordinates are
`(a,a,-a)`, so retaining `a` with weight 3 gives exact six-dimensional
coordinates. Their cardinalities are 126 and 72. The exceptional family exists
at ranks 6, 7, and 8, which is a mathematical family definition.

The operator checks point count, dimension and all-pair work budgets before
generating any coordinate list. The shared verifier then independently checks
integer entries, positive weights, distinctness, every norm, every pair, and
exact span rank. `dimension` and `ambient_dimension` are both the requested
rank, including `A_n`, `E_6`, and `E_7`; no uncomputed rank claim is inherited.
Large A/D systems may exceed the default explicit point or pair budgets. Such
an execution refusal is a resource limit, not a mathematical rank restriction.

Recipes `roots01.json` through `roots08.json` replay A1, A2, A3, D4, D5,
E6, E7, and E8. They are classical construction controls with cardinalities
2, 6, 12, 24, 40, 72, 126, and 240. This operation does not certify global
optimality, literature priority, or coverage of every root or lattice family.
