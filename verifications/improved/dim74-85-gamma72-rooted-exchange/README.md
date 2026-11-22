# Rooted exchange for the Γ72 norm-ten latitude

The latitude vertex sets and backbone stay unchanged. The new proof couples
the three harmonic counting systems for a fixed first body. This lowers
the maximum allowed triangle count in the conflict graph and increases twelve
kissing-number lower bounds. The independent sets are proved to exist, not
materialized as coordinate files.

| Dimension | Preceding | Verified lower bound | Gain |
|---:|---:|---:|---:|
| 74 | 6,218,213,549 | **6,218,213,672** | +123 |
| 75 | 6,218,331,159 | **6,218,331,258** | +99 |
| 76 | 6,218,635,448 | **6,218,636,373** | +925 |
| 77 | 6,219,013,264 | **6,219,014,956** | +1,692 |
| 78 | 6,219,716,626 | **6,219,721,604** | +4,978 |
| 79 | 6,220,590,652 | **6,220,599,336** | +8,684 |
| 80 | 6,221,713,348 | **6,221,725,998** | +12,650 |
| 81 | 6,222,197,896 | **6,222,202,604** | +4,708 |
| 82 | 6,224,427,690 | **6,224,435,960** | +8,270 |
| 83 | 6,227,374,552 | **6,227,388,086** | +13,534 |
| 84 | 6,231,140,058 | **6,231,161,058** | +21,000 |
| 85 | 6,235,824,404 | **6,235,855,650** | +31,246 |

## Rooted counting

Fix an oriented norm-ten body v. For j=4,5,6, let N_j be the number of
oriented norm-ten bodies w with v·w=j, and let f_j(a,b) be the number of
oriented norm-ten u with v·u=a and w·u=b, averaged over those w. The
values of N_j are independent of v: 12,497,420,000, 250,483,456 and
1,056,125, respectively. We likewise average the accompanying norm-eight
profiles over the same w. All identities from the parent harmonic proof
are linear, with coefficients depending only on j, so every averaged
profile satisfies those identities. At degree fourteen, M_2=0 gives
vanishing separately on both shells; the replay also includes the
norm-eight version of this identity.

For each j,a in {4,5,6}, swapping the roles of w and u gives the exact
double count

    N_j f_j(a,b) = N_a f_a(j,b).

Any exceptional head, such as w or -w, is treated in its known bin, not an
unknown. This averaging is over neighbors of each fixed v separately.
Consequently every rooted triangle count is bounded, rather than only a
global average. No transitivity on the norm-ten bodies is needed.

The resulting system has 328 unknown nonexceptional antipodal bins.
The verifier recomputes a rank of 325 modulo 49,999 and checks a rational
affine space with three independent free coordinates satisfies every
identity. These two checks prove that the affine space contains every
actual rooted profile. For each of dimensions 77–85 the triangle argument adds
nonnegative weights on three nonnegative profile bins to cancel the three
free-coordinate coefficients. Hence the remaining constant is an exact
upper bound on the corresponding rooted cost.

The integer tail palette and angular comparisons are regenerated exactly.
For a fixed root tail the cost first sums over the other tail choices;
the first tail can be any palette vector because signed permutations
act transitively on this palette. One representative per body line is used
in this graph. The signed third-body sum is twice the projective sum by
tail sign symmetry, and each triangle contributes twice to the ordered
neighbor count. Together these give the factor four between the replay's
cost accumulator and its triangle cap; the same-body contributions are
included with a safe overcount. The inherited blocked Boßecker–Rautenbach
recurrence then gives the selected latitude population.

Dimensions 74–76 instead use the fully oriented graph and the parent
13, 60 and 246 tail-ray palettes. The independent tail costs are computed
exactly for every first ray; there is no division by two for a choice of
body orientation. Each triangle still contributes twice to the ordered
neighbor count. Nonnegative rational profile multipliers bound each
root-tail cost; any remaining affine coefficient is bounded by the
norm-ten shell population, since all three free coordinates are counts
of norm-ten profiles. The verifier takes the worst rooted triangle cap
over all first rays and uses each ray's exact degree.

## Exact replay

Run:

    python verifications/improved/dim74-85-gamma72-rooted-exchange/verify.py

The package requires NumPy and SymPy, and does not need a numerical
optimizer or network access. The following sibling files are pinned by
SHA-256 in the verifier: the harmonic system generator, shell and graph
recurrence routines, and preceding results ledger. The original exact
construction is [Nebe's Γ72 lattice](https://arxiv.org/abs/1008.2862);
the [parent latitude proof](../dim74-85-gamma72-harmonic-latitudes/) gives
the inherited geometry and shell-design sources. See also
[Boßecker–Rautenbach](https://doi.org/10.1016/j.disc.2010.05.026).
