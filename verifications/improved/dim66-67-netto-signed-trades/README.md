# Signed Netto quadruple trades in 66D and 67D

| Dimension | Prior | Inserted | Removed | New lower bound | Gain |
|---:|---:|---:|---:|---:|---:|
| 66 | 342,389,558 | 162 | 146 | **342,389,574** | +16 |
| 67 | 355,195,391 | 77 | 64 | **355,195,404** | +13 |

The exterior BCH caps and nonary layers are completely unchanged. This trades
some of the surviving signed Netto four-support points for new signed
four-support points of the same norm. For such a new vector, the signed
support consists of four coordinates of absolute value 4, giving squared
norm 64. Its dot products with the old weight-sixteen layer are at most
4·4·2=32, with the old dense equatorial layer at most 4·4=16, with the
axes at most 4·8=32, and with the existing BCH caps at most
max(16s,12s+4r)<32 for s=sqrt(37)/8 and r=3sqrt(3).

A distinct signed four-support point violates the dot bound only if it has
a different support and agrees on a signed triple. A nine-support point
violates it only if it contains the entire four-support with all signs
matching. Each chosen new point is checked against both rules; distinct
signings at the same new support have dot at most 32.

Run `OPENBLAS_NUM_THREADS=1 python verify.py`. The verifier expands the
old Netto support orbit and reconstructs all actual pentacode and partial
nine-support signings from pinned sibling packages. It verifies the complete
prior lists of retired Netto signings, then projects the nine-support phases
against the new support list, computes every new point's exact blockers,
and checks all simultaneous substitutions. Each line of a `win*.txt`
gives four ascending zero-based coordinates, a four-bit negative mask in
that order, and its old blockers encoded by `16 × support index + mask`.
Net gain is the number of new points minus the size of the union of old
blockers. Search used columns with at most three original quadruple blockers;
the displayed bounds do not depend on the search algorithm's optimality.
