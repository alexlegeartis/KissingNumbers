# 60D: one matched norm-12 layer with shared blocker retirement

This construction raises the exact 60D kissing-number lower bound from **58,663,722** to **58,664,259**. It adds **626** new rays arranged in pairs on **123** norm-12 lattice heads, and charges at most **89** old blockers. The certified net gain is at least **537**. It reuses the completed [global-owner and tetrahedral-tail certificate](../dim50-60-p48-global-owner-refill/) and adds one latitude, not a separate packing assumption.

## New rays

For every listed signed row `Z` of squared norm 10 and its assigned integral P48p head `w` of squared norm 12, include the two points

`(w/(4√2), Z/4)` and `(−w/(4√2), Z/4)`.

The squared head and tail lengths are 3/8 and 5/8. Two different rows are compatible exactly when `|w·w′| + 2 Z·Z′ ≤ 16`; the two sign mates of a row have product 1/4. All equatorial directions are safe by the P48p minimum: any squared-norm-six lattice head `x` satisfies `|w·x| ≤ 6`, giving an inner product at most `6/(8√3) < 1/2`. The old 840 poles have the form `D/2`, and the replay checks `|Z·D| ≤ 4`. Every new head is orthogonal to the 16-head cube frame; each new tail has ℓ¹ norm at most 8, so the unupgraded cube layer is safe. The replay separately checks the actual sparse tetrahedra and the upgraded cube subcosets, including the sign pattern of the latter. The coordinates, thresholds and cardinalities are exact.

## Label permutation and one-point cost

The old body has 420 distinct antipodal direction labels. Assign them to its 420 classes by the recorded permutation. This preserves all old inequalities: within a class the label is shared; distinct classes use different labels with absolute mutual product at most 1/2; all interactions of the old body with the equator, poles and tetrahedral layers use bounds uniform in the label. No old owner is reassigned to a different class.

For an old signed head `x` in the class given label `D`, its two lifts are `(x/(2√2), ±D/4)`. A new pair conflicts with at least one of those two lifts precisely when `|w·x| + |Z·D| > 8`. Take the union of obstructing old owners, charging each only once even when it blocks many new rays. Remove both lifts of any obstructing actual owner and restore the old unit equatorial ray `x/√6`. The ancestor equator had omitted it and its norm-six ownership is unique; it is also disjoint from the reserved equatorial heads of the parent layers. Thus the net cost per removed owner is one. The restored rays are compatible with the whole old construction and the new layer by the same old bound and the norm-12 argument above.

To obtain a conservative exact charge, the checker uses an overinclusive stock: every signed ancestor row, even rows retired in the earlier exchanges, plus all 62,937 rows inserted or reassigned by the parent migration and refills. Therefore the computed blocker count is an upper bound on the actual cost. It needs no claim about the optimality of the matching or of the lower bound.

## Replay

From the repository root with NumPy installed:

```sh
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim60-p48-n12-blockers/verify.py
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim60-p48-n12-blockers/audit_angles.py
```

The first checker uses bounded exact integer-valued binary64 products; the optional second checker independently uses int64 for owner alignment and rational-square identities for the sparse and upgraded angles. Both accept `--siblings`, `--body`, and `--refill` to locate the pinned sources in an offline exported layout. They write no files.
