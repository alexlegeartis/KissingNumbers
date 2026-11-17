# Compatible support-eight layer alongside RS9 support-ten

Keep the independently audited RS9 support-ten construction with72 nonlinear signs and all D_n roots. Its coordinate groups are indexed0,...,9, first nine groups having9 coordinates and final group having k=min(9,n−81). Add the following NEW support-eight points, each of norm²45.

Select four group indices whose bitwise XOR is0. There are18 such four-subsets of0,...,9; any two intersect in at most2 groups, because three entries determine the fourth. Exactly4 contain group9.

On the nine local coordinates of each group use the nine near-perfect matchings of K9. Round c consists of pairs{c+r,c−r} mod9, r=1,...,4. Give its four pairs arbitrary fixed labels0,...,3. In each selected group choose one pair from the SAME round c, with the four pair labels having XOR0. There are64 such label tuples and any two differ in at least2 positions. Retain only supports within the available k coordinates of the final group.

Two supports with different four-group sets intersect in at most4 coordinates. With the same group set but different rounds, pairs in each group intersect in at most1 coordinate (each edge of K9 belongs to exactly one round), hence at most4 overall. Within the same group set and round, different label tuples agree in at most2 groups, hence intersect in at most4. Thus all distinct supports intersect in at most4.

On each support take the128 even-parity sign vectors with entries ±sqrt(45/8). Same-support signs differ in at least2 places, giving inner product≤22.5; distinct-support overlap≤4 gives the same bound. Against any RS9 support-ten point, overlap is at most4, since the new support uses only4 groups and RS9 uses one coordinate per group. The cross-inner-product is therefore≤4*45/sqrt80<22.5. Against D_n roots it is≤2*sqrt(45/2)*sqrt(45/8)=22.5.

Each new support has at most2 head coordinates (all head coordinates lie in the final group). Dense cross-inner-products are therefore at most(6+2sqrt(3/2))*sqrt(45/8)<22.5. Other coordinates have magnitude≤1. Support sizes2,8,10 and at least42 distinguish roots,new layer,RS9 layer,and dense code, so there are no duplicates between families.

For each of14 group-quadruples not involving the final group, there are9*64=576 supports. For each of4 involving it, the sum of available pairs over all nine rounds is C(k,2), and every allowed pair label has16 extensions to the XOR-zero four-label code. Total supports are14*576+64*C(k,2). The script independently enumerates and deduplicates them. Multiply by128, and ADD this layer to the previously audited nonlinear72 RS9 total.

This new additive layer is pending independent audit. Its purpose is further refinement, not a claimed large-factor advance.

## Stronger final version: Steiner quadruple system

This version REPLACES the18 XOR-zero group quadruples above by30 blocks of a Steiner quadruple system S(3,4,10). They are all PGL2(F9) images of{infinity,0,1,2} on P1(F9). The verifier generates the30 blocks and checks directly that every triple occurs exactly once, so any two blocks meet in at most2 groups. This finite verification, not a cited design theorem, is sufficient.

There are18 blocks not containing group9 and12 containing it. Thus the FINAL support count is18*576+192*C(k,2), replacing the earlier14*576+64*C(k,2). All matching, XOR pair-label, even-sign, and cross-family angle proofs are unchanged. Use certificate.json totals, not the superseded18-block totals. This is one replacement layer; do not union both versions.
