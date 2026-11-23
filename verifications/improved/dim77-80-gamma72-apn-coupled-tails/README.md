# Complete D5, E6, E7 and E8 tails over Γ72

| Dimension | Exact lower bound | Increase over preceding bound | Antipodal root directions |
|---:|---:|---:|---:|
| 77 | 6,223,234,328 | 35,584 | 20 |
| 78 | 6,226,571,472 | 757,272 | 36 |
| 79 | 6,232,148,054 | 5,513,576 | 63 |
| 80 | 6,243,307,240 | 14,875,528 | 120 |

The [APN depth-six construction](../dim81-85-gamma72-apn-multipartite/) gives 384 phases of depth six as three colours over the seven-dimensional B shift space. We assign one distinct phase to every antipodal root direction of the complete exceptional tail. The D5 stock consists of all signed two-coordinate roots. The E8 stock is represented in doubled Euclidean coordinates, which are twice the actual vertical coordinates: two nonzero coordinates each ±2, or eight signs ±1 of even parity. E7 and E6 are the sections with the last two or last three coordinates equal, respectively. The checker verifies all vertical dots and the cardinalities directly. The pure roots are multiplied by two in actual coordinates.

Different colours give horizontal phase difference of depth six. Equal colours, different phases, give a nonzero B-only difference, of depth at least four. Since the vertical absolute dot of distinct rays is zero or one, a conflict is possible only for equal phase colours and nonorthogonal directions. For every such pair the later ray clips against the entire earlier fibre; no clipping is needed against its own antipode. A B-only short difference can only be raw squared norm eight on one Leech block. The checker counts all surviving shortest contacts in each clipped ray, covering profiles 4+4+4, 4+8 and 6+6. For 4+4+4, membership in the shifted norm-four class is the exact conflict test; the two-block profiles use the dot-four test against the 48 signed representatives of the blocker class. Cyclically permuting the Leech blocks rotates among the three nonzero colours. Vertical dot bounds within and against the pure root shell follow from the verified tail root dots. The original 6,218,175,600 core points are all retained.

## Coupled restoration

The directed whole-fibre clip is sufficient but not always necessary. Two triple contacts 4+4+4 conflict only when exactly one B label is shifted by the class difference, leaving the other two B labels and all A labels unchanged. Because the shifted distinct norm-four local rays have dot zero, this label test covers all incompatible sign choices. A triple cannot conflict with a two-block point: the only short difference is in a single block, and the other two unchanged blocks already force the changed component to keep norm four. The checker reconstructs every conflicting pair of triples by sorted 48-bit label lookup, verifies that the established directed clip agrees in total, then restores any missing triple with no retained conflicting neighbour in reverse fibre order. It directly rechecks the final pair list. The recovered groups add 592, 10,224, 37,552 and 161,552 points, respectively.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` using NumPy and SymPy. The checked finite integer-sum bound below 2^53 makes the integer-valued BLAS products exact. These are constructive lower bounds, not optimality claims.

## Extra norm-four tails

A totally singular four-dimensional subspace of E8/2E8 gives norm-four vectors with even dot products. The shortest vectors in a single class are antipodal or orthogonal; they use one common nonzero B-only shift. Distinct classes have dot at most two, so use distinct shifts. The verifier enumerates the singular subspaces and restricts to D5, E6, E7 or E8, choosing 26, 30, 84 or 240 tails, respectively. Each B-only shift gives exactly 144 horizontal contacts. The doubled-coordinate checks certify that every new tail also has dot at most two with any old norm-two tail; the old mixed phases have nonzero colour and horizontal squared depth at least six against these shifts. No old point is retired.
