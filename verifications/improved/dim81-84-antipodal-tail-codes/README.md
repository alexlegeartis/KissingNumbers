# 81D and 84D: antipodal support codes as Γ72 tails

| Dimension | Verified lower bound | Increase over preceding bound |
|---:|---:|---:|
| 81 | **6,251,523,810** | 2,744,114 |
| 84 | **6,291,169,688** | 9,327,824 |

Use Nebe's Γ72 norm-eight core with a new exact orthogonal tail.
The short code in the tail uses the coordinate axes and signed
weight-four supports. The support lists have pairwise intersection
at most two, so with the same normalization these are antipodal
kissing configurations. Take 153 rays in 9 dimensions from the
[Leech–Sloane (P_{9a}) construction](https://doi.org/10.4153/CJM-1971-081-3)
and 420 in 12 dimensions from the [Best weight-four support system](https://ir.cwi.nl/pub/6831/6831D.pdf).
The integer lists themselves are supplied; their exact pairwise Gram
checks are part of the verification and do not assume the source's counts.

Divide each norm-four vector by sqrt(2) to give a tilted tail of
squared norm two. 81D uses 153 distinct rays for the tilted layer,
and 84D uses 384 of 420 rays for it. All 153 or 420 rays also contribute
both antipodes in a pure norm-eight tail layer. Use the existing
horizontal 4+4+4, 4+8 and 6+6 profile graph separation from the
[half-sign package](../dim81-85-gamma72-half-sign/). Its criteria only
need the exact zero/nonzero tail inner products and the bound one on
absolute tail dot; the verifier recomputes both directly from the
integer Gram matrices.

Run `OPENBLAS_NUM_THREADS=1 python3 verify.py`. Dependencies are
NumPy, SymPy and little-endian C++17. These replace, not adjoin,
the preceding tail stock. They are exact constructive lower bounds
using Nebe's core, not optimality claims.
