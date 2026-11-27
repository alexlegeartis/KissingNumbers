# Independent exact replay of Philip Carder's 370,568-point 33D code

Audited 28 November 2026 Moscow time (27 November UTC). The author of the 370,568 point construction is **Philip Carder**, not Kravatskiy. It extends Kravatskiy's earlier 370,528 point stock (+40), and exceeds the newer independently registered Kravatskiy 370,556 witness by 12.

- Title: *A 370,568-point kissing configuration in dimension 33*, Zenodo, 27 November 2026, version 1.0.
- Cohn cites concept DOI https://doi.org/10.5281/zenodo.24577265; actual independently read immutable version DOI https://doi.org/10.5281/zenodo.24577266, created 2026-11-27T05:12:42.717164Z.
- Original https://zenodo.org/api/records/24577266/files/kissing-33-370568-submission.zip/content ; 22,928-byte ZIP from metadata. MD5 56e3e0ac2a39196a4522ffc76c052c82, matching metadata; SHA256 134bca4f8b814fd6f4d77b31ae1c3d19be4c870e02dcf8c48da18bb227cced16. Internal witness.npz SHA256 fa20e5fc12c7212a34f7fe574615915888c5c080cd2d554bc9598e5069ae199d. Third-party verify.py was not executed.
- Exact result: 131,072 distinct binary codewords using 32 distinct valid coordinates, minimum nonzero linear weight 8; 1,837 distinct 8-coordinate masks with pairwise overlap at most 4 and one valid sign-parity bit apiece; 4,360 individually listed integer vectors of squared norm 32 and largest sparse-sparse dot exactly 16 over all 9,502,620 unordered distinct pairs. The sparse absolute sum over the 32 dense axes has max 16. Across every sparse row and each octad support the absolute sum has max 8.
- Proof from these exact integers: dense vectors have squared norm 32 and dot ≤32−2·8=16. Eight ±2 entries give the same norm; within a support parity gives at least two sign flips and dot≤16, and distinct supports contribute at most 4·4. Dense-to-octad dot≤8·2. Sparse-to-dense dot≤16 by the absolute sum, and sparse-to-octad dot≤2·8. Sparse-to-sparse checked directly. All the layer unions are disjoint since duplicates would give dot 32, excluded by those same cross checks. Total 131,072+1,837·128+4,360=370,568, with every distinct normalized dot ≤1/2. Mathematically verified independently, not merely a report in Cohn's table, and not an optimality theorem.
- Program used NumPy integer matrix products with int32 values, not float arithmetic. For all sparse products the exact bound 32 on the absolute partial sums follows already from the norm screen and Cauchy–Schwarz; no overflow. The other enumerations use Python unbounded integers or integer axis-sums. No external program was run.

## Independent replay code
Pass the downloaded Zenodo ZIP as argv[1]. This variant uses explicit stop checks rather than optimized-away assertions.

```python
import sys
from hashlib import sha256
from io import BytesIO
from itertools import combinations
from zipfile import ZipFile
import numpy as np

def need(condition, label):
    if not condition:
        raise ValueError(label)

blob = open(sys.argv[1], "rb").read()
need(sha256(blob).hexdigest() == "134bca4f8b814fd6f4d77b31ae1c3d19be4c870e02dcf8c48da18bb227cced16", "archive")
with ZipFile(BytesIO(blob)) as z:
    raw = z.read("witness.npz")
need(sha256(raw).hexdigest() == "fa20e5fc12c7212a34f7fe574615915888c5c080cd2d554bc9598e5069ae199d", "witness")
with np.load(BytesIO(raw), allow_pickle=False) as w:
    masks = list(map(int, w["supports"]))
    parity = list(map(int, w["parity"]))
    gens = list(map(int, w["dense_generator"]))
    pos = list(map(int, w["dense_positions"]))
    S = w["sparse"].astype(np.int32)
    shift = int(w["dense_shift"])
need(len(set(pos)) == len(pos) == 32 and all(0 <= i < 33 for i in pos), "axes")
need(all(0 <= i < 2**32 for i in gens + [shift]), "dense bits")
C = [0]
for g in gens:
    C += [x ^ g for x in C]
need(len(set(C)) == len(C) == 131072, "dense cardinality")
dmin = min(x.bit_count() for x in C if x)
need(dmin >= 8, "dense distance")
need(len(masks) == len(set(masks)) == len(parity) == 1837, "distinct octads")
need(all(0 <= m < 2**33 and m.bit_count() == 8 for m in masks) and set(parity) <= {0, 1}, "octad masks/parities")
overlap = max((x & y).bit_count() for x, y in combinations(masks, 2))
need(overlap <= 4, "octad overlaps")
need(S.shape == (4360, 33) and np.all((S*S).sum(axis=1) == 32), "sparse norms")
dense_l1 = int(abs(S[:, pos]).sum(axis=1).max())
need(dense_l1 <= 16, "dense/sparse")
cols = np.array([[k for k in range(33) if (m >> k) & 1] for m in masks], dtype=np.int64)
octad_l1 = max(int(abs(S[i:i+160, cols]).sum(axis=-1).max()) for i in range(0, len(S), 160))
need(octad_l1 <= 8, "octad/sparse")
wmax = -99
pairs = 0
for i in range(0, len(S), 128):
    W = S[i:i+128] @ S.T
    for j, row in enumerate(W):
        if i + j + 1 < len(S):
            wmax = max(wmax, int(row[i+j+1:].max()))
            pairs += len(S) - (i + j + 1)
need(wmax <= 16 and pairs == 9502620, "sparse/sparse")
print("PASS", len(C) + len(masks)*128 + len(S), dmin, overlap, dense_l1, octad_l1, pairs, wmax)
```
