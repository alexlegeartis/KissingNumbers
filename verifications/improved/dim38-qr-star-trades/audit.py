"""Independent star/root enumeration; bound all other new pair types."""
from itertools import combinations, product
import json
from pathlib import Path
import numpy as np

P = Path(__file__).resolve().parent
rows = [tuple(map(int, s.split())) for s in (P/'stars.txt').read_text().splitlines()]
A = np.zeros((len(rows), 38), dtype=np.int16)
for i, (c, w, s, z) in enumerate(rows):
    assert w.bit_count() == 4 and not (w >> c & 1) and z & w == z
    assert 0 <= c < 38 and w < 1 << 38 and s in (-1, 1)
    for j in range(38):
        A[i,j] = (2*s if j == c else (-1 if z >> j & 1 else 1) if w >> j & 1 else 0)
assert len(np.unique(A, axis=0)) == len(A)
assert np.all(np.sum(A*A, axis=1) == 8)
G = A @ A.T
assert int(G[np.triu_indices(len(A),1)].max()) <= 4
R = np.zeros((2812, 38), dtype=np.int16)
for i, ((c, j), (s, t)) in enumerate(product(combinations(range(38), 2), product((-1, 1), repeat=2))):
    R[i,c], R[i,j] = s, t
B = A @ R.T
remove = np.any(B > 2, axis=0)
S = np.array([list(map(int,s)) for s in (P.parent/'dim37-38-qr-merge'/'supports_3088.txt').read_text().splitlines()], dtype=np.int16)
assert int((S @ np.abs(A).T).max()) <= 4
for h in (0,3,-3):
    weights = np.full(38,2,dtype=np.int16)
    weights[37] = abs(h)
    peak = int((np.abs(A) @ weights).max())
    assert 4*peak*peak <= 8*(148+h*h)
n = 652028+len(A)-int(remove.sum())
print(json.dumps(dict(status='PASS', dimension=38, points=n, new_stars=len(A),
    removed_roots=int(remove.sum()), tested_star_root_products=B.size,
    tested_star_octad_supports=S.shape[0]*len(A)),sort_keys=True))
