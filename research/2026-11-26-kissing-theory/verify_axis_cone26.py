#!/usr/bin/env python3
"""Usage: python verify_axis_cone26.py /path/to/Takhanov/version2/certificates"""
import sys
from pathlib import Path
from itertools import product
from collections import Counter
import hashlib
import json
import numpy as np

root = Path(sys.argv[1])
src = root / 'data'
sys.dont_write_bytecode = True
sys.path.insert(0, str(root))
from common import leech_shell

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(src / 'd26_blocks.json') == '5018425664146ee2ace418021537b8b13863fc57c7dc64376ab6f9e811aa5162'
assert sha(src / 'd26_certificate.npz') == '35fd8143cf0c7e6b0146719485e00c86952cbd5e404c7881c5f6b92fa25b21ec'
blocks = json.loads((src / 'd26_blocks.json').read_text())
I = blocks['kernel_coordinate_indices']
J = sorted(set(range(24)) - set(I))
F = np.array(blocks['frame_rows'], dtype=np.int64)
C = [list(x) for x in blocks['five_coordinate_supports']]
with np.load(src / 'd26_certificate.npz', allow_pickle=False) as z:
    leech = leech_shell()
    keep = np.ones(len(leech), dtype=bool)
    keep[z['deleted_leech_indices']] = False
    L = np.pad(leech[keep], ((0, 0), (0, 2)))
    A = z['added_integer'].copy()
L, A = [x[np.any(x[:, J] != 0, axis=1)] for x in (L, A)]
assert (len(L), len(A)) == (191904, 8960)
def proj(x):
    q = x[:, I].astype(np.int64) @ F.T
    assert np.all(q % 32 == 0)
    return np.unique(q // 32, axis=0)
PL, PA = proj(L), proj(A)
P = np.unique(np.concatenate([PL, PA]), axis=0)
assert np.all(PL[:, 8:] == 0)
assert (len(PA), len(P)) == (328, 2601)
D = {tuple(row) for row in P.tolist()}
for chart in C:
    for pool in [PA, P]:
        S = np.abs(pool[:, chart])
        assert int(np.sum(S,axis=1).max()) == 4
        assert int(S[:,-1].max()) == 2
    for signs in product((-1,1),repeat=4):
        for sig in (-1,1):
            n = [0]*10
            for i,x in zip(chart[:4],signs): n[i] = x
            n[17-chart[-1]] = 2*sig
            assert tuple(n) in D
T = blocks['four_coordinate_supports']
pairs = [sorted(set(v)-{8,9}) for v in T if 8 in v and 9 in v]
heads = [set(v) for v in T if 8 not in v and 9 not in v]
assert len(pairs)==4 and len(heads)==10
nativ = {tuple(sorted(t[:4])) for t in C}
surv = {tuple(sorted(x)) for x in product(*pairs)
        if all(len(set(x)&h) <= 2 for h in heads)}
assert len(nativ)==4 and surv == nativ
print('PASS: exterior rows 191904+8960; 328 added/2601 total normals; four maxima (4,2); 128 literal signed normals; four surviving dominant supports of 16.')
