"""Check all new interactions with the separately certified affine-shell code."""
import argparse
import json
from pathlib import Path
import numpy as np


def verify(root):
    assert __debug__, "Do not run the verifier with Python -O"
    Q = np.array(json.loads((root / 'gram.json').read_text()), dtype=np.int64)
    R = np.loadtxt(root / 'minimum_shell_refill.txt', dtype=np.int64, ndmin=2)
    assert Q.shape == (40, 40) and np.array_equal(Q, Q.T)
    assert R.shape == (8, 40) and np.max(np.abs(R)) <= 128
    assert np.max(np.abs(Q)) <= 128
    # These bounds make the matrix products and all accumulation int64-safe.
    L = Q @ R.T
    assert ((R @ Q * R).sum(1) == 6).all()
    self_dot = int((R @ L)[np.triu_indices(len(R), 1)].max())
    assert self_dot <= 3
    a = np.load(root / 'core.npz', allow_pickle=False)
    b = np.load(root / 'auxiliary.npz', allow_pickle=False)
    totals = 0
    mx = []
    for rows, N, n in [(a['coefficients'], 54, 1084590), (b['norm36'], 36, 40590), (b['norm30'], 30, 1134)]:
        assert rows.dtype == np.int8 and rows.shape == (n, 40)
        d = max(int((block.astype(np.int64) @ L).max()) for block in np.array_split(rows, max(1, len(rows) // 32768)))
        assert d <= 0 or 4*d*d <= 6*N
        mx.append(d)
        totals += len(rows)
    return dict(status='PASS', dimension=40, base=totals, new=len(R), total=totals+len(R),
                max_new_new=self_dot, max_core_54=mx[0], max_aux_36=mx[1], max_aux_30=mx[2],
                exact_cross_pairs=totals*len(R), exact_new_pairs=len(R)*(len(R)-1)//2)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('source', nargs='?', default=Path(__file__).resolve().parent, type=Path)
    args=parser.parse_args()
    print(json.dumps(verify(args.source), indent=2))
