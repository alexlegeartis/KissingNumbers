"""Count the exhibited direction pairs and selected disjoint head classes."""
from pathlib import Path
import json
import numpy as np
ROOT = Path(__file__).resolve().parent

def results():
    with np.load(ROOT / 'tail840_exact.npz', allow_pickle=False) as tail:
        A, B = tail['A'], tail['B']
    assert A.shape == B.shape == (840, 12)
    keys = {tuple(np.r_[a, b]) for a, b in zip(A, B)}
    assert len(keys) == 840 and all(tuple(-np.r_[a, b]) in keys for a, b in zip(A, B))
    pairs = len(keys) // 2
    with np.load(ROOT / 'optimized420.npz', allow_pickle=False) as heads:
        bits = np.unpackbits(heads['retained_bits'], axis=1)
        assert len(bits) == pairs == 420 and not bits[:, 7783:].any()
        assert np.array_equal(bits[:, :7783].sum(axis=1), heads['sizes'])
        selected = int(heads['sizes'].sum())
    poles = sum(bool(line.strip()) for line in (ROOT / 'source/coordinates841.txt').read_text().splitlines())
    total = 52416000 + 2 * selected + poles
    assert selected == 3063342 and poles == 841 and total == 58543525
    return [{'dimension': 60, 'total': total}]

if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
