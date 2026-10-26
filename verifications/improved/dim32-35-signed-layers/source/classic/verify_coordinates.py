"""Independent coordinate-based checker (Python + NumPy).

Expands every middle-layer point and all roots, and multiplies integer
coordinate arrays for each exceptional pair type.  Dense/dense distances
use linearity; sparse pairs on <=4 shared coordinates have a direct bound.
Does not import the search code or verify.py.
"""
import itertools
import json
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent


def verify(path):
    started = time.monotonic()
    c = json.loads(path.read_text())
    n = c['dimension']
    groups = {}
    for mask in c['full_even_supports']:
        positions = [j for j in range(n) if mask >> j & 1]
        assert len(positions) == 8 and mask not in groups
        vectors = []
        for signs in itertools.product((-2, 2), repeat=8):
            if sum(x < 0 for x in signs) % 2:
                continue
            v = [0] * n
            for j, s in zip(positions, signs):
                v[j] = s
            vectors.append(v)
        groups[mask] = vectors
    for mask, minus in c['exceptional_antipodal_lines']:
        assert mask >= 0 and mask < 1 << n and minus & mask == minus
        v = [(-2 if minus >> j & 1 else 2) if mask >> j & 1 else 0 for j in range(n)]
        groups.setdefault(mask, []).extend((v, [-x for x in v]))
    group_arrays = {}
    for mask, vectors in groups.items():
        x = np.asarray(vectors, dtype=np.int32)
        assert len(np.unique(x, axis=0)) == len(x)
        assert np.all(np.sum(x*x, axis=1) == 32)
        g = x @ x.T
        np.fill_diagonal(g, -32)
        assert g.max() <= 16
        group_arrays[mask] = x
    checked = 0
    for a, b in itertools.combinations(group_arrays, 2):
        if bin(a & b).count('1') <= 4:
            continue
        dots = group_arrays[a] @ group_arrays[b].T
        assert dots.max() <= 16
        checked += dots.size
    middle = np.concatenate(list(group_arrays.values()))
    roots = []
    extension = c.get('heavy_pair_extension')
    heavy = tuple(extension['coordinates']) if extension else None
    for i, j in itertools.combinations(range(n), 2):
        for si, sj in itertools.product((-4, 4), repeat=2):
            if heavy == (i, j) and si == -sj:
                continue
            v = [0] * n
            v[i], v[j] = si, sj
            roots.append(v)
    roots = np.asarray(roots, dtype=np.int32)
    assert np.all(np.sum(roots*roots, axis=1) == 32)
    rg = roots @ roots.T
    np.fill_diagonal(rg, -32)
    assert rg.max() <= 16
    # Any root/middle dot is bounded by 2*4*2=16, and a dense/middle dot
    # by 8*1*2=16, from the verified coordinate shapes.
    assert np.all(np.count_nonzero(middle, axis=1) == 8)
    assert np.max(np.abs(middle)) == 2
    generators = np.asarray(c['dense_generator'], dtype=np.uint32)
    words = np.array([0], dtype=np.uint32)
    for generator in generators:
        words = np.concatenate((words, words ^ generator))
    assert len(np.unique(words)) == len(words) == 131072
    lut = np.asarray([bin(i).count('1') for i in range(256)], dtype=np.uint8)
    weights = lut[words.view(np.uint8).reshape(-1, 4)].sum(axis=1)
    assert weights[1:].min() >= 8
    cube = 1-2*((words[:, None] ^ np.uint32(c.get('dense_coset_shift', 0))) >> np.arange(32, dtype=np.uint32) & 1).astype(np.int32)
    heads = np.empty((0, n), dtype=np.int32)
    mixed_checks = None
    if heavy:
        vectors = []
        for body, tail in itertools.product((-1, 1), repeat=2):
            v = [body] * n
            v[heavy[0]], v[heavy[1]] = 3*tail, -3*tail
            vectors.append(v)
        heads = np.asarray(vectors, dtype=np.int32)
        hn = np.sum(heads*heads, axis=1)
        assert np.all(hn == n+16)
        maxima = {}
        for label, points in [('middle', middle), ('roots', roots), ('dense', cube)]:
            q = heads[:, :32] if label == 'dense' else heads
            dots = q @ points.T
            assert np.all(4*np.maximum(dots, 0).astype(np.int64)**2 <= 32*hn[:, None])
            maxima[label] = int(dots.max())
        hg = heads @ heads.T
        np.fill_diagonal(hg, -(n+16))
        assert np.all(2*hg <= n+16)
        mixed_checks = dict(maximum_integer_dots=maxima, head_head=int(hg.max()))
    total = len(words)+len(middle)+len(roots)+len(heads)
    assert total == c['points']
    return dict(dimension=n, verified_points=total, expanded_middle_points=len(middle),
                expanded_roots=len(roots), expanded_dense_points=len(words), expanded_heads=len(heads),
                exceptional_cross_dot_products=checked, mixed_norm_checks=mixed_checks,
                seconds=time.monotonic()-started, all_checks_pass=True)


if __name__ == '__main__':
    results = []
    for path in sorted(ROOT.glob('certificate*.json')):
        result = verify(path)
        results.append(result)
        print(json.dumps(result), flush=True)
    recorded = [{k: v for k, v in r.items() if k != 'seconds'} for r in results]
    (ROOT/'coordinate_verification.json').write_text(json.dumps(recorded, indent=2)+'\n')
