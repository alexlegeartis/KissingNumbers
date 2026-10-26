"""Independent NumPy check expanding all vectors in the new 35D witness.

No imports from the constructor or the symbolic verifier. Integer matrix
products check every new-head interaction with every base-layer point.
"""
import itertools
import json
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent


def verify(path):
    start = time.monotonic()
    c = json.loads(path.read_text())
    n = c['dimension']
    signed_supports = c['full_parity_supports']
    assert len(signed_supports) == len({mask for mask, _ in signed_supports})
    sign_words = list(itertools.product((-2, 2), repeat=8))
    parities = [sum(v < 0 for v in word) % 2 for word in sign_words]
    patterns = [np.asarray([word for word, p in zip(sign_words, parities) if p == parity],
                           dtype=np.int32) for parity in (0, 1)]
    middle = np.zeros((128*len(signed_supports), n), dtype=np.int32)
    for k, (mask, parity) in enumerate(signed_supports):
        assert 0 <= mask < 1 << n and parity in (0, 1)
        coordinates = [j for j in range(n) if mask >> j & 1]
        assert len(coordinates) == 8
        middle[k*128:(k+1)*128, coordinates] = patterns[parity]
        gram = patterns[parity] @ patterns[parity].T
        np.fill_diagonal(gram, -32)
        assert gram.max() <= 16
    for (mask, _), (other, _) in itertools.combinations(signed_supports, 2):
        assert bin(mask & other).count('1') <= 4
    assert np.all(np.sum(middle*middle, axis=1) == 32)
    # Distinct supports share at most four nonzero entries of absolute value 2.
    assert np.max(np.abs(middle)) == 2
    removed = {tuple(x) for x in c['removed_root_lines']}
    assert len(removed) == len(c['removed_root_lines'])
    roots = []
    for a, b in itertools.combinations(range(n), 2):
        for sa, sb in itertools.product((-4, 4), repeat=2):
            if (a, b, sa*sb//16) in removed:
                continue
            v = [0]*n
            v[a], v[b] = sa, sb
            roots.append(v)
    roots = np.asarray(roots, dtype=np.int32)
    assert len(roots) == 2*n*(n-1)-2*len(removed)
    assert np.all(np.sum(roots*roots, axis=1) == 32)
    gram = roots @ roots.T
    np.fill_diagonal(gram, -32)
    assert gram.max() <= 16
    words = np.array([0], dtype=np.uint32)
    for generator in c['dense_generator']:
        assert 0 <= generator < 1 << 32
        words = np.concatenate((words, words ^ np.uint32(generator)))
    assert len(words) == len(np.unique(words)) == 131072
    popcounts = np.asarray([bin(x).count('1') for x in range(256)], dtype=np.uint8)
    weights = popcounts[words.view(np.uint8).reshape(-1, 4)].sum(axis=1)
    assert weights[1:].min() >= 8
    coordinates = c['dense_coordinates']
    assert len(coordinates) == len(set(coordinates)) == 32
    assert all(0 <= j < n for j in coordinates)
    dense = np.zeros((len(words), n), dtype=np.int32)
    dense[:, coordinates] = 1-2*((words[:, None] ^ np.uint32(c['dense_coset_shift']))
                                 >> np.arange(32, dtype=np.uint32) & 1).astype(np.int32)
    assert np.all(np.sum(dense*dense, axis=1) == 32)
    # Coordinate shapes prove root/middle <=16, dense/middle <=16,
    # and dense/root <=8. Dense/dense follows from the independently
    # enumerated minimum nonzero codeword weight.
    representatives = np.asarray(c['head_antipodal_lines'], dtype=np.int32)
    heads = np.concatenate((representatives, -representatives))
    norms = np.sum(heads.astype(np.int64)**2, axis=1)
    assert np.all(norms > 0)
    maxima = {}
    products = 0
    for label, layer in [('dense', dense), ('middle', middle), ('roots', roots)]:
        dots = heads @ layer.T
        assert np.all(4*np.maximum(dots, 0).astype(np.int64)**2 <= 32*norms[:, None])
        maxima[label] = int(dots.max())
        products += dots.size
    gram = heads @ heads.T
    np.fill_diagonal(gram, 0)
    assert np.all(4*np.maximum(gram, 0).astype(np.int64)**2 <= norms[:, None]*norms[None, :])
    maxima['heads'] = int(gram.max())
    count = len(dense)+len(middle)+len(roots)+len(heads)
    assert count == c['points']
    return dict(dimension=n, points=count, verified=True,
                explicit_head_base_dot_products=products,
                maximum_dots=maxima, seconds=round(time.monotonic()-start, 3))


if __name__ == '__main__':
    results = [verify(path) for path in sorted(ROOT.glob('certificate*.json'))]
    (ROOT/'coordinate_verification.json').write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))
