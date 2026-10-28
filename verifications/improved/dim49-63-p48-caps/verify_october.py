"""Check the historical October family and both exact refilled-family replays."""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest


def canon(rows):
    assert np.all(np.any(rows != 0, axis=1))
    first = np.argmax(rows != 0, axis=1)
    return rows * np.sign(rows[np.arange(len(rows)), first])[:, None]


def product(a, b):
    # Bound every partial sum before using signed 64-bit integer products.
    bound = a.shape[1] * int(np.abs(a).max()) * int(np.abs(b).max())
    assert bound < 2**63
    return a @ b


def main():
    check_manifest(ROOT)
    with np.load(ROOT / 'october/family_55_61.npz', allow_pickle=False) as z:
        data = {key: z[key] for key in z.files}
    assert all(v.dtype.kind in 'iu' for v in data.values())
    gram = data['gram'].astype(np.int64)
    assert np.array_equal(gram, np.load(ROOT / 'data/p48p_gram.npy'))
    base = data['base'].astype(np.int64)
    assert base.shape == (7387, 48)
    assert len(np.unique(canon(base), axis=0)) == len(base)
    weighted = product(base, gram)
    assert np.all(np.sum(weighted * base, axis=1) == 6)
    histogram = np.zeros(7, dtype=np.int64)
    for start in range(0, len(base), 384):
        dots = np.abs(product(weighted[start:start+384], base.T))
        for i in range(len(dots)):
            off = np.r_[dots[i, :start+i], dots[i, start+i+1:]]
            assert np.all(off <= 2)
            histogram += np.bincount(dots[i, start+i+1:], minlength=7)
    expected = json.loads((ROOT / 'october/family_55_61.verified.json').read_text())
    assert histogram.tolist() == expected['base_certificate']['pair_histogram']
    sizes, offsets = data['sizes'], data['offsets']
    transforms, retained = data['transforms'], data['retained']
    assert len(sizes) == len(transforms) == 577 and transforms.shape[1:] == (48, 48)
    assert offsets[0] == 0 and offsets[-1] == len(retained)
    assert np.array_equal(np.diff(offsets), sizes) and np.all(sizes > 0)
    assert data['projection'].shape == (48,)
    hashes = []
    for i, transform in enumerate(transforms):
        transform = transform.astype(np.int64)
        assert np.array_equal(product(product(transform, gram), transform.T), gram)
        indices = retained[offsets[i]:offsets[i+1]]
        assert np.all((indices >= 0) & (indices < len(base)))
        rows = canon(product(base[indices], transform))
        assert np.all(np.sum(product(rows, gram) * rows, axis=1) == 6)
        # Hash arithmetic is deliberately modulo 2^64. Distinct hashes imply
        # distinct canonical lines; a collision would cause rejection.
        hashes.append(rows @ data['projection'])
    hashes = np.concatenate(hashes)
    assert len(np.unique(hashes)) == len(hashes) == expected['distinct_mother_lines']
    assert sizes.tolist() == expected['class_sizes']
    print(json.dumps(dict(status='PASS', base_lines=len(base),
                          pairs=int(histogram.sum()), classes=len(sizes),
                          distinct_mother_lines=len(hashes))))
    from verify_refilled import main as verify_refilled
    verify_refilled()


if __name__ == '__main__':
    main()
