"""Two complete exact replays of the refilled P48 class family.

The first uses int64 Gram products and canonical line hashes. The second
uses binary64 products with an explicit exactness bound and compares the
full canonical coordinates, without hashes. No search state is trusted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest
from published import TAU
from october_family import load_refilled


def canon(rows):
    assert np.all(np.any(rows != 0, axis=1))
    first = np.argmax(rows != 0, axis=1)
    return rows * np.sign(rows[np.arange(len(rows)), first])[:, None]


def product(a, b):
    assert a.dtype == b.dtype == np.int64
    bound = a.shape[1] * int(np.abs(a).max()) * int(np.abs(b).max())
    assert bound < 2**63
    return a @ b


def exact_float_product(a, b):
    # Bound every product and every partial sum, including cancelled terms.
    bound = a.shape[1] * int(np.abs(a).max()) * int(np.abs(b).max())
    assert bound < 2**53
    return a.astype(np.float64) @ b.astype(np.float64)


def verify():
    path = ROOT / 'october/family_refilled.npz'
    data = load_refilled()
    assert all(value.dtype.kind in 'iu' for value in data.values())
    gram = data['gram'].astype(np.int64)
    assert np.array_equal(gram, np.load(ROOT / 'data/p48p_gram.npy'))
    base = data['base'].astype(np.int64)
    assert base.ndim == 2 and base.shape[1] == 48 and len(base) > 7678
    assert len(np.unique(canon(base), axis=0)) == len(base)
    weighted = product(base, gram)
    assert np.all(np.sum(weighted * base, axis=1) == 6)
    histogram = np.zeros(7, dtype=np.int64)
    for start in range(0, len(base), 384):
        dots = np.abs(product(weighted[start:start+384], base.T))
        for i in range(len(dots)):
            assert np.all(dots[i, :start+i] <= 2)
            assert np.all(dots[i, start+i+1:] <= 2)
            histogram += np.bincount(dots[i, start+i+1:], minlength=7)
    assert int(histogram.sum()) == len(base)*(len(base)-1)//2

    sizes, offsets = data['sizes'], data['offsets']
    transforms, retained = data['transforms'], data['retained']
    assert len(sizes) == len(transforms) == 577
    assert transforms.shape[1:] == (48, 48) and offsets.shape == (578,)
    assert retained.ndim == 1 and offsets[0] == 0 and offsets[-1] == len(retained)
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
        # Modulo 2^64 arithmetic is intentional here. Distinct hash values
        # imply distinct lines; any collision rejects the certificate.
        hashes.append(rows @ data['projection'])
    hashes = np.concatenate(hashes)
    assert len(np.unique(hashes)) == len(hashes) == int(sizes.sum())

    # Independent base-pair calculation and full-coordinate disjointness.
    # Every binary64 product below is an exact integer by the bound above.
    with threadpool_limits(limits=4):
        other_weighted = exact_float_product(base, gram)
        for start in range(0, len(base), 512):
            dots = np.abs(exact_float_product(other_weighted[start:start+512], base.T))
            dots[np.arange(len(dots)), start+np.arange(len(dots))] = 0
            assert np.all(dots <= 2)
        coordinates = []
        for i, transform in enumerate(transforms):
            assert np.array_equal(exact_float_product(
                exact_float_product(transform, gram), transform.T), gram)
            indices = retained[offsets[i]:offsets[i+1]]
            rows = canon(exact_float_product(base[indices], transform))
            assert np.all(rows == np.rint(rows)) and np.abs(rows).max() < 128
            coordinates.append(rows.astype(np.int8))
        coordinates = np.concatenate(coordinates)
        # Compare all 48 bytes, not a hash or a search cache.
        keys = np.ascontiguousarray(coordinates).view(np.dtype((np.void, 48))).ravel()
        assert len(np.unique(keys)) == len(keys) == len(hashes)

    old = np.load(ROOT / 'october/family_refilled_7678.npz', allow_pickle=False)['sizes']
    old, new = np.sort(old)[::-1], np.sort(sizes)[::-1]
    rows = []
    for k in range(1, 14):
        count = 378 if k == 12 else TAU[k]//2
        previous = 52416000 + 2*int(old[:count].sum()) + TAU[k]
        value = 52416000 + 2*int(new[:count].sum()) + TAU[k]
        assert value > previous
        rows.append(dict(dimension=48+k, classes=count,
                         selected_lines=int(new[:count].sum()),
                         previous_cap_bound=previous, cap_bound=value,
                         gain=value-previous))
    return dict(status='PASS: two complete exact geometric replays',
                base_certificate=dict(lines=len(base), norm_squared=6,
                                      max_abs_off_diagonal=2,
                                      pair_histogram=histogram.tolist()),
                number_of_classes=len(sizes), distinct_mother_lines=len(hashes),
                class_sizes=sizes.tolist(), rows=rows,
                base_pairs_checked_twice=int(histogram.sum()),
                disjointness='int64 canonical hashes and independent full 48-coordinate comparison',
                certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='49-60D registered bounds; 61D cap construction remains superseded by code layers')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true')
    args = parser.parse_args()
    if not __debug__:
        raise SystemExit('Run without -O: verification assertions must be enabled.')
    if not args.write_report:
        check_manifest(ROOT)
    report = verify()
    path = ROOT / 'october/family_refilled.verified.json'
    if args.write_report:
        path.write_text(json.dumps(report, indent=2)+'\n')
    else:
        assert report == json.loads(path.read_text())
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
