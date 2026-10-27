"""Deterministically regenerate the disjoint family from its exhibited base.

    python scripts/refill_family.py --output /tmp/family_refilled.npz

The base class is explicit proof data in the installed certificate. This
script reproduces the family selection, rather than rerunning the timed
local search that discovered the base. Verification is separate.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]


def build(base):
    with np.load(ROOT / 'october/family_55_61.npz', allow_pickle=False) as data:
        transforms = list(data['transforms'])
        projection = data['projection'].copy()
        gram = data['gram'].copy()
    generators = np.load(ROOT / 'data/p48p_gens.npy').astype(np.int64)
    seed = 20261028
    rng = np.random.default_rng(seed)
    seen = {transform.tobytes() for transform in transforms}
    while len(transforms) < 2577:
        transform = np.eye(48, dtype=np.int64)
        for _ in range(30):
            generator = generators[int(rng.integers(len(generators)))]
            assert 48*int(np.abs(transform).max())*int(np.abs(generator).max()) < 2**63
            transform = transform @ generator
        if transform.tobytes() not in seen:
            seen.add(transform.tobytes())
            transforms.append(transform)
    transforms = np.array(transforms)
    hashes = []
    with threadpool_limits(limits=2):
        for transform in transforms:
            # Each product and partial sum is an integer smaller than 2^24.
            assert 48*int(np.abs(base).max())*int(np.abs(transform).max()) < 2**24
            rows = (base.astype(np.float32) @ transform.astype(np.float32)).astype(np.int64)
            first = np.argmax(rows != 0, axis=1)
            rows *= np.sign(rows[np.arange(len(rows)), first])[:, None]
            hashed = rows @ projection  # Intentional arithmetic modulo 2^64.
            assert len(np.unique(hashed)) == len(base)
            hashes.append(hashed)
    hashes = np.array(hashes)
    unique, ids = np.unique(hashes, return_inverse=True)
    ids = ids.reshape(hashes.shape).astype(np.int32)
    used = np.zeros(len(unique), dtype=bool)
    counts = np.full(len(ids), len(base), dtype=np.int64)
    chosen, retained, sizes = [], [], []
    for _ in range(577):
        index = int(np.argmax(counts))
        keep = ~used[ids[index]]
        fresh = ids[index, keep]
        used[fresh] = True
        chosen.append(index)
        retained.append(np.flatnonzero(keep).astype(np.int32))
        sizes.append(len(fresh))
        counts = np.count_nonzero(~used[ids], axis=1)
        counts[chosen] = -1
    return dict(base=base, transforms=transforms[chosen],
                retained=np.concatenate(retained), offsets=np.r_[0, np.cumsum(sizes)],
                projection=projection, sizes=np.array(sizes), gram=gram)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    installed = ROOT / 'october/family_refilled.npz'
    if args.output.resolve() == installed:
        parser.error('Choose a separate output so the installed proof is preserved.')
    sys.path.insert(0, str(ROOT))
    from october_family import load_refilled
    expected = load_refilled()
    base = expected['base'].copy()
    rebuilt = build(base)
    assert rebuilt.keys() == expected.keys()
    assert all(np.array_equal(rebuilt[key], expected[key]) for key in rebuilt)
    np.savez_compressed(args.output, **rebuilt)
    print('PASS: every regenerated array equals the installed certificate')


if __name__ == '__main__':
    main()
