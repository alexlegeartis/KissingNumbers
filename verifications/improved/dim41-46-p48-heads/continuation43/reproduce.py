"""Rebuild the exact 43D point certificate from the retained predecessor and delta.

No optimizer, numerical search or network is needed. This checks coordinate
assembly, not angles: run both mathematical verifiers as well.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def assembled():
    with np.load(HERE/'previous_points.npz', allow_pickle=False) as source:
        before = dict(source)
    with np.load(HERE/'delta.npz', allow_pickle=False) as source:
        delta = dict(source)
    assert all(v.dtype.kind in 'iu' for v in before.values())
    assert all(v.dtype.kind in 'iu' for v in delta.values())
    denominator = int(delta['head_denominator'])
    previous_denominator = int(before['head_denominator'])
    assert denominator > 0 and previous_denominator > 0
    assert denominator % previous_denominator == 0
    deleted = delta['deleted_indices']
    assert np.array_equal(deleted, np.unique(deleted))
    assert np.all((deleted >= 0) & (deleted < len(before['body_indices'])))
    kept = np.setdiff1d(np.arange(len(before['body_indices'])), deleted)
    heads = np.r_[before['head_numerators'][kept].astype(object) *
                  (denominator//previous_denominator), delta['added_head_numerators']]
    assert all(abs(int(v)) < 2**63 for v in heads.flat)
    return dict(body_indices=np.r_[before['body_indices'][kept], delta['added_body_indices']],
                head_numerators=heads.astype(np.int64),
                head_denominator=np.array(denominator, dtype=np.int64))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    expected = assembled()
    path = HERE.parent/'partial43/explicit_points_43.npz'
    with np.load(path, allow_pickle=False) as current:
        assert set(current.files) == set(expected)
        assert all(np.array_equal(current[k], v) for k,v in expected.items())
    if args.output:
        if args.output.resolve() == path.resolve():
            parser.error('choose an output outside the registered certificate')
        np.savez_compressed(args.output, **expected)
    with np.load(HERE/'previous_points.npz') as before:
        previous_count = len(before['body_indices'])
    with np.load(HERE/'delta.npz') as delta:
        deleted_count = len(delta['deleted_indices'])
        added_count = len(delta['added_body_indices'])
    print(json.dumps(dict(status='PASS', previous_tiny_points=previous_count,
                         deleted_points=deleted_count, inserted_points=added_count,
                         explicit_tiny_points=len(expected['body_indices']),
                         gain=len(expected['body_indices'])-previous_count,
                         certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest()),
                     indent=2))


if __name__ == '__main__':
    main()
