#!/usr/bin/env python3
"""Exact consistency replay, conditional on the inherited body populations.

No optimizer or python-flint is needed for this replay. Optional
--replay-populations runs the original older population certificate in a
temporary copy, requiring python-flint and g++. See README.md for its scope.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from final import ROOT, ledger, results


def exact_rank(rows):
    """Rational row reduction; the supplied heads have at most eight columns."""
    basis = {}
    for row in np.unique(rows, axis=0):
        vector = [Fraction(int(value)) for value in row]
        for pivot, prior in sorted(basis.items()):
            multiple = vector[pivot]
            if multiple:
                vector = [a - multiple * b for a, b in zip(vector, prior)]
        pivot = next((i for i, value in enumerate(vector) if value), None)
        if pivot is not None:
            divisor = vector[pivot]
            basis[pivot] = [value / divisor for value in vector]
    return len(basis)


def verify_manifest():
    manifest = ROOT / 'data/SHA256SUMS.txt'
    listed = set()
    for line in manifest.read_text(encoding='utf-8').splitlines():
        expected, name = line.split('  ', 1)
        path = (ROOT / name).resolve()
        assert ROOT in path.parents and path != manifest
        assert name not in listed
        listed.add(name)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, name
    actual = {path.relative_to(ROOT).as_posix() for path in ROOT.rglob('*')
              if path.is_file() and '__pycache__' not in path.parts and path != manifest}
    assert listed == actual, ('manifest coverage', sorted(listed ^ actual))
    return len(listed)


def verify_populations(classes, radii, weights, geometry, tiny):
    for class_id, radius, weight in zip(classes, radii, weights):
        class_id, radius, weight = int(class_id), int(radius), int(weight)
        assert 0 <= class_id < 6561 and weight > 0
        mu = int(geometry['mu'][class_id])
        if radius == 10:
            assert mu % 24 == 8
            expected = int(tiny[class_id])
        elif radius == 12:
            assert mu % 24 == 0 and class_id != 0
            expected = 495
        elif radius == 14:
            assert mu % 24 == 16
            expected = 9180
        elif radius == 16:
            assert mu % 24 == 8
            expected = 116235 - 40 * int(tiny[class_id])
        elif radius == 18:
            assert mu % 24 == 0
            expected = 1092000 if class_id == 0 else 1084590
        else:
            raise AssertionError(('unsupported radius', radius))
        assert weight == expected, (class_id, radius, weight, expected)


def body_bounds(classes, radii, geometry, start, end):
    differences = ((geometry['digits'][classes[start:end], None, :] -
                    geometry['digits'][classes][None, :, :]) % 3) @ geometry['powers']
    bounds = radii[start:end, None] + radii - 18 + geometry['mu'][differences] // 4
    same = ((classes[start:end, None] == classes) &
            (radii[start:end, None] == radii))
    bounds = np.where(same, 2 * radii[start:end, None], bounds)
    bounds[np.arange(end - start), np.arange(start, end)] = 2 * radii[start:end] - 18
    return bounds


def verify_current(record, geometry, tiny):
    with np.load(ROOT / record['candidate'], allow_pickle=False) as source:
        model = dict(source)
    classes, radii, weights, heads = (model[key] for key in
                                     ('classes', 'r3', 'weights', 'head_numerators'))
    denominator = int(model['head_denominator'])
    pure = model.get('pure_directions', model.get(
        'pure_head_numerators', np.empty((0, 8), dtype=np.int64)))
    count = len(classes)
    assert count == len(radii) == len(weights) and heads.shape == (count, 8)
    assert pure.ndim == 2 and pure.shape[1] == 8
    assert all(array.dtype.kind in 'iu' for array in (classes, radii, weights, heads, pure))
    assert denominator > 0
    assert len(np.unique(np.c_[classes, radii, heads], axis=0)) == count
    verify_populations(classes, radii, weights, geometry, tiny)

    # Check the integer range BEFORE using fixed-width NumPy products/squares.
    head_max = max(abs(int(heads.min())), abs(int(heads.max())))
    norm_bound = 18 * denominator ** 2 + 8 * head_max ** 2
    assert (2 * norm_bound) ** 2 < 2 ** 63
    norms = radii * denominator ** 2 + np.sum(heads * heads, axis=1)
    assert np.all(norms > 0)
    for start in range(0, count, 64):
        end = min(count, start + 64)
        bounds = body_bounds(classes, radii, geometry, start, end)
        products = bounds * denominator ** 2 + 2 * heads[start:end] @ heads.T
        assert max(abs(int(products.min())), abs(int(products.max()))) ** 2 < 2 ** 63
        assert not np.any((products > 0) &
                          (products * products > norms[start:end, None] * norms))

    rank = exact_rank(heads)
    assert rank == exact_rank(np.r_[heads, pure]) == record['dimension'] - 40
    report = dict(dimension=record['dimension'], points=record['total'], families=count,
                  pure_points=len(pure), head_rank=rank, family_angles='PASS',
                  inherited_population_formulas='PASS', fresh_body_enumeration=False)
    if len(pure):
        pure_max = max(abs(int(pure.min())), abs(int(pure.max())))
        pure_bound = 8 * pure_max ** 2
        assert max((2 * pure_bound) ** 2,
                   (16 * pure_max * head_max) ** 2,
                   pure_bound * norm_bound) < 2 ** 63
        lengths = np.sum(pure * pure, axis=1)
        assert np.all(lengths > 0)
        cross = pure @ heads.T
        doubled = 2 * cross
        assert not np.any((doubled > 0) &
                          (doubled * doubled > lengths[:, None] * norms))
        pair = 2 * pure @ pure.T
        np.fill_diagonal(pair, 0)
        assert not np.any((pair > 0) & (pair * pair > lengths[:, None] * lengths))
        maximum = max((Fraction(int(value) ** 2, int(lengths[i]) * int(norms[j]))
                       for i, row in enumerate(cross) for j, value in enumerate(row)
                       if value > 0), default=Fraction(0))
        report.update(pure_angles='PASS', maximum_positive_cross_squared=str(maximum))
        if record['dimension'] in (45, 46):
            dim45 = record['dimension'] == 45
            assert len(pure) == (16 if dim45 else 27)
            assert set(map(int, lengths)) == ({80} if dim45 else {48})
            gram = pure @ pure.T
            assert set(map(int, gram[~np.eye(len(pure), dtype=bool)])) == (
                {-48, 16} if dim45 else {-24, 12})
            assert maximum == (Fraction(8, 45) if dim45 else Fraction(1, 6))
        else:
            assert int(model['pure_head_denominator']) == denominator
            assert np.all(lengths == 18 * denominator ** 2)
    return report


def verify_older_43(geometry, tiny):
    """Retain the independent check of the historical recovered seed."""
    entry = next(row for row in ledger()['constructions'] if row['dimension'] == 43)
    assert entry['status'] == 'candidate_supplied'
    with np.load(ROOT / entry['older_seed'], allow_pickle=False) as source:
        seed = dict(source)
    classes = geometry['old_to_new'][seed['body_classes']]
    radii, weights = seed['body_norm_times3'], seed['weights']
    verify_populations(classes, radii, weights, geometry, tiny)
    bounds = body_bounds(classes, radii, geometry, 0, len(classes))
    assert np.array_equal(bounds, seed['body_dot_bound_times6'])
    # The recovered seed has denominator 10^10: int64 arithmetic would overflow.
    heads = seed['head_numerators'].astype(object)
    denominator = int(seed['head_denominator'])
    assert denominator > 0 and exact_rank(seed['head_numerators']) == 3
    norms = radii.astype(object) * denominator ** 2 + np.sum(heads * heads, axis=1)
    products = bounds.astype(object) * denominator ** 2 + 2 * heads @ heads.T
    assert not np.any((products > 0) & (products * products > norms[:, None] * norms))
    assert sum(map(int, weights)) == entry['older_seed_total']
    return dict(dimension=43, older_seed_points=entry['older_seed_total'],
                latest_working_total=entry['total'], latest_certificate='partial43/verify.py',
                older_seed_angles='PASS', fresh_body_enumeration=False)


def verify_portable(record):
    """Replay the supplied arbitrary-size integer certificates in isolation.

    The upstream 43D script writes replay.json. Run in a temporary copy so
    verification preserves every installed source and pinned report.
    """
    source = ROOT / record['portable_replay']
    if record['dimension'] == 43:
        with np.load(ROOT / record['explicit_points'], allow_pickle=False) as points:
            indices, heads = points['body_indices'], points['head_numerators']
            assert indices.dtype.kind in 'iu' and heads.dtype.kind in 'iu'
            with np.load(ROOT / record['candidate'], allow_pickle=False) as model:
                width = model['head_numerators'].shape[1]
            assert heads.shape == (len(indices), width)
            assert len(np.unique(np.c_[indices, heads], axis=0)) == len(indices)
    else:
        with np.load(ROOT / record['candidate'], allow_pickle=False) as model:
            medium = set(map(int, model['classes'][model['r3'] == 12]))
            assert {77, 202, 239} <= medium
    with tempfile.TemporaryDirectory(prefix='kissing_p48_portable_') as directory:
        target = Path(directory) / source.name
        shutil.copytree(source, target, ignore=shutil.ignore_patterns('__pycache__'))
        replay = subprocess.run([sys.executable, '-B', str(target / 'verify.py')],
                                cwd=target, check=True, capture_output=True, text=True)
        report = json.loads(replay.stdout)
    assert report['status'] == 'PASS'
    assert report['dimension'] == record['dimension'] and report['points'] == record['total']
    report['fresh_body_enumeration'] = False
    return report


def replay_populations(jobs):
    """Preserve the pinned reports: the upstream replay writes its outputs."""
    with tempfile.TemporaryDirectory(prefix='kissing41_46_population_') as directory:
        target = Path(directory) / 'population_replay'
        shutil.copytree(ROOT / 'population_replay', target)
        subprocess.run([sys.executable, '-B', str(target / 'run_verify.py'),
                        '--jobs', str(jobs)], cwd=target, check=True)
    print('Older 41D/42D/44D population replay: PASS. Latest-family counts still use '
          'the inherited tables and POPULATIONS.txt; see README.md.')


def main():
    if not __debug__:
        raise SystemExit('Run without -O: verification assertions must be enabled.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay-populations', action='store_true')
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    assert args.jobs > 0
    print('SHA-256 manifest: %d files PASS' % verify_manifest(), flush=True)
    with np.load(ROOT / 'data/e8_geometry.npz', allow_pickle=False) as source:
        geometry = dict(source)
    tiny = np.load(ROOT / 'data/population_table.npz', allow_pickle=False)['tiny']
    assert geometry['digits'].shape == (6561, 8) and tiny.shape == (6561,)
    assert sum(map(int, tiny)) == 30240
    for record in results():
        report = (verify_portable(record) if 'portable_replay' in record
                  else verify_current(record, geometry, tiny))
        print(json.dumps(report), flush=True)
    print(json.dumps(verify_older_43(geometry, tiny)), flush=True)
    if args.replay_populations:
        replay_populations(args.jobs)
    print('PASS: current 41D-46D constructions, including partial 43D and coordinated '
          '45D; complete-family populations inherit the documented identities. '
          'No new 47D/48D claim.')


if __name__ == '__main__':
    main()
