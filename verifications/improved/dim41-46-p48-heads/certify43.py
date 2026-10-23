#!/usr/bin/env python3
"""Fresh population and exact angle certificate for the recovered 43D seed.

This certifies 2,878,155 points, not the unavailable 2,911,990-point extension.
Published inputs: P48n has minimum 6, 52,416,000 minimal vectors, and its
minimal shell is an 11-design. Requires NumPy, python-flint and g++ (C++17).
All generated files go into a temporary directory unless --output is given.
"""
import argparse
import concurrent.futures
import contextlib
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from final import ROOT
from verify import exact_rank, verify_manifest, verify_older_43


def certify(work, jobs):
    from flint import fmpz_mat

    replay = work / 'population_replay'
    shutil.copytree(ROOT / 'population_replay', replay)
    sys.path.insert(0, str(replay))
    import geometry
    import verify_masses

    g = geometry.geometry()
    with (work / 'mass_replay.log').open('w') as stream, contextlib.redirect_stdout(stream):
        verify_masses.main()
    print('Parent embedding and exact population identities: PASS', flush=True)

    compact = dict(np.load(ROOT / 'data/e8_geometry.npz', allow_pickle=False))
    ai = np.array((4 * fmpz_mat(g['p']['A'].tolist()).inv()).tolist(), dtype=np.int64)
    assert np.array_equal(g['p']['A'] @ ai, 4 * np.eye(8, dtype=np.int64))
    assert not np.any(g['u'] @ ai % 4)
    powers = 3 ** np.arange(8, dtype=np.int64)
    ids = ((g['u'] @ ai // 4) % 3) @ powers
    oldnew = np.empty(6561, dtype=np.int64)
    oldnew[g['uclasses']] = ids
    assert len(set(map(int, oldnew))) == 6561
    assert np.array_equal(oldnew[g['uclasses']], ids)
    assert np.array_equal(oldnew, compact['old_to_new'])
    assert np.array_equal(compact['powers'], powers)
    assert np.array_equal(compact['digits'], (np.arange(6561)[:, None] // powers) % 3)
    assert np.array_equal(compact['mu'][oldnew], g['mu12'])

    seed = dict(np.load(ROOT / 'data/old43_seed.npz', allow_pickle=False))
    classes, radii, weights = (seed[key] for key in
                              ('body_classes', 'body_norm_times3', 'weights'))
    assert classes.shape == radii.shape == weights.shape == (113,)
    assert all(array.dtype.kind in 'iu' for array in (classes, radii, weights))
    assert set(map(int, radii)) == {12, 14, 16, 18}
    assert exact_rank(seed['head_numerators']) == 3
    expected = {}
    for old, radius, weight in zip(classes, radii, weights):
        key = (int(old), int(radius))
        assert key not in expected or expected[key] == int(weight)
        expected[key] = int(weight)

    # The verified design identities give the first shells at normal norms
    # 4 and 6 uniformly: 9,180 and 495 points. No modular-form lookup is used.
    for (old, radius), weight in expected.items():
        if radius == 12:
            assert g['mu12'][old] == 24 and weight == 495
        elif radius == 14:
            assert g['mu12'][old] == 16 and weight == 9180

    # Recount ALL thirteen larger shells used here, including the norm-six
    # centre. Construct a rigorously enclosing integer ellipsoid for each.
    # For Q=q, E=P^3 Q-M D M^t and eta=max row sum |E|,
    # ||w||^2 <= target*tr(Q^-1) when wQw^t=target. Consequently
    # w M D M^t w^t <= target*(P^3+eta*tr(Q^-1)).
    gram = fmpz_mat(g['q'].tolist())
    inverse = gram.inv()
    trace = sum((Fraction(str(inverse[i, i])) for i in range(40)), Fraction(0))
    scale = 2 ** 24
    chol = np.linalg.cholesky(g['q'].astype(float))
    lower = chol / np.diag(chol)[None, :]
    matrix = np.rint(lower * scale).astype(np.int64)
    diagonal = np.rint(np.diag(chol) ** 2 * scale).astype(np.int64)
    assert np.all(diagonal > 0) and not np.any(np.triu(matrix, 1))
    assert np.all(np.diag(matrix) == scale)
    error = (scale ** 3 * gram - fmpz_mat(matrix.tolist()) *
             fmpz_mat(np.diag(diagonal).tolist()) * fmpz_mat(matrix.tolist()).transpose())
    eta = max(sum(abs(int(error[i, j])) for j in range(40)) for i in range(40))
    inputs = work / 'count_inputs'
    inputs.mkdir()
    shells = sorted(key for key in expected if key[1] >= 16)
    assert len(shells) == 13
    for old, radius in shells:
        bound = 3 * radius * (scale ** 3 + eta * trace)
        enclosure = (bound.numerator + bound.denominator - 1) // bound.denominator
        center = (g['labels'][old] + 1) % 3 - 1
        with (inputs / ('%d_%d.input' % (old, radius))).open('w') as stream:
            stream.write('40 %d 3 1\n' % scale)
            for array in (diagonal[None, :], matrix, g['q']):
                for row in array:
                    stream.write(' '.join(map(str, row)) + '\n')
            stream.write('%d %d %d\n' % (old, 3 * radius, enclosure))
            stream.write(' '.join(map(str, center)) + '\n')
    exe = work / ('exact_count.exe' if sys.platform == 'win32' else 'exact_count')
    subprocess.run(['g++', '-O3', '-std=c++17', str(replay / 'exact_count.cpp'),
                    '-o', str(exe)], check=True)

    def recount(key):
        old, radius = key
        process = subprocess.run([str(exe), str(inputs / ('%d_%d.input' % key))],
                                 capture_output=True, text=True, check=True)
        result = json.loads(process.stdout)
        assert result['class_id'] == old
        assert result['body_norm_numerator'] == radius
        assert result['body_norm_denominator'] == 3
        assert result['count'] == expected[key], (key, result, expected[key])
        result['compact_class'] = int(oldnew[old])
        return result

    counts = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        for result in pool.map(recount, shells):
            counts.append(result)
            print('Fresh shell:', result['class_id'], result['count'], 'PASS', flush=True)

    tiny = np.load(ROOT / 'data/population_table.npz', allow_pickle=False)['tiny']
    angles = verify_older_43(compact, tiny)
    assert angles['older_seed_angles'] == 'PASS'
    assert sum(map(int, weights)) == 2878155
    sources = [ROOT / name for name in
               ('certify43.py', 'verify.py', 'final.py', 'bounds.json',
                'data/old43_seed.npz', 'data/e8_geometry.npz', 'data/population_table.npz',
                'population_replay/geometry.py', 'population_replay/verify_masses.py',
                'population_replay/old_e8_first_shells.py', 'population_replay/exact_count.cpp')]
    sources.extend(path for path in (ROOT / 'population_replay/data').rglob('*') if path.is_file())
    report = dict(
        status='PASS', dimension=43, points=2878155, families=113,
        head_rank=3, angle_arithmetic='arbitrary-size integers',
        family_pairs=113 * 114 // 2,
        large_shells='fresh exhaustive integer enumeration', counts=counts,
        smaller_shells='exact 11-design identities: 41 families of 9180 and 59 of 495',
        inherited_modular_population_lemma_required=False,
        published_inputs=['P48n minimum squared norm 6', '52416000 minimal vectors',
                          'minimal shell is an 11-design'],
        larger_working_total=2911990, larger_certificate='MISSING; NOT CERTIFIED',
        inputs_sha256={path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in sorted(sources)})
    (work / 'certificate43.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('PASS: 43D, 2,878,155 points; larger 2,911,990 certificate still missing.', flush=True)


def main():
    if not __debug__:
        raise SystemExit('Run without -O: verification assertions must be enabled.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--output', type=Path, help='New directory for replay outputs; otherwise temporary')
    args = parser.parse_args()
    assert args.jobs > 0
    verify_manifest()
    if args.output:
        args.output.mkdir(parents=True, exist_ok=False)
        certify(args.output.resolve(), args.jobs)
    else:
        with tempfile.TemporaryDirectory(prefix='kissing43_certificate_') as directory:
            certify(Path(directory), args.jobs)


if __name__ == '__main__':
    main()
