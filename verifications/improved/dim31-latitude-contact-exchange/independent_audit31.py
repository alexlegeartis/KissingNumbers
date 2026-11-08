"""Separate exact full-pair audit of the 8 November 31-dimensional candidate.

Run with Python and NumPy, g++ and AVX-512F/BW (or --scalar elsewhere):
  taskset -c 7 env OPENBLAS_NUM_THREADS=1 python independent_audit31.py SOURCE.npz TARGET.npz
No discovery or incumbent geometry-checker implementation is imported.
The optional reuse flags consume this audit's completed full C++ scan and
its exact binary input. All other checks are still repeated.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import subprocess
import tempfile
import time

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'data/baseline31_251196.npz'
TARGET = HERE / 'data/coordinates31_251200.npz'
EXPECTED_SOURCE = '02cddb863829f7298a1b80779603f3963f0c0db3eef0ff2a9050abe601e56127'
EXPECTED_TARGET = 'b75a8e9d958017fe03926265f5ee9d42f24b16c9ec848414b4a4e7eba0d329ee'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path, expected):
    require(sha(path) == expected, 'Frozen file hash: ' + str(path))
    with np.load(path, allow_pickle=False) as archive:
        require(set(archive.files) == {'points', 'squared_norm', 'inner_product_limit',
                                      'coordinate_denominator'}, 'Archive fields')
        for key, value in [('squared_norm', 128), ('inner_product_limit', 64),
                           ('coordinate_denominator', 2)]:
            require(archive[key].shape == () and archive[key].dtype.kind in 'iu'
                    and int(archive[key]) == value, 'Scalar integer metadata: ' + key)
        points = archive['points'].copy()
    require(points.ndim == 2 and points.shape[1] == 31 and points.dtype == np.int16,
            'Exact stored coordinate type and dimension')
    require(np.all((points >= -8) & (points <= 8)), 'Coefficient range')
    require(np.all(np.sum(points.astype(np.int64)**2, axis=1) == 128), 'Integer norms')
    keys = [tuple(map(int, row)) for row in points]
    require(len(set(keys)) == len(keys), 'Distinct complete coordinates')
    return points, keys


def int_categories(left, right, same=False):
    # Independent direct-int64 dot products. Bounded by 31*8*8 = 1984.
    matrix = left.astype(np.int64) @ right.astype(np.int64).T
    values = matrix[np.triu_indices(len(left), 1)] if same else matrix.ravel()
    require(len(values) > 0 and int(values.max()) <= 64, 'Independent integer angle')
    return dict(pairs=int(values.size), contacts=int(np.sum(values == 64)),
                violations=int(np.sum(values > 64)), minimum=int(values.min()),
                maximum=int(values.max()))


def pad(points, path):
    np.pad(points, ((0, 0), (0, 1))).astype('<i2').tofile(path)


def run(binary, data, count, retained, check=True):
    return subprocess.run([str(binary), str(data), str(count), str(retained)],
                          text=True, capture_output=True, check=check, cwd=binary.parent)


def controls(binary, directory, target):
    sample = target[np.linspace(0, len(target)-1, 128, dtype=int)]
    data = directory / 'sample.bin'
    pad(sample, data)
    actual = json.loads(run(binary, data, len(sample), len(sample)).stdout)
    require(actual['categories'][0] == int_categories(sample, sample, same=True),
            'SIMD and independent int64 sample results agree')
    # Same norm, distinct, dot 112: a real angle violation, not a norm failure.
    v = np.zeros((2, 31), dtype=np.int16)
    v[0, :2] = [8, 8]
    v[1, :6] = [8, 6, 4, 2, 2, 2]
    pad(v, data)
    result = run(binary, data, 2, 2, check=False)
    require(result.returncode == 1 and json.loads(result.stdout)['categories'][0]
            ['maximum'] == 112, 'Angle-violation negative control rejected')
    v[1] = v[0]
    pad(v, data)
    result = run(binary, data, 2, 2, check=False)
    require(result.returncode == 1 and json.loads(result.stdout)['categories'][0]
            ['maximum'] == 128, 'Duplicate negative control rejected')
    v[1, 0] = 7
    pad(v, data)
    result = run(binary, data, 2, 2, check=False)
    require(result.returncode != 0 and 'norm' in result.stderr,
            'Wrong-norm negative control rejected')
    return {'sample_unordered_pairs': 8128, 'simd_vs_int64': 'PASS',
            'same_norm_distinct_dot112_rejected': True,
            'duplicate_dot128_rejected': True, 'wrong_norm_rejected': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, nargs='?', default=SOURCE)
    parser.add_argument('target', type=Path, nargs='?', default=TARGET)
    parser.add_argument('--outdir', type=Path, default=HERE)
    parser.add_argument('--repository-root', type=Path)
    parser.add_argument('--reuse-full-pairs', type=Path)
    parser.add_argument('--reuse-input-bin', type=Path)
    parser.add_argument('--scalar', action='store_true')
    args = parser.parse_args()
    require(bool(args.reuse_full_pairs) == bool(args.reuse_input_bin), 'Both reuse paths required')
    args.outdir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    source, source_keys = load(args.source, EXPECTED_SOURCE)
    target, target_keys = load(args.target, EXPECTED_TARGET)
    require(len(source) == 251196 and len(target) == 251200, 'Recorded point counts')
    source_set, target_set = set(source_keys), set(target_keys)
    removed_ids = [i for i, row in enumerate(source_keys) if row not in target_set]
    added_ids = [i for i, row in enumerate(target_keys) if row not in source_set]
    retained_ids = [i for i, row in enumerate(source_keys) if row in target_set]
    require(len(removed_ids) == 4 and len(added_ids) == 8 and len(retained_ids) == 251192,
            'Set-derived 8-add / 4-remove delta')
    kept, removed, added = source[retained_ids], source[removed_ids], target[added_ids]
    require(np.array_equal(target, np.concatenate([kept, added])),
            'Entire candidate is retained source in order followed by exact additions')
    require(added_ids == list(range(len(kept), len(target))), 'C++ pair-category split')
    a, b, c, d = np.zeros((4, 31), dtype=np.int16)
    a[[0, 2, 6, 23]] = [-4, 4, 4, -4]
    b[[8, 9, 13, 14, 15, 17, 18, 22]] = [-2, 2, 2, 2, -2, 2, 2, 2]
    c[[24, 25]] = [-4, 4]
    d[[27, 29]] = [4, -4]
    blocks = np.array([a, b, c, d], dtype=np.int64)
    require(np.array_equal(blocks @ blocks.T, np.diag([64, 32, 32, 32])),
            'Four orthogonal algebraic blocks and their squared norms')
    formula_added = {tuple(s*a+t*b+u*c) for s,t,u in itertools.product((-1, 1), repeat=3)}
    formula_removed = {tuple(s*a+t*(c+d)) for s,t in itertools.product((-1, 1), repeat=2)}
    require(formula_added == set(map(tuple, added)), 'Eight-point formula')
    require(formula_removed == set(map(tuple, removed)), 'Four-point formula')
    all_old_dots = added.astype(np.int64) @ source.astype(np.int64).T
    blockers = [np.flatnonzero(row > 64).tolist() for row in all_old_dots]
    require(set(itertools.chain.from_iterable(blockers)) == set(removed_ids),
            'Deleted points exactly equal the union of inserted-point source blockers')
    repository_provenance = None
    if args.repository_root:
        # Optional source review is separate from the self-contained geometry proof.
        root = args.repository_root.resolve()
        original_source = root / 'verifications/improved/dim27-31-algebraic-completion/data/coordinates31_251196.npz'
        package = original_source.parents[1]
        manifest = package / 'data/SHA256SUMS.txt'
        manifest_files = 0
        for line in manifest.read_text().splitlines():
            if line:
                digest, relative = line.split('  ', 1)
                require(sha(package / relative) == digest, 'Source package hash: ' + relative)
                manifest_files += 1
        head = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
        source_relative = str(original_source.relative_to(root))
        git_bytes = subprocess.check_output(['git', '-C', str(root), 'show', 'HEAD:' + source_relative])
        require(hashlib.sha256(git_bytes).hexdigest() == EXPECTED_SOURCE, 'Current HEAD source bytes')
        provenance = json.loads((package / 'data/provenance.json').read_text())
        source_record = next(r for r in provenance['records'] if r['dimension'] == 31)
        require(source_record['total'] == 251196 and source_record['certificate_sha256'] == EXPECTED_SOURCE,
                'Current repository baseline provenance')
        repository_provenance = {'source': source_relative, 'git_head': head,
                                 'source_matches_HEAD_blob': True,
                                 'baseline_package_manifest_files_checked': manifest_files}

    cpp = HERE / 'independent_integer_audit.cpp'
    full_path = args.outdir / 'independent_full_pairs.json'
    with tempfile.TemporaryDirectory(prefix='independent_audit31_') as temporary:
        directory = Path(temporary)
        binary, data = directory / 'all_pairs', directory / 'points.bin'
        architecture_flags = [] if args.scalar else ['-mavx512f', '-mavx512bw']
        subprocess.run(['g++', '-O3', *architecture_flags, '-Wall', '-Wextra', '-Werror',
                        str(cpp), '-o', str(binary)], check=True)
        pad(target, data)
        padded_hash = sha(data)
        negative_controls = controls(binary, directory, target)
        if args.reuse_full_pairs:
            # First audit scan's frozen input remains available for byte comparison.
            require(args.reuse_input_bin.read_bytes() == data.read_bytes(),
                    'Previously scanned binary exactly matches frozen target bytes')
            full = json.loads(args.reuse_full_pairs.read_text())
        else:
            result = run(binary, data, len(target), len(kept))
            full_path.write_text(result.stdout)
            full = json.loads(result.stdout)
    require(full['rows'] == 251200 and full['retained'] == 251192, 'Full scan dimensions')
    rr, ra, aa = full['categories']
    expected_counts = [len(kept)*(len(kept)-1)//2, len(kept)*len(added),
                       len(added)*(len(added)-1)//2]
    require([r['pairs'] for r in [rr, ra, aa]] == expected_counts, 'All pair categories covered')
    require(all(r['violations'] == 0 and r['maximum'] <= 64 for r in [rr, ra, aa]),
            'Every candidate pair passes exact integer test')
    require(ra == int_categories(kept, added), 'Independent retained/new dot products')
    require(aa == int_categories(added, added, same=True), 'Independent new/new dot products')
    rd, dd = int_categories(kept, removed), int_categories(removed, removed, same=True)
    source_pairs = rr['pairs'] + rd['pairs'] + dd['pairs']
    target_pairs = sum(expected_counts)
    require(source_pairs == len(source)*(len(source)-1)//2, 'Complete source pair coverage')
    require(target_pairs == len(target)*(len(target)-1)//2, 'Complete candidate pair coverage')
    result = {
        'status': 'PASS', 'dimension': 31, 'baseline_points': 251196, 'candidate_points': 251200,
        'gain': 4, 'removed': 4, 'added': 8,
        'source_sha256': EXPECTED_SOURCE, 'candidate_sha256': EXPECTED_TARGET,
        'candidate_raw_points_sha256': hashlib.sha256(target.astype('<i2').tobytes()).hexdigest(),
        'padded_scan_input_sha256': padded_hash, 'repository_provenance': repository_provenance,
        'all_norms': 128, 'all_distinct': True, 'coordinate_denominator': 2,
        'stored_coordinate_abs_bound': 8, 'integer_accumulation_abs_bound': 1984,
        'candidate_unordered_pairs_checked': target_pairs,
        'baseline_unordered_pairs_checked': source_pairs,
        'pair_categories': {'retained_retained': rr, 'retained_added': ra, 'added_added': aa,
                            'retained_removed': rd, 'removed_removed': dd},
        'removed_source_indices_zero_based': removed_ids,
        'added_target_indices_zero_based': added_ids,
        'removed_points': removed.tolist(), 'added_points': added.tolist(),
        'algebraic_formula_checked': True, 'orthogonal_block_squared_norms': [64, 32, 32, 32],
        'added_point_source_blockers': blockers,
        'negative_controls': negative_controls, 'cpp_sha256': sha(cpp),
        'python_sha256': sha(__file__), 'full_integer_scan_seconds': full['seconds'],
        'audit_other_seconds': time.monotonic()-started,
        'integer_kernel': 'scalar int32' if args.scalar else 'AVX-512 signed int16 multiply / int32 sum',
        'arithmetic': 'All final unordered pairs individually enumerated using exact integer '
                      'multiplication and int32 summation. Independent changed-pair '
                      'checks use NumPy int64. No tolerance, antipodal compression or geometric '
                      'property of the source is assumed.',
        'mathematical_conclusion': '251200 distinct equal-norm vectors in R31 with all distinct '
                                   'pair inner products <= one half their squared norm; '
                                   'therefore kissing number K(31) >= 251200.',
        'limitations': ['This is a lower-bound certificate, not an optimality proof.',
                       'This checks the committed local 251196 baseline, not a fresh worldwide '
                       'literature or author-repository comparison.',
                       'Search-pool completeness and numerical optimizer optimality are not '
                       'needed and are not certified.']}
    (args.outdir / 'independent_audit31_report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
