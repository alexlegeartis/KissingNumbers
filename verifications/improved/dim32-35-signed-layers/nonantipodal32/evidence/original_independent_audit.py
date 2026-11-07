#!/usr/bin/env python3
"""Standalone exact audit; imports no construction, search, or earlier checker code.

From the workspace root:
  taskset -c 5 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
    research/independent_audit_nonantipodal32_20261107/audit.py

Requires Python 3.10+ and NumPy. Certificate masks use bit j for coordinate j.
An explicit [support, negatives] entry is ONE point, never an antipodal line.
"""
import argparse
import collections
import copy
import hashlib
import itertools
import json
import os
import pathlib
import platform
import time

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
HERE = pathlib.Path(__file__).resolve().parent
SOURCES = {
    'fixed_body': ('research/release_full_octad32_20261107/certificate32_352992_fixed_body_explicit.json', 'ab06963895d074d8dd6162bea2298796391dc876330ae997a5324adfb6638c5e'),
    'discovery': ('research/release_full_octad32_20261107/certificate32_352992_explicit.json', 'a21b8ee4086cdbbc90651bcab378b888d0c3f10643aa0d63b5d208cafe471d1a'),
    'baseline': ('research/expand_signed_exceptional32_20261107/certificate32_352832.json', '69c3f368a4aae5956ca75a8569f761224148831e8a85edd5972a0b3b057af646'),
    'delta': ('research/release_full_octad32_20261107/fixed_body_delta.json', '4ee5f24fa852acda415028c4e4c6ec974b1225271428cdb858134eead4f7d37d'),
}

class Invalid(ValueError):
    pass


def require(truth, message):
    if not truth:
        raise Invalid(message)


def valid_mask(value):
    return type(value) is int and 0 <= value < (1 << 32)


def row_for(pair):
    support, negative = pair
    return [(-2 if negative & (1 << j) else 2) if support & (1 << j) else 0 for j in range(32)]


def point_product(a, b):
    shared = a[0] & b[0]
    return 4 * (shared.bit_count() - 2 * (shared & (a[1] ^ b[1])).bit_count())


def unpack_points(document):
    """Expand legacy lines only; the explicit format has no implicit points."""
    if document['format'] == 'signed_octad_kissing_explicit_v1':
        entries = document['exceptional_signed_points']
        antipodal = False
    elif document['format'] == 'signed_octad_kissing_v1':
        entries = document['exceptional_antipodal_lines']
        antipodal = True
    else:
        raise Invalid('certificate format')
    require(type(entries) is list, 'exceptional collection')
    result = []
    for entry in entries:
        require(type(entry) is list and len(entry) == 2, 'exceptional entry shape')
        support, negative = entry
        require(valid_mask(support) and support.bit_count() == 8, 'exceptional support mask')
        require(valid_mask(negative) and negative & ~support == 0, 'exceptional negative mask')
        result.append((support, negative))
        if antipodal:
            result.append((support, support ^ negative))
    require(len(result) == len(set(result)), 'duplicate exceptional point')
    return result


def code_words(generator):
    require(type(generator) is list and len(generator) == 17 and all(valid_mask(g) for g in generator), 'dense generator shape')
    pivots = {}
    for entry in generator:
        word = entry
        while word:
            pivot = word.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = word
                break
            word ^= pivots[pivot]
    require(len(pivots) == 17, 'dense generator rank')
    # Consecutive Gray labels change precisely the row indexed by v_2(i).
    words = [0]
    current = 0
    for i in range(1, 1 << 17):
        current ^= generator[(i & -i).bit_length() - 1]
        words.append(current)
    require(len(set(words)) == 1 << 17, 'dense enumeration distinctness')
    weights = collections.Counter(w.bit_count() for w in words)
    require(min(w for w in weights if w) >= 8, 'dense minimum distance')
    return words, dict(sorted(weights.items()))


def histogram(values):
    keys, amounts = np.unique(values, return_counts=True)
    return {int(k): int(v) for k, v in zip(keys, amounts)}


def complete_check(document, materialize=True):
    require(document.get('dimension') == 32 and type(document.get('dimension')) is int, 'dimension')
    points = unpack_points(document)
    full = document['full_even_supports']
    require(type(full) is list and all(valid_mask(s) and s.bit_count() == 8 for s in full), 'full support mask')
    require(len(full) == len(set(full)), 'duplicate full support')
    claimed = 131072 + 128 * len(full) + len(points) + 1984
    require(type(document.get('points')) is int and document['points'] == claimed, 'claimed total')
    words, weights = code_words(document['dense_generator'])
    full_intersections = collections.Counter()
    for i, support in enumerate(full):
        for previous in full[:i]:
            overlap = (support & previous).bit_count()
            full_intersections[overlap] += 1
            require(overlap <= 4, 'full/full overlap')
    exceptional_supports = sorted({s for s, n in points})
    mixed_intersections = collections.Counter()
    for support in exceptional_supports:
        for other in full:
            overlap = (support & other).bit_count()
            mixed_intersections[overlap] += 1
            require(overlap <= 4, 'full/exceptional overlap')
    # Two separately expressed arithmetic routes on EVERY explicit point pair.
    pair_hist = collections.Counter()
    negative_witness = None
    for i, point in enumerate(points):
        for j in range(i):
            product = point_product(point, points[j])
            require(product <= 16, 'exceptional ordinary product')
            pair_hist[product] += 1
            if product < -16 and product > -32 and negative_witness is None:
                negative_witness = {'i': i, 'j': j, 'point_i': list(point), 'point_j': list(points[j]), 'dot': product}
    exceptional = np.array([row_for(p) for p in points], dtype=np.int64).reshape((-1, 32))
    require(np.all(np.sum(exceptional * exceptional, axis=1) == 32), 'exceptional norms')
    require(len({r.tobytes() for r in exceptional}) == len(points), 'exceptional coordinate distinctness')
    gram = exceptional @ exceptional.T
    upper = gram[np.triu_indices(len(points), 1)]
    require(np.all(upper <= 16), 'exceptional coordinate ordinary product')
    require(histogram(upper) == dict(pair_hist), 'independent product histograms agree')
    # Reconstruct all 128 even sign patterns using a product alphabet, not masks.
    patterns = np.array([s for s in itertools.product((-2, 2), repeat=8) if sum(x < 0 for x in s) % 2 == 0], dtype=np.int64)
    require(patterns.shape == (128, 8), 'full parity layer count')
    require(len({r.tobytes() for r in patterns}) == 128, 'full parity layer distinctness')
    require(np.all(np.sum(patterns * patterns, axis=1) == 32), 'full parity layer norms')
    full_same = (patterns @ patterns.T)[np.triu_indices(128, 1)]
    require(int(full_same.max()) <= 16, 'full parity internal product')
    roots = []
    for a in range(32):
        for b in range(a + 1, 32):
            for x, y in itertools.product((-4, 4), repeat=2):
                row = [0] * 32
                row[a], row[b] = x, y
                roots.append(row)
    roots = np.array(roots, dtype=np.int64)
    require(roots.shape == (1984, 32), 'root count')
    require(len({r.tobytes() for r in roots}) == 1984, 'root distinctness')
    require(np.all(np.sum(roots * roots, axis=1) == 32), 'root norms')
    root_upper = (roots @ roots.T)[np.triu_indices(1984, 1)]
    require(int(root_upper.max()) <= 16, 'root/root product')
    point_set = set(points)
    unpaired = sorted((s, n) for s, n in points if (s, s ^ n) not in point_set)
    counts = {'dense': len(words), 'full_even': 128 * len(full), 'exceptional': len(points), 'roots': len(roots)}
    # All ten unordered family categories are covered. Shape bounds are pointwise,
    # so neither an untested lattice minimum nor a full gigantic Gram is needed.
    category_bounds = {
        'dense/dense': {'upper_bound': 32 - 2 * min(w for w in weights if w), 'reason': 'Binary linear differences range over the exhaustively enumerated nonzero codewords.'},
        'dense/full': {'upper_bound': 16, 'reason': 'Eight products of magnitude 2.'},
        'dense/exceptional': {'upper_bound': 16, 'reason': 'Eight products of magnitude 2, independently of signs or parity.'},
        'dense/root': {'upper_bound': 8, 'reason': 'Two products of magnitude 4.'},
        'full/full': {'upper_bound': 16, 'reason': 'Same support: exact 128-pattern Gram; different supports: every overlap at most 4, each contribution magnitude 4.'},
        'full/exceptional': {'upper_bound': 16, 'reason': 'Every support overlap checked at most 4, each contribution magnitude 4.'},
        'full/root': {'upper_bound': 16, 'reason': 'At most two products of magnitude 8.'},
        'exceptional/exceptional': {'upper_bound': int(upper.max()), 'reason': 'Every unordered pair checked by scalar signed masks and independently by int64 coordinate products.'},
        'exceptional/root': {'upper_bound': 16, 'reason': 'At most two products of magnitude 8.'},
        'root/root': {'upper_bound': int(root_upper.max()), 'reason': 'Every unordered pair checked in int64 coordinates.'},
    }
    require(max(v['upper_bound'] for v in category_bounds.values()) <= 16, 'all family bounds')
    report = {
        'format': document['format'], 'counts': counts, 'total': sum(counts.values()),
        'dense_generator_rank': 17, 'dense_words_enumerated': len(words), 'dense_minimum_distance': min(w for w in weights if w),
        'dense_weight_distribution': weights, 'full_support_count': len(full), 'exceptional_support_count': len(exceptional_supports),
        'full_full_support_intersection_histogram': dict(sorted(full_intersections.items())),
        'full_exceptional_support_intersection_histogram': dict(sorted(mixed_intersections.items())),
        'full_same_support_pattern_pair_histogram': histogram(full_same),
        'exceptional_pair_count_checked_twice': len(points) * (len(points) - 1) // 2,
        'exceptional_pair_product_histogram': dict(sorted(pair_hist.items())),
        'exceptional_signed_parity_histogram': dict(sorted(collections.Counter(n.bit_count() % 2 for s, n in points).items())),
        'exceptional_support_population_histogram': dict(sorted(collections.Counter(collections.Counter(s for s, n in points).values()).items())),
        'exceptional_paired_points': len(points) - len(unpaired), 'exceptional_antipodal_pairs': (len(points) - len(unpaired)) // 2,
        'exceptional_unpaired_points': len(unpaired), 'non_antipodal_negative_product_witness': negative_witness,
        'root_pair_count': len(roots) * (len(roots) - 1) // 2, 'root_pair_product_histogram': histogram(root_upper),
        'pair_category_bounds': category_bounds, 'squared_norm': 32, 'maximum_distinct_pair_product': 16,
        'minimum_squared_distance': 32, 'pair_scope': 'Complete proof by all ten pair categories, not a full-configuration brute-force Gram sweep.'
    }
    if materialize:
        dense = np.array([[(-1 if w & (1 << j) else 1) for j in range(32)] for w in words], dtype=np.int8)
        full_rows = np.zeros((128 * len(full), 32), dtype=np.int8)
        for i, support in enumerate(full):
            columns = [j for j in range(32) if support & (1 << j)]
            full_rows[128 * i:128 * (i + 1), columns] = patterns
        matrix = np.concatenate((dense, full_rows, exceptional.astype(np.int8), roots.astype(np.int8)))
        require(matrix.shape == (claimed, 32), 'materialized shape')
        require(np.all(np.sum(matrix.astype(np.int64) ** 2, axis=1) == 32), 'materialized norms')
        row_set = {r.tobytes() for r in matrix}
        distinct = len(row_set)
        require(distinct == claimed, 'materialized global distinctness')
        global_unpaired = sum((-r).tobytes() not in row_set for r in matrix)
        require(global_unpaired == len(unpaired), 'materialized global antipode count')
        report['materialized_rows'] = len(matrix)
        report['materialized_distinct_rows'] = distinct
        report['materialized_global_unpaired_points'] = global_unpaired
        report['coordinates_sha256_int8_rowmajor'] = hashlib.sha256(matrix.tobytes()).hexdigest()
        report['coordinate_order'] = 'Dense Gray-generator order; full support file order and lexicographic (-2,+2) even patterns; explicit file order (legacy line then antipode); roots lexicographic coordinates and (-4,+4) signs.'
    return report


def full_octad_points(document):
    result = set(unpack_points(document))
    for support in document['full_even_supports']:
        bits = [1 << j for j in range(32) if support & (1 << j)]
        for label in range(256):
            if label.bit_count() % 2 == 0:
                result.add((support, sum(bits[k] for k in range(8) if label & (1 << k))))
    return result


def audit_delta(baseline, fixed, discovery, advertised):
    require(baseline['dense_generator'] == fixed['dense_generator'] == discovery['dense_generator'], 'unchanged dense generator')
    require(baseline['full_even_supports'] == fixed['full_even_supports'], 'unchanged full-even body')
    old, new = set(unpack_points(baseline)), set(unpack_points(fixed))
    removed, added = sorted(old - new), sorted(new - old)
    unpaired = sorted((s, n) for s, n in new if (s, s ^ n) not in new)
    require(advertised['parent_sha256'] == SOURCES['baseline'][1], 'delta parent hash')
    require(advertised['candidate_sha256'] == SOURCES['fixed_body'][1], 'delta candidate hash')
    require(advertised['removed_points'] == [list(p) for p in removed], 'delta removed list')
    require(advertised['added_points'] == [list(p) for p in added], 'delta added list')
    require(advertised['unpaired_exceptional_points'] == [list(p) for p in unpaired], 'delta unpaired list')
    require(advertised['gain'] == len(added) - len(removed) == 160, 'delta gain')
    require((len(removed), len(added), len(old), len(new)) == (1316, 1476, 1664, 1824), 'expected delta counts')
    old_supports, new_supports = {s for s, n in old}, {s for s, n in new}
    require((len(old_supports & new_supports), len(new_supports - old_supports)) == (15, 14), 'support turnover')
    discovery_points, normalized_points = full_octad_points(discovery), full_octad_points(fixed)
    only_discovery, only_fixed = sorted(discovery_points - normalized_points), sorted(normalized_points - discovery_points)
    require(len(only_discovery) == len(only_fixed) == 128, 'discovery normalization population')
    changed_supports = sorted({s for s, n in only_discovery + only_fixed})
    require(len(changed_supports) == 1, 'discovery normalization support')
    require(all(n.bit_count() % 2 == 1 for s, n in only_discovery), 'discovery changed odd parity')
    require(all(n.bit_count() % 2 == 0 for s, n in only_fixed), 'fixed changed even parity')
    recovered = sorted(set(fixed['full_even_supports']) - set(discovery['full_even_supports']))
    require(len(recovered) == 3 and set(discovery['full_even_supports']) < set(fixed['full_even_supports']), 'three recovered full layers')
    return {'baseline_exceptional_points': len(old), 'fixed_exceptional_points': len(new),
        'removed_points': len(removed), 'added_points': len(added), 'retained_points': len(old & new), 'net_gain': len(new) - len(old),
        'old_exceptional_supports': len(old_supports), 'new_exceptional_supports': len(new_supports),
        'retained_exceptional_supports': len(old_supports & new_supports), 'added_exceptional_supports': len(new_supports - old_supports),
        'removed_exceptional_supports': len(old_supports - new_supports), 'dense_and_full_body_unchanged': True,
        'advertised_delta_matches_exactly': True, 'discovery_normalization': {
            'same_count_not_same_coordinate_set': True, 'three_recovered_supports': recovered,
            'odd_to_even_support': changed_supports[0], 'removed_odd_points': len(only_discovery), 'added_even_points': len(only_fixed)}}


def negative_controls(original):
    results = []
    def attempt(name, mutation, reason, witness=None):
        document = copy.deepcopy(original)
        mutation(document)
        try:
            complete_check(document, materialize=False)
        except Invalid as error:
            require(str(error) == reason, f'control {name}: unexpected rejection: {error}')
            results.append({'name': name, 'rejected': True, 'reason': str(error), **({'witness': witness} if witness else {})})
        else:
            raise Invalid(f'control unexpectedly accepted: {name}')
    attempt('duplicate explicit point', lambda d: d['exceptional_signed_points'].__setitem__(0, d['exceptional_signed_points'][1]), 'duplicate exceptional point')
    attempt('wrong declared point count', lambda d: d.__setitem__('points', d['points'] + 1), 'claimed total')
    attempt('zero dense generator row', lambda d: d['dense_generator'].__setitem__(0, 0), 'dense generator rank')
    attempt('weight-one dense generator', lambda d: d['dense_generator'].__setitem__(0, 1), 'dense minimum distance')
    attempt('out-of-range generator mask', lambda d: d['dense_generator'].__setitem__(0, 1 << 32), 'dense generator shape')
    attempt('invalid explicit support weight', lambda d: d['exceptional_signed_points'][0].__setitem__(0, 1), 'exceptional support mask')
    s, n = original['exceptional_signed_points'][0]
    outside = next(1 << j for j in range(32) if not s & (1 << j))
    attempt('explicit negative bit outside support', lambda d: d['exceptional_signed_points'][0].__setitem__(1, n | outside), 'exceptional negative mask')
    attempt('out-of-range explicit support', lambda d: d['exceptional_signed_points'][0].__setitem__(0, s | (1 << 32)), 'exceptional support mask')
    attempt('duplicate full layer', lambda d: d['full_even_supports'].__setitem__(0, d['full_even_supports'][1]), 'duplicate full support')
    f = original['full_even_supports'][0]
    moved = (f ^ (f & -f)) | next(1 << j for j in range(32) if not f & (1 << j))
    def add_bad_full(d):
        d['full_even_supports'].append(moved)
        d['points'] += 128
    attempt('new full support with overlap seven', add_bad_full, 'full/full overlap')
    def add_conflicting_full(d):
        d['full_even_supports'].append(s)
        d['points'] += 128
    attempt('explicit support reused as a full layer', add_conflicting_full, 'full/exceptional overlap')
    pts = unpack_points(original)
    point_set = set(pts)
    corruption = None
    for i, (support, negative) in enumerate(pts):
        for j in range(32):
            if not support & (1 << j):
                continue
            replacement = (support, negative ^ (1 << j))
            if replacement in point_set:
                continue
            offender = next((k for k, other in enumerate(pts) if k != i and point_product(replacement, other) > 16), None)
            if offender is not None:
                corruption = (i, replacement, offender)
                break
        if corruption is not None:
            break
    require(corruption is not None, 'find genuine incompatible corruption')
    i, replacement, offender = corruption
    witness = {'replaced_index': i, 'replacement': list(replacement), 'other_index': offender, 'other_point': list(pts[offender]), 'dot': point_product(replacement, pts[offender])}
    attempt('actual nonduplicate incompatible sign corruption', lambda d: d['exceptional_signed_points'].__setitem__(i, list(replacement)), 'exceptional ordinary product', witness)
    # Find a valid pair with dot -24, not merely the trivial antipodal -32.
    semantic = None
    for i, p in enumerate(pts):
        antipode = (p[0], p[0] ^ p[1])
        if antipode in point_set:
            continue
        for j, q in enumerate(pts):
            value = point_product(p, q)
            if -32 < value < -16:
                semantic = (p, q, antipode, value)
                break
        if semantic:
            break
    require(semantic is not None, 'non-antipodal semantic witness')
    p, q, antipode, product = semantic
    def add_antipode(d):
        d['exceptional_signed_points'].append(list(antipode))
        d['points'] += 1
    witness = {'valid_point': list(p), 'other_valid_point': list(q), 'original_dot': product, 'missing_antipode': list(antipode), 'new_dot': -product}
    attempt('adding a missing antipode is invalid', add_antipode, 'exceptional ordinary product', witness)
    results.append({'name': 'absolute-product rule would falsely reject a valid non-antipodal pair', 'valid_ordinary_product': product, 'incorrect_absolute_product': abs(product), 'threshold': 16, 'witness': witness, 'demonstrated': True})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=pathlib.Path, default=ROOT)
    parser.add_argument('--output', type=pathlib.Path, default=HERE / 'report.json')
    args = parser.parse_args()
    started = time.monotonic()
    source_docs, source_meta = {}, {}
    for name, (relative, expected_hash) in SOURCES.items():
        raw = (args.root / relative).read_bytes()
        actual = hashlib.sha256(raw).hexdigest()
        require(actual == expected_hash, f'{name} SHA256')
        source_docs[name] = json.loads(raw)
        source_meta[name] = {'path': relative, 'sha256': actual, 'bytes': len(raw)}
    reports = {}
    for name in ('fixed_body', 'discovery', 'baseline'):
        reports[name] = complete_check(source_docs[name])
        print(json.dumps({'stage': name, 'status': 'PASS', 'points': reports[name]['total'], 'seconds': round(time.monotonic() - started, 3)}), flush=True)
    require(reports['fixed_body']['counts'] == {'dense': 131072, 'full_even': 218112, 'exceptional': 1824, 'roots': 1984}, 'fixed expected populations')
    require(reports['discovery']['counts'] == {'dense': 131072, 'full_even': 217728, 'exceptional': 2208, 'roots': 1984}, 'discovery expected populations')
    require((reports['fixed_body']['exceptional_unpaired_points'], reports['fixed_body']['exceptional_antipodal_pairs']) == (380, 722), 'fixed pairing populations')
    require((reports['discovery']['exceptional_unpaired_points'], reports['discovery']['exceptional_antipodal_pairs']) == (380, 914), 'discovery pairing populations')
    delta = audit_delta(source_docs['baseline'], source_docs['fixed_body'], source_docs['discovery'], source_docs['delta'])
    controls = negative_controls(source_docs['fixed_body'])
    # Re-read originals after every control to confirm no frozen input changed.
    for name, (relative, expected_hash) in SOURCES.items():
        require(hashlib.sha256((args.root / relative).read_bytes()).hexdigest() == expected_hash, f'{name} input immutable')
    report = {'status': 'PASS', 'source_files': source_meta, 'audit_script_sha256': hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), 'certificates': reports, 'delta': delta, 'negative_controls': controls,
        'independence': 'New standalone implementation. No producer search/verifier functions or prior checker functions are imported or called. An earlier independent audit was read only for format context.',
        'scientific_scope': 'Both 352992-point configurations are proved valid kissing configurations in R^32; each is a constructive lower bound, not an optimality proof. Exact comparisons to the frozen 352832 certificate are established. No worldwide-priority claim.',
        'limitations': 'Pair-category proof rather than a 62-billion-pair full Gram sweep. NumPy int64 products and Python exact integers are used; arithmetic magnitudes are at most 32. No external theorem about a lattice minimum/population is assumed. Public-source freshness is recorded separately.',
        'runtime': {'python': platform.python_version(), 'numpy': np.__version__, 'cpu_affinity': sorted(os.sched_getaffinity(0)), 'thread_environment': {v: os.environ[v] for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')}},
        'elapsed_seconds': time.monotonic() - started}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'totals': {k: v['total'] for k, v in reports.items()}, 'gain': delta['net_gain'], 'negative_controls': len(controls), 'elapsed_seconds': report['elapsed_seconds'], 'report': str(args.output)}, indent=2), flush=True)

if __name__ == '__main__':
    main()
