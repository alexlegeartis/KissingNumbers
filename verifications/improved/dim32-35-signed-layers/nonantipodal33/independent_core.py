#!/usr/bin/env python3
"""Independent finite lower-bound audit. No producer modules are imported.

Reconstructs signed integer rows in R^33, checks every pair category by exact
arithmetic or a proved shape bound, and replays the actual old/new set delta.
Requires Python >=3.10 and NumPy. Run from any directory; --workspace relocates
the self-contained package. No network operations and no remote writes.
"""
import argparse
import copy
import hashlib
import itertools
import json
import os
from pathlib import Path
import time

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np

HERE = Path(__file__).resolve().parent
FILES = {
    'candidate': ('research/transfer_nonantipodal32_to33_20261107/certificate33_370392_explicit.json', 'a3bda380ae2f4aff5ab5429fd1ec5bae9b4e06f51b3b1525331519a325bbb53d'),
    'delta': ('research/transfer_nonantipodal32_to33_20261107/delta_370392.json', '27cde348c511f1d0a7b8c24b845bf7b23fe7aa3848d009f7203a201f7e39a4da'),
    'old_trade': ('research/publish43_2920901_repo/verifications/improved/dim32-35-signed-layers/root_trades/certificate33.json', '1c11d1585922b2e96ba70c05e342e2593bd6defc6ab9e17a7a9a5b38e847bb3b'),
    'classic': ('research/publish43_2920901_repo/verifications/improved/dim32-35-signed-layers/source/classic/certificate33.json', 'e87a39623fe881340c84c4bba8c89960b9e9eca90fe5452e169d64c0511a9751'),
    'source32': ('research/release_full_octad32_20261107/certificate32_352992_fixed_body_explicit.json', 'ab06963895d074d8dd6162bea2298796391dc876330ae997a5324adfb6638c5e'),
}


# Frozen source identifiers remain unchanged for exact delta provenance.
LOCAL_FILES = {'candidate': 'certificate33_370392_explicit.json', 'delta': 'delta_370392.json', 'old_trade': 'sources/certificate33_370376.json', 'classic': 'sources/certificate33_370368.json', 'source32': 'sources/certificate32_352992_fixed_body_explicit.json'}

def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked_json(data, expected):
    require(hashlib.sha256(data).hexdigest() == expected, 'SHA-256 mismatch')
    return json.loads(data)


def mask(x, width, label):
    require(type(x) is int and 0 <= x < (1 << width), label + ': invalid mask')


def positions(s, width=33):
    return [i for i in range(width) if (s >> i) & 1]


def negative_mask(s, local):
    return sum(1 << c for i, c in enumerate(positions(s)) if (local >> i) & 1)


def mapped_mask(s, cmap):
    return sum(1 << j for i, j in enumerate(cmap) if (s >> i) & 1)


def schema(d):
    n = d['dimension']
    require(n in (32, 33), 'dimension must be 32 or 33')
    require(type(d['points']) is int and d['points'] > 0, 'invalid point count')
    full = d['full_even_supports']
    require(len(full) == len(set(full)), 'duplicate full support')
    for s in full:
        mask(s, n, 'full support')
        require(s.bit_count() == 8, 'full support does not have weight eight')
    ek = 'exceptional_signed_points' if 'exceptional_signed_points' in d else 'exceptional_antipodal_lines'
    ep = d[ek]
    require(len(ep) == len(set(map(tuple, ep))), 'duplicate exceptional point or line')
    for s, m in ep:
        mask(s, n, 'exceptional support')
        mask(m, n, 'exceptional negative mask')
        require(s.bit_count() == 8 and m & ~s == 0, 'exceptional sign mask outside octad')
    removed = d.get('removed_root_lines', [])
    require(len(removed) == len(set(map(tuple, removed))), 'duplicate removed root line')
    for i, j, s in removed:
        require(type(i) is int and type(j) is int and 0 <= i < j < n and s in (-1, 1), 'bad root line')
    stars = d.get('stars', [])
    canonical = []
    for st in stars:
        a, bs = st['center'], st['leaves']
        require(type(a) is int and 0 <= a < n, 'bad star center')
        require(len(bs) == 4 and all(type(b) is int and 0 <= b < n for b in bs), 'bad star leaves')
        require(len(set(bs)) == 4 and a not in bs, 'star center/leaves not distinct')
        canonical.append((a, tuple(sorted(bs))))
    require(len(canonical) == len(set(canonical)), 'duplicate star')
    coords = d.get('dense_coordinates', list(range(32)))
    require(len(coords) == 32 and len(set(coords)) == 32, 'dense support must have exactly 32 coordinates')
    require(all(type(c) is int and 0 <= c < n for c in coords), 'invalid dense coordinate')
    mask(d.get('dense_coset_shift', 0), 32, 'dense coset shift')
    return ek


def dense(d):
    gens = d['dense_generator']
    require(len(gens) == 17, 'dense generator count must be 17')
    pivots = {}
    for g in gens:
        mask(g, 32, 'dense generator')
        w = g
        while w:
            p = w.bit_length() - 1
            if p in pivots:
                w ^= pivots[p]
            else:
                pivots[p] = w
                break
    require(len(pivots) == 17, 'dense code has deficient rank')
    words = [0]
    for g in gens:
        words += [w ^ g for w in words]
    require(len(set(words)) == 131072, 'dense code size is not 131072')
    hist = {}
    for w in words:
        hist[w.bit_count()] = hist.get(w.bit_count(), 0) + 1
    distance = min(k for k in hist if k)
    require(distance >= 8, 'dense minimum distance below eight')
    coords = d.get('dense_coordinates', list(range(32)))
    require(len(coords) == 32 and len(set(coords)) == 32, 'dense support must have exactly 32 coordinates')
    shift = d.get('dense_coset_shift', 0)
    vals = np.asarray([w ^ shift for w in words], dtype=np.uint64)
    result = np.zeros((len(words), 33), dtype=np.int8)
    result[:, coords] = 1 - 2 * ((vals[:, None] >> np.arange(32, dtype=np.uint64)) & 1).astype(np.int8)
    return result, {'binary_rank': len(pivots), 'size': len(words), 'minimum_distance': distance,
                    'weight_enumerator': hist, 'actual_coordinate_support': coords,
                    'coset_shift': shift, 'pair_maximum': 32 - 2 * distance,
                    'antipodal': (2**32 - 1) in set(words)}


def sparse(d, cmap=None):
    """Canonical descriptors: (kind, support, negative mask, heavy coordinate).

    kind 0: magnitude two octad; kind 1: magnitude four root;
    kind 2: four leaves of magnitude two and center of magnitude four.
    Source IDs are reconstructed, never imported from producer code.
    """
    ek = schema(d)
    n = d['dimension']
    cmap = list(range(n)) if cmap is None else cmap
    require(len(cmap) == n and len(set(cmap)) == n and all(0 <= x < 33 for x in cmap), 'invalid sparse coordinate map')
    result = {}
    categories = {'full': [], 'explicit': [], 'root': [], 'star': []}

    def add(desc, sid, cat):
        require(desc not in result, 'duplicate reconstructed sparse row')
        result[desc] = sid
        categories[cat].append(desc)

    for index, s in enumerate(d['full_even_supports']):
        for local in range(256):
            if local.bit_count() % 2 == 0:
                neg = negative_mask(s, local)
                add((0, mapped_mask(s, cmap), mapped_mask(neg, cmap), 0), ['full', index, local], 'full')
    for index, (s, m) in enumerate(d[ek]):
        if ek == 'exceptional_signed_points':
            add((0, mapped_mask(s, cmap), mapped_mask(m, cmap), 0), ['exceptional_point', index], 'explicit')
        else:
            for sign in (1, -1):
                mm = m if sign == 1 else s ^ m
                add((0, mapped_mask(s, cmap), mapped_mask(mm, cmap), 0), ['exceptional_line', index, sign], 'explicit')
    removed = set(map(tuple, d.get('removed_root_lines', [])))
    for i, j in itertools.combinations(range(n), 2):
        for a, b in itertools.product((1, -1), repeat=2):
            if (i, j, a*b) not in removed:
                ss = (1 << cmap[i]) | (1 << cmap[j])
                mm = (0 if a == 1 else 1 << cmap[i]) | (0 if b == 1 else 1 << cmap[j])
                add((1, ss, mm, 0), ['root', i, j, a, b], 'root')
    for index, st in enumerate(d.get('stars', [])):
        center, leaves = st['center'], st['leaves']
        ss = sum(1 << cmap[x] for x in [center] + leaves)
        for local in range(16):
            if local.bit_count() % 2 == 0:
                leafneg = sum(1 << cmap[x] for k, x in enumerate(leaves) if (local >> k) & 1)
                for sign in (1, -1):
                    mm = leafneg | (0 if sign == 1 else 1 << cmap[center])
                    add((2, ss, mm, cmap[center]), ['star', index, sign, local], 'star')
    return result, categories


def row(desc):
    kind, supp, neg, center = desc
    z = np.zeros(33, dtype=np.int8)
    for c in positions(supp):
        magnitude = 4 if kind == 1 or (kind == 2 and c == center) else 2
        z[c] = -magnitude if (neg >> c) & 1 else magnitude
    return z


def rows(descs):
    return np.asarray([row(d) for d in descs], dtype=np.int8).reshape((-1, 33))


def reconstruct(d, cmap=None):
    schema(d)
    dd, dproof = dense(d)
    sm, cats = sparse(d, cmap)
    arrays = {'dense': dd}
    arrays.update({c: rows(ds) for c, ds in cats.items()})
    allrows = np.concatenate(list(arrays.values()))
    require(len(allrows) == d['points'], 'claimed count differs from materialized count')
    require(np.all(np.sum(allrows.astype(np.int64)**2, axis=1) == 32), 'row squared norm is not 32')
    strings = {r.tobytes() for r in allrows}
    require(len(strings) == len(allrows), 'global duplicate row')
    require(set(np.flatnonzero(np.any(dd != 0, axis=0)).tolist()) == set(dproof['actual_coordinate_support']), 'dense coordinate support inconsistent')
    return {'data': d, 'sparse': sm, 'descriptors': cats, 'arrays': arrays,
            'rows': allrows, 'rowset': strings, 'dense_proof': dproof}


def pair_products(A, B=None, label='pairs'):
    self_pairs = B is None
    B = A if B is None else B
    A, B = A.astype(np.int64), B.astype(np.int64)
    maximum, minimum, count = -1000, 1000, 0
    witness = None
    for begin in range(0, len(A), 128):
        block = A[begin:begin+128] @ B.T
        if self_pairs:
            for ii in range(len(block)):
                vals = block[ii, begin+ii+1:]
                if len(vals):
                    mx, mn = int(vals.max()), int(vals.min())
                    if mx > maximum:
                        witness = [begin+ii, begin+ii+1+int(vals.argmax()), mx]
                    maximum, minimum = max(maximum, mx), min(minimum, mn)
                    count += len(vals)
                    require(mx <= 16, f'{label}: forbidden ordinary product {mx}')
        elif block.size:
            mx, mn = int(block.max()), int(block.min())
            if mx > maximum:
                ii, jj = np.unravel_index(block.argmax(), block.shape)
                witness = [begin+int(ii), int(jj), mx]
            maximum, minimum = max(maximum, mx), min(minimum, mn)
            count += block.size
            require(mx <= 16, f'{label}: forbidden ordinary product {mx}')
    return {'ordinary_maximum': maximum, 'ordinary_minimum': minimum,
            'pairs_checked': count, 'maximum_witness': witness, 'arithmetic': 'int64 coordinates'}


def full_geometry(full):
    maxoverlap = max(((s&t).bit_count() for s, t in itertools.combinations(full, 2)), default=0)
    require(maxoverlap <= 4, 'full/full support intersection exceeds four')
    patterns = [p for p in range(256) if p.bit_count() % 2 == 0]
    minhd = min((s^t).bit_count() for s, t in itertools.combinations(patterns, 2))
    require(len(patterns) == 128 and minhd >= 2, 'full/full parity distance invalid')
    return {'support_pairs_checked': len(full)*(len(full)-1)//2,
            'maximum_distinct_support_overlap': maxoverlap,
            'same_support_patterns': len(patterns), 'same_support_minimum_hamming_distance': minhd,
            'distinct_support_product_bound': 4*maxoverlap,
            'same_support_product_bound': 32-8*minhd}


def explicit_geometry(descs, full):
    mx, mn = -1000, 1000
    for x, y in itertools.combinations(descs, 2):
        _, s, a, _ = x
        _, t, b, _ = y
        dot = 4*((s&t).bit_count() - 2*((a^b)&s&t).bit_count())
        mx, mn = max(mx, dot), min(mn, dot)
        require(dot <= 16, f'explicit/explicit: forbidden ordinary product {dot}')
    supps = {d[1] for d in descs}
    overlap = max(((s&t).bit_count() for s in supps for t in full), default=0)
    require(overlap <= 4, 'explicit/full overlap exceeds four')
    return {'scalar_pairs_checked': len(descs)*(len(descs)-1)//2,
            'scalar_ordinary_maximum': mx, 'scalar_ordinary_minimum': mn,
            'explicit_support_count': len(supps),
            'explicit_full_support_pairs_checked': len(supps)*len(full),
            'explicit_full_maximum_overlap': overlap,
            'explicit_full_product_bound': 4*overlap}


def rank_mod101(A):
    pivots = {}
    for source in A:
        w = source.astype(np.int64) % 101
        for col, p in sorted(pivots.items()):
            w = (w - w[col]*p) % 101
        nz = np.flatnonzero(w)
        if len(nz):
            c = int(nz[0])
            pivots[c] = (w*pow(int(w[c]), -1, 101)) % 101
            if len(pivots) == 33:
                break
    return len(pivots)


def geometry(rec, label):
    a, ds = rec['arrays'], rec['descriptors']
    # Distinct sparse magnitude patterns certify all category shape premises.
    require(np.all(np.sum(np.abs(a['root']) == 4, axis=1) == 2), 'root shape invalid')
    for cat in ('full', 'explicit'):
        require(np.all(np.sum(np.abs(a[cat]) == 2, axis=1) == 8), 'octad shape invalid')
    full = [d[1] for d in ds['full'][::128]]
    proof = {'counts': {k: len(v) for k, v in a.items()}, 'total': len(rec['rows']),
             'dense': rec['dense_proof'], 'full_full': full_geometry(full)}
    proof['explicit_masks'] = explicit_geometry(ds['explicit'], full)
    proof['explicit_coordinates'] = pair_products(a['explicit'], label='explicit/explicit')
    require(proof['explicit_masks']['scalar_ordinary_maximum'] == proof['explicit_coordinates']['ordinary_maximum'], 'scalar/coordinate exceptional maxima disagree')
    require(proof['explicit_masks']['scalar_ordinary_minimum'] == proof['explicit_coordinates']['ordinary_minimum'], 'scalar/coordinate exceptional minima disagree')
    proof['root_root'] = pair_products(a['root'], label='root/root')
    proof['automatic_category_bounds'] = {
        'dense_full': 'Each octad has eight entries of absolute value 2; dense entries have absolute value at most 1, so product <= 16.',
        'dense_explicit': 'Same coordinatewise magnitude bound: product <= 16.',
        'dense_root': 'Two root entries of absolute value 4 against entries of absolute value at most 1: product <= 8.',
        'full_root': 'Only two root coordinates contribute, each by at most 8: product <= 16.',
        'explicit_root': 'Only two root coordinates contribute, each by at most 8: product <= 16.'}
    if len(a['star']):
        require(np.all(np.sum(np.abs(a['star']) == 4, axis=1) == 1) and np.all(np.sum(np.abs(a['star']) == 2, axis=1) == 4), 'star magnitude shape invalid')
        proof['star_categories'] = {k: pair_products(a['star'], v, 'star/'+k) for k, v in a.items() if k != 'star'}
        proof['star_categories']['star'] = pair_products(a['star'], label='star/star')
    proof['rank_mod101'] = rank_mod101(a['root'])
    proof['missing_antipodes'] = sum((-r).tobytes() not in rec['rowset'] for r in rec['rows'])
    proof['row_sha256'] = hashlib.sha256(rec['rows'].tobytes()).hexdigest()
    print(label + ': geometry PASS ' + str(proof['total']), flush=True)
    return proof


def check_delta(candidate, old, left, delta):
    c, o, l = candidate['sparse'], old['sparse'], left['sparse']
    cs, os_, ls = set(c), set(o), set(l)
    added, removed = cs-os_, os_-cs
    require(len(added) == 1648 and len(removed) == 1632, 'wrong exact delta sizes')
    require(all(d[0] == 0 for d in added|removed), 'delta changes a root or star')
    require(np.array_equal(candidate['arrays']['dense'], old['arrays']['dense']), 'dense rows changed')
    for cat in ('root', 'star'):
        require(set(candidate['descriptors'][cat]) == set(old['descriptors'][cat]), cat + ' rows changed')
    require(delta['coordinate_map'] == list(range(1, 33)), 'delta coordinate map differs')
    require(delta['certificate_sha256'] == FILES['candidate'][1], 'delta certificate hash differs')
    for source in delta['sources']:
        require((source['file'], source['sha256']) in FILES.values(), 'delta source identifier differs')
    require(len(delta['sources']) == 3 and len({x['file'] for x in delta['sources']}) == 3, 'delta sources missing or repeated')
    leftonly, rightonly = sorted(ls-os_), sorted(os_-ls)
    recorded = {}
    for kind, expected, only, source in [('added', added, leftonly, l), ('removed', removed, rightonly, o)]:
        seen = set()
        for graph_id, source_id, raw in delta[kind]:
            desc = tuple(raw)
            require(type(graph_id) is int and 0 <= graph_id < len(only), 'bad graph-side ID')
            require(desc == only[graph_id], 'graph-side ID does not reproduce descriptor')
            require(desc in source and source[desc] == source_id, 'original source ID does not reproduce descriptor')
            require(desc not in seen, 'duplicate delta point')
            seen.add(desc)
        require(seen == expected, kind + ' record differs from exact set difference')
        recorded[kind] = len(seen)
    require(cs == (os_-removed)|added, 'exact sparse set reconstruction differs')
    newrows, oldrows = candidate['rowset'], old['rowset']
    require(len(newrows-oldrows) == 1648 and len(oldrows-newrows) == 1632, 'coordinate delta differs')
    return dict(recorded, net_gain=len(newrows)-len(oldrows), old_total=len(oldrows), new_total=len(newrows),
                retained_roots=len(candidate['arrays']['root']), retained_stars=len(candidate['arrays']['star']),
                unchanged_dense=len(candidate['arrays']['dense']),
                source_common=131072+len(os_&ls), source_left_only=len(leftonly), source_right_only=len(rightonly),
                source_id_checks=len(added)+len(removed),
                scope='Exact geometry, source identities, set delta and union counts only; no matching/cover or union-optimality audit.')


def semantic_controls(d, rec, original_bytes):
    outcomes = []
    def reject(name, fn, expected=None):
        try:
            fn()
        except ValueError as e:
            require(expected is None or expected in str(e), name + ' failed for an unintended reason: ' + str(e))
            outcomes.append({'name': name, 'rejected': True, 'reason': str(e)})
            return
        raise ValueError(name + ': corruption was incorrectly accepted')
    def mutate(fn):
        x = copy.deepcopy(d)
        fn(x)
        return x
    reject('corrupted frozen byte hash', lambda: checked_json(original_bytes+b' ', FILES['candidate'][1]), 'SHA-256')
    x = mutate(lambda x: x['full_even_supports'].append(x['full_even_supports'][0]))
    reject('duplicate full support', lambda: schema(x), 'duplicate full')
    x = mutate(lambda x: x['exceptional_signed_points'].append(x['exceptional_signed_points'][0]))
    reject('duplicate explicit signed point', lambda: schema(x), 'duplicate exceptional')
    x = mutate(lambda x: x['stars'].append(x['stars'][0]))
    reject('duplicate star family', lambda: schema(x), 'duplicate star')
    x = mutate(lambda x: x.__setitem__('points', x['points']+1))
    reject('incorrect total count', lambda: reconstruct(x), 'claimed count')
    x = mutate(lambda x: x['dense_generator'].__setitem__(0, x['dense_generator'][1]))
    reject('deficient dense code rank', lambda: dense(x), 'deficient rank')
    x = mutate(lambda x: x['dense_generator'].__setitem__(0, 1))
    reject('dense code minimum distance one', lambda: dense(x), 'minimum distance')
    x = mutate(lambda x: x['dense_coordinates'].append(32))
    reject('dense support erroneously all 33 coordinates', lambda: schema(x), 'exactly 32')
    x = mutate(lambda x: x['dense_coordinates'].__setitem__(0, 1))
    reject('duplicate dense coordinate', lambda: schema(x), 'exactly 32')
    x = mutate(lambda x: x['full_even_supports'].__setitem__(0, x['full_even_supports'][0] | (1 << 33)))
    reject('out-of-range full support bit', lambda: schema(x), 'invalid mask')
    x = mutate(lambda x: x['exceptional_signed_points'][0].__setitem__(1, 1 << 32))
    reject('negative sign bit outside its exceptional support', lambda: schema(x), 'outside octad')
    x = mutate(lambda x: x['stars'][0]['leaves'].__setitem__(0, x['stars'][0]['center']))
    reject('star leaf duplicates its center', lambda: schema(x), 'not distinct')
    # Semantic geometry corruption with valid masks and norms, not a hash check.
    rootline = d['removed_root_lines'][0]
    i, j, s = rootline
    restored = np.zeros((1, 33), dtype=np.int8)
    restored[0, i], restored[0, j] = 4, 4*s
    reject('restore a forbidden root', lambda: pair_products(restored, rec['arrays']['star'], 'restored root/star'), '24')
    badstar = rec['arrays']['star'][0].copy()
    badstar[d['stars'][0]['leaves'][0]] *= -1
    reject('odd star leaf-sign parity at unchanged norm', lambda: pair_products(badstar[None, :], rec['arrays']['star'], 'bad star/star'), '24')
    badfull = rec['arrays']['full'][0].copy()
    badfull[positions(d['full_even_supports'][0])[0]] *= -1
    reject('odd full-octad sign parity at unchanged norm', lambda: pair_products(badfull[None, :], rec['arrays']['full'][:128], 'bad full/full'), '24')
    ep = rec['descriptors']['explicit']
    eset = set(ep)
    witness = None
    for ii, a in enumerate(ep):
        anti = (a[0], a[1], a[1]^a[2], a[3])
        if anti in eset:
            continue
        for jj, b in enumerate(ep):
            dot = int(row(anti).astype(np.int64) @ row(b).astype(np.int64))
            if dot == 24:
                witness = {'original_explicit_index': ii, 'conflict_explicit_index': jj,
                           'original': list(a), 'missing_antipode': list(anti),
                           'conflict_point': list(b), 'ordinary_product_after_completion': dot,
                           'ordinary_product_before_completion': -dot}
                break
        if witness:
            break
    require(witness is not None, 'no required missing-antipode +24 witness found')
    anti = tuple(witness['missing_antipode'])
    reject('append a missing antipode', lambda: pair_products(row(anti)[None, :], rec['arrays']['explicit'], 'antipode/explicit'), '24')
    require(int(row(tuple(witness['original'])).astype(np.int64) @ row(tuple(witness['conflict_point'])).astype(np.int64)) == -24, 'negative-product witness not reproduced')
    return {'rejection_count': len(outcomes), 'results': outcomes, 'missing_antipode_witness': witness,
            'semantic_note': 'The valid configuration contains an ordinary product -24. An absolute-product test would reject it incorrectly; completing antipodes creates +24 and is invalid.'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workspace', type=Path, default=HERE)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--skip-controls', action='store_true')
    args = ap.parse_args()
    start = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)
    loaded, hashes = {}, {}
    for key, (path, expected) in FILES.items():
        data = (args.workspace/LOCAL_FILES[key]).read_bytes()
        loaded[key] = checked_json(data, expected)
        hashes[key] = {'file': path, 'sha256': expected}
    d, oldtrade, classic = loaded['candidate'], loaded['old_trade'], loaded['classic']
    require(d['format'] == 'signed_octad_star_explicit_v1', 'wrong candidate format')
    require(d['dimension'] == 33 and d['points'] == 370392, 'wrong frozen target')
    require(d['dense_coordinates'] == list(range(32)) and d['dense_coset_shift'] == 0, 'frozen dense placement changed')
    require(oldtrade['baseline_sha256'] == FILES['classic'][1], 'old trade baseline hash mismatch')
    require(oldtrade['baseline'] == 'source/classic/certificate33.json', 'unexpected old baseline path')
    old = copy.deepcopy(classic)
    old.update({k: oldtrade[k] for k in ('points', 'stars', 'removed_root_lines')})
    candidate = reconstruct(d)
    cproof = geometry(candidate, 'candidate')
    require(cproof['missing_antipodes'] == 372, 'missing-antipode count differs')
    require(cproof['rank_mod101'] == 33, 'candidate does not span dimension 33')
    require(cproof['counts'] == {'dense':131072, 'full':235136, 'explicit':2064, 'root':2056, 'star':64}, 'category counts differ')
    oldrec = reconstruct(old)
    oproof = geometry(oldrec, 'actual old33 control')
    classicrec = reconstruct(classic)
    # Direct old geometry checks all classic octad and dense pairs. Restoring
    # classic roots is safe by the root-root and magnitude bounds; nonetheless
    # replay its complete geometry too, with no reliance on a producer proof.
    bproof = geometry(classicrec, 'classic33 baseline')
    left = reconstruct(loaded['source32'], list(range(1, 33)))
    lproof = geometry(left, 'shifted source32 with old dense support')
    delta = check_delta(candidate, oldrec, left, loaded['delta'])
    result = {'status':'PASS', 'dimension':33, 'lower_bound':370392,
              'hashes':hashes, 'candidate':cproof, 'old_control':oproof,
              'classic_baseline':bproof, 'shifted_source32':lproof, 'delta':delta,
              'independence':'No producer proof, search, graph, matching, or verification module was imported or called.',
              'scope':'Direct finite lower-bound audit and exact +16 comparison with the actual old33 repository control. No worldwide-priority, global-optimality, or matching/cover claim.'}
    (args.output/'geometry_pass.json').write_text(json.dumps(result, indent=2)+'\n')
    print('MATHEMATICAL PASS: 370392 distinct norm-32 rows in R^33; ordinary products <=16; exact old-control gain +16.', flush=True)
    if not args.skip_controls:
        result['controls'] = semantic_controls(d, candidate, (args.workspace/LOCAL_FILES['candidate']).read_bytes())
    result['elapsed_seconds'] = time.monotonic()-start
    result['verifier_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.output/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print('AUDIT PASS including ' + str(result.get('controls', {}).get('rejection_count', 0)) + ' controls; report ' + str(args.output/'report.json'), flush=True)


if __name__ == '__main__':
    main()
