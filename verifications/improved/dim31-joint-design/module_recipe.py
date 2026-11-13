"""Reconstruct the Leech--Sloane module and replay the exact 31D exchange.

The small construction is rebuilt from a Steiner quadruple system and two
affine binary codes. Every module pair, the supplied rational isometry, source
identity, and every pair involving an added point are checked with int64
arithmetic. The separate verify_geometry.py checks all final pairs without
assuming validity of the external source. No optimizer or external data is
needed by this checker.
"""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np

ROOT = Path(__file__).resolve().parent
NORMAL = (0, 2, 4, 5, 6, 10, 11, 23)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer_array(value, label, bound=16):
    value = np.asarray(value)
    require(value.dtype.kind in 'iu', label + ': integer storage')
    require(np.all(value >= -bound) and np.all(value <= bound),
            label + ': safe coefficient bounds')
    return value.astype(np.int64)


def row_keys(points):
    """Small-coordinate integer rows, encoded without native-endian ambiguity."""
    return [row.tobytes() for row in np.asarray(points, dtype='<i2')]


def canonical_hash(points):
    ordered = np.unique(points, axis=0)
    return hashlib.sha256(np.asarray(ordered, dtype='<i2').tobytes()).hexdigest()


def reconstruct(recipe):
    """Build all 2564 points in the weighted model, then map to Euclidean space."""
    blocks = [tuple(block) for block in recipe['transversal_blocks']]
    require(len(blocks) == 70 and len(set(blocks)) == 70,
            'Exactly 70 distinct transversal blocks')
    for block in blocks:
        require(len(block) == 4 and len(set(block)) == 4
                and all(type(j) is int and 0 <= j < 14 for j in block)
                and len({j // 2 for j in block}) == 4,
                'Each transversal uses four different coordinate pairs')
    pair_unions = [(2*i, 2*i+1, 2*j, 2*j+1)
                   for i, j in combinations(range(7), 2)]
    quads = [tuple(sorted(block)) for block in pair_unions + blocks]
    triples = Counter(triple for quad in quads
                      for triple in combinations(quad, 3))
    require(len(quads) == 91 and len(set(quads)) == 91,
            'Exactly 91 distinct quadruples')
    require(set(triples) == set(combinations(range(14), 3))
            and set(triples.values()) == {1},
            'Every 3-subset of 14 coordinates belongs to exactly one quadruple')

    rows = []
    for j in range(14):
        for sign in (-4, 4):
            row = [0] * 15
            row[j] = sign
            rows.append(row)
    for support in quads:
        for signs in product((-2, 2), repeat=4):
            row = [0] * 15
            for j, sign in zip(support, signs):
                row[j] = sign
            rows.append(row)

    layers = recipe['glue_layers']
    require(len(layers) == 2 and {layer['height'] for layer in layers} == {-1, 1},
            'Two opposite-height glue layers')
    for layer in layers:
        origin, basis = layer['origin'], layer['basis']
        require(type(origin) is int and 0 <= origin < 2**14
                and len(basis) == 9
                and all(type(b) is int and 0 < b < 2**14 for b in basis),
                'Affine binary code format')
        words = {origin}
        for generator in basis:
            words |= {word ^ generator for word in words}
        require(len(words) == 512 and layer['count'] == 512,
                'Each affine glue code has 512 distinct words')
        for word in sorted(words):
            rows.append([1 - 2*((word >> j) & 1) for j in range(14)]
                        + [layer['height']])
    for pair in range(7):
        for a, b, height in product((-2, 2), repeat=3):
            row = [0] * 15
            row[2*pair], row[2*pair+1], row[14] = a, b, height
            rows.append(row)
    weighted = np.asarray(rows, dtype=np.int64)
    require(weighted.shape == (2564, 15), 'Weighted construction count')
    require(len(set(row_keys(weighted))) == 2564, 'Weighted points are distinct')
    metric = np.diag([1]*14 + [2]).astype(np.int64)
    require(np.all(np.sum((weighted @ metric) * weighted, axis=1) == 16),
            'Weighted construction norms')

    transform = np.zeros((15, 15), dtype=np.int64)
    for j in range(7):
        transform[j, 2*j] = transform[j, 2*j+1] = 1
        transform[8+j, 2*j], transform[8+j, 2*j+1] = 1, -1
    transform[7, 14] = 2
    require(np.array_equal(transform.T @ transform, 2*metric),
            'Exact weighted-to-Euclidean isometry with scale squared 2')
    module = weighted @ transform.T
    require(np.all(np.sum(module * module, axis=1) == 32), 'Euclidean module norms')
    digest = canonical_hash(module)
    require(digest == recipe['module_sorted_int16_sha256'],
            'Pinned canonical construction identity')
    require(recipe['module_points'] == len(module), 'Recipe count')
    gram = module @ module.T
    np.fill_diagonal(gram, -32)
    require(int(gram.max()) <= 16, 'Every module pair by direct int64 products')

    tail_norm = np.sum(module[:, 8:] ** 2, axis=1)
    histogram = {int(k): int(v) for k, v in
                 zip(*np.unique(tail_norm, return_counts=True))}
    require(histogram == {0: 240, 16: 2240, 32: 84}, 'Module layer counts')
    equator = []
    for i, j in combinations(range(8), 2):
        for a, b in product((-4, 4), repeat=2):
            row = [0]*8
            row[i], row[j] = a, b
            equator.append(row)
    equator.extend(signs for signs in product((-2, 2), repeat=8)
                   if sum(sign < 0 for sign in signs) % 2 == 0)
    require(set(row_keys(module[tail_norm == 0, :8])) == set(row_keys(equator)),
            'The equator is the explicit 240-point E8 root configuration')
    tails = []
    for i, j in combinations(range(7), 2):
        for a, b in product((-4, 4), repeat=2):
            row = [0]*7
            row[i], row[j] = a, b
            tails.append(row)
    require(set(row_keys(module[tail_norm == 32, 8:])) == set(row_keys(tails)),
            'The polar layer is the explicit 84-point D7 root configuration')
    return module, digest, histogram, int(gram.max())


def replay(package=ROOT):
    started = time.monotonic()
    package = Path(package)
    data = package / 'data'
    provenance = json.loads((data / 'provenance.json').read_text())
    recipe = json.loads((data / 'construction.json').read_text())
    dimension = provenance['dimension']
    require(dimension == 31, 'This certificate embeds the full 15D module in 31D')
    full, module_hash, layers, module_maximum = reconstruct(recipe)
    with np.load(data / 'module.npz', allow_pickle=False) as z:
        r = integer_array(z['transformation_numerator'], 'Isometry numerator', 2**20)
        denominator = int(z['transformation_denominator'])
        oriented = integer_array(z['points'], 'Oriented module')
    require(r.shape == (15, 15) and 0 < denominator <= 2**20,
            'Exact rational-isometry format')
    require(np.array_equal(r @ r.T, denominator**2 * np.eye(15, dtype=np.int64)),
            'Exact rational isometry')
    numerator = full @ r.T
    require(not np.any(numerator % denominator), 'Integral transformed module')
    require(oriented.shape == (2564, 15)
            and set(row_keys(numerator // denominator)) == set(row_keys(oriented)),
            'Stored full module equals the exact transformed construction')

    name = provenance['coordinate_file']
    require(isinstance(name, str) and Path(name).name == name, 'Portable coordinate filename')
    with np.load(data / name, allow_pickle=False) as z:
        require(int(z['squared_norm']) == 128 and int(z['inner_product_limit']) == 64,
                'Full-coordinate normalization')
        require(int(z['coordinate_denominator']) == 2, 'Full-coordinate denominator')
        final128 = integer_array(z['points'], 'Full coordinates')
    require(final128.shape == (provenance['target_count'], dimension)
            and not np.any(final128 % 2), 'Full-coordinate shape and scaling')
    final = final128 // 2
    with np.load(data / 'exchange.npz', allow_pickle=False) as z:
        added128 = integer_array(z['added'], 'Added coordinates')
        removed128 = integer_array(z['removed'], 'Removed coordinates')
    require(added128.shape == (provenance['added_count'], dimension)
            and removed128.shape == (provenance['removed_count'], dimension)
            and not np.any(added128 % 2) and not np.any(removed128 % 2),
            'Minimal exchange shape and normalization')
    added, removed = added128 // 2, removed128 // 2
    require(np.all(np.sum(final * final, axis=1) == 32), 'All final norms')
    require(np.all(np.sum(removed * removed, axis=1) == 32), 'All removed norms')
    final_rows = row_keys(final)
    fk, ak, rk = set(final_rows), set(row_keys(added)), set(row_keys(removed))
    require(len(fk) == len(final) and len(ak) == len(added) and len(rk) == len(removed),
            'Distinct final points and distinct exchange points')
    require(ak <= fk and not (rk & fk) and not (ak & rk), 'Minimal exchange membership')
    source_keys = (fk - ak) | rk
    require(len(source_keys) == provenance['source_count'], 'Reconstructed source count')
    source = np.frombuffer(b''.join(source_keys), dtype='<i2').reshape(-1, dimension)
    source_hash = canonical_hash(source)
    require(source_hash == provenance['source_canonical_norm32_sha256'],
            'Exact pinned source identity after undoing the minimal exchange')
    require(len(final) - len(source) == len(added) - len(removed), 'Exchange count identity')

    normal = tuple(provenance.get('normal_octad', NORMAL))
    require(len(normal) == 8 and len(set(normal)) == 8
            and all(type(j) is int and 0 <= j < 24 for j in normal),
            'Normal-octad embedding')
    embedding = normal + tuple(range(24, 31))
    embedded = np.zeros((2564, dimension), dtype=np.int64)
    embedded[:, embedding] = oriented
    mk = set(row_keys(embedded))
    require(ak <= mk - source_keys, 'Every addition is a new transformed module point')
    whole_module = provenance.get('full_module_retained', False)
    require(type(whole_module) is bool, 'Whole-module retention flag')
    if whole_module:
        require(mk <= fk and ak == mk - source_keys,
                'The complete transformed module is retained with its exact new-point set')

    # All coefficients are at most 8 after dividing the stored coordinates by 2.
    # Thus each dot-product partial sum has magnitude at most 31*8^2. The
    # following loop uses signed integer multiplication and addition only.
    columns = np.ascontiguousarray(final.T)
    own_indices = {key: i for i, key in enumerate(final_rows) if key in ak}
    maximum = -32
    comparisons = 0
    for point, key in zip(added, row_keys(added)):
        dots = np.zeros(len(final), dtype=np.int64)
        for j in np.flatnonzero(point):
            dots += int(point[j]) * columns[j]
        own = own_indices[key]
        require(dots[own] == 32, 'Added-point self product')
        dots[own] = -32
        current = int(dots.max())
        require(current <= 16, 'Every changed pair by direct int64 coordinate sums')
        maximum = max(maximum, current)
        comparisons += len(final) - 1
    return dict(status='PASS', dimension=dimension, source_count=len(source),
                total=len(final), removed=len(removed), added=len(added),
                gain=len(final)-len(source), source_canonical_norm32_sha256=source_hash,
                reconstructed_module_count=len(full), reconstructed_layer_counts=layers,
                module_sorted_int16_sha256=module_hash,
                module_unordered_pairs_checked=len(full)*(len(full)-1)//2,
                maximum_module_product=module_maximum,
                retained_module_points=len(mk & fk), full_module_retained=mk <= fk,
                available_new_module_points=len(mk - source_keys),
                omitted_new_module_points=len((mk - source_keys) - fk),
                changed_point_comparisons=comparisons, maximum_changed_product=maximum,
                arithmetic='Direct int64 products and sums; no tolerance or floating matrix products',
                scope='Every changed pair; unchanged pairs are covered by the separate complete geometry replay',
                seconds=time.monotonic()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', type=Path, default=ROOT)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    result = replay(args.package)
    if args.report:
        args.report.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
