"""Search compatible unions with exact sextet isometries of the Leech coordinates.

SciPy is needed for this optional discovery script, but not for certificate replay.
Every accepted output receives a complete independent geometry check. A blocker
cutoff restricts the candidate pool; it does not certify an unrestricted optimum.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import argparse
import itertools
import json
from pathlib import Path
import tempfile
import numpy as np
from scipy.sparse import csr_matrix, vstack
from scipy.sparse.csgraph import maximum_bipartite_matching
from verify_geometry import require, verify

OCTAD = (0, 2, 4, 5, 6, 10, 11, 23)


def frames():
    generator = np.array([1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1], dtype=np.uint8)
    matrix = np.zeros((12, 23), dtype=np.uint8)
    for i in range(12):
        matrix[i, i:i+12] = generator
    bits = ((np.arange(4096)[:, None] >> np.arange(12)) & 1).astype(np.uint8)
    code = np.zeros((4096, 24), dtype=np.uint8)
    code[:, :23] = (bits @ matrix) % 2
    code[:, 23] = code[:, :23].sum(axis=1) % 2
    octads = sorted(sum(1 << int(j) for j in np.flatnonzero(word))
                    for word in code[code.sum(axis=1) == 8])
    require(len(octads) == 759, 'Golay octads')
    seen = set()
    result = []
    for tetrad in itertools.combinations(OCTAD, 4):
        mask = sum(1 << j for j in tetrad)
        blocks = [list(tetrad)] + [[j for j in range(24) if (o ^ mask) >> j & 1]
                                   for o in octads if o & mask == mask]
        require(len(blocks) == 6 and all(len(b) == 4 for b in blocks)
                and sorted(sum(blocks, [])) == list(range(24)), 'Sextet partition')
        key = tuple(sorted(tuple(b) for b in blocks))
        if key not in seen:
            seen.add(key)
            result.append(blocks)
    require(len(result) == 35, 'Sextets compatible with the normal octad')
    return result


def transform(x, blocks):
    # On each tetrad S = (J-2I)/2. Since J^2=4J, S^T S=I.
    # Negating the first tetrad supplies the Leech-lattice sign correction.
    h = np.ones((4, 4), dtype=np.int64)-2*np.eye(4, dtype=np.int64)
    require(np.array_equal(h @ h.T, 4*np.eye(4, dtype=np.int64)), 'Isometry')
    if any(np.any(x[:, b].sum(axis=1) % 2) for b in blocks):
        return None  # The integer-coordinate search does not round fractional images.
    y = x.copy()
    for b in blocks:
        y[:, b] = x[:, b].sum(axis=1)[:, None]//2-x[:, b]
    y[:, blocks[0]] *= -1
    require(np.all(np.sum(y.astype(np.int64)**2, axis=1) == 128), 'Image norms')
    return y


def refill(x, y, blocker_limit=None):
    original = {row.tobytes(): i for i, row in enumerate(x)}
    transformed = {row.tobytes(): i for i, row in enumerate(y)}
    common = original.keys() & transformed.keys()
    ai = np.array([i for k, i in original.items() if k not in common], dtype=np.int32)
    bi = np.array([i for k, i in transformed.items() if k not in common], dtype=np.int32)
    if not len(bi):
        return x, dict(gain=0, candidates=0, removed=0, added=0)
    a, b = x[ai], y[bi]
    require(x.shape[1]*16**2 < 2**24 and np.max(np.abs(y)) <= 16,
            'Exact integer BLAS bound')
    af = a.astype(np.float32)
    if blocker_limit is not None:
        ids = np.arange(len(b))
        bf = b.astype(np.float32)
        counts = np.zeros(len(b), dtype=np.int32)
        order = np.random.default_rng(281060).permutation(len(a))
        for at in range(0, len(a), 256):
            counts += (bf @ af[order[at:at+256]].T > 64).sum(axis=1, dtype=np.int32)
            good = counts <= blocker_limit
            ids, bf, counts = ids[good], bf[good], counts[good]
            if not len(ids):
                return x, dict(gain=0, candidates=0, removed=0, added=0)
        b = b[ids]
    bf = b.astype(np.float32)
    graph = vstack([csr_matrix(bf[at:at+256] @ af.T > 64)
                    for at in range(0, len(b), 256)], format='csr')
    matching = maximum_bipartite_matching(graph, perm_type='column')
    gain = len(b)-int(np.count_nonzero(matching >= 0))
    if gain == 0:
        return x, dict(gain=0, candidates=len(b), removed=0, added=0)
    reverse = np.full(len(a), -1, dtype=np.int64)
    for i, j in enumerate(matching):
        if j >= 0:
            reverse[j] = i
    left = set(map(int, np.flatnonzero(matching < 0)))
    right = set()
    queue = list(left)
    for i in queue:
        for j in graph.indices[graph.indptr[i]:graph.indptr[i+1]]:
            j = int(j)
            if matching[i] == j or j in right:
                continue
            right.add(j)
            if reverse[j] >= 0 and int(reverse[j]) not in left:
                left.add(int(reverse[j]))
                queue.append(int(reverse[j]))
    keep = np.ones(len(x), dtype=bool)
    keep[ai[sorted(right)]] = False
    result = np.unique(np.concatenate((x[keep], b[sorted(left)])), axis=0)
    require(len(result) == len(x)+gain, 'Matching augmentation count')
    return result, dict(gain=gain, candidates=len(b), removed=len(right), added=len(left))


def save(path, x):
    np.savez_compressed(path, points=x, squared_norm=128,
                        inner_product_limit=64, coordinate_denominator=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rounds', type=int, default=1)
    parser.add_argument('--frames', type=int, default=35)
    parser.add_argument('--max-blockers', type=int)
    args = parser.parse_args()
    require(args.rounds > 0 and 1 <= args.frames <= 35, 'Search limits')
    require(args.max_blockers is None or args.max_blockers >= 0, 'Blocker cutoff')
    with np.load(args.input, allow_pickle=False) as z:
        denominator = int(z['coordinate_denominator'])
        require(denominator in (1, 2) and int(z['squared_norm']) == 32*denominator**2,
                'Input normalization')
        raw = z['points']
        require(raw.dtype.kind in 'iu' and np.all(raw >= -16) and np.all(raw <= 16),
                'Input coefficients')
        x = raw.astype(np.int16)*(2//denominator)
    with tempfile.TemporaryDirectory() as directory:
        normalized = Path(directory)/'input.npz'
        save(normalized, x)
        verify(normalized)
    records = []
    for sweep in range(1, args.rounds+1):
        before = len(x)
        for index, blocks in enumerate(frames()[:args.frames], 1):
            y = transform(x, blocks)
            if y is None:
                print(json.dumps(dict(sweep=sweep, frame=index, skipped='fractional image')),
                      flush=True)
                continue
            x, record = refill(x, y, args.max_blockers)
            record.update(sweep=sweep, frame=index, total=len(x), blocks=blocks)
            records.append(record)
            print(json.dumps(record), flush=True)
        if len(x) == before:
            break
    save(args.output, x)
    proof = verify(args.output)
    args.output.with_suffix('.json').write_text(json.dumps(
        dict(search=records, verification=proof, blocker_limit=args.max_blockers), indent=2)+'\n')
    print(json.dumps(proof), flush=True)


if __name__ == '__main__':
    main()
