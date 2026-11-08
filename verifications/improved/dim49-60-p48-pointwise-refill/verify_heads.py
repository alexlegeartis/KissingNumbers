"""Read-only local-frame int64 proof of the antipodal pointwise-refill heads.

Independent of discovery ownership maps and cached conflict lists. The second
checker independently materializes global coordinates and uses exactly bounded
binary64 arithmetic instead of these int64 local-frame products.
"""
from pathlib import Path
import hashlib
import io
import json
import numpy as np

GRAM_SHA256 = '01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer_array(a):
    require(a.dtype.kind in 'iu', 'Certificate arrays must be integers')
    require(not a.size or (int(a.min()) >= -2**31 and int(a.max()) < 2**31),
            'Certificate integer exceeds the checked input range')
    return a.astype(np.int64)


def magnitude(a):
    return max(abs(int(a.min())), abs(int(a.max()))) if a.size else 0


def product(a, b):
    require(a.ndim == b.ndim == 2 and a.shape[1] == b.shape[0], 'Product shape')
    require(a.shape[1] * magnitude(a) * magnitude(b) < 2**63, 'int64 product bound')
    return a @ b


def norms(x, gram):
    xg = product(x, gram)
    require(x.shape[1] * magnitude(xg) * magnitude(x) < 2**63, 'int64 norm bound')
    require(np.all(np.sum(xg * x, axis=1) == 6), 'A head is not a norm-six lattice vector')


def canonical(x):
    require(np.all(np.any(x != 0, axis=1)), 'Zero head')
    j = np.argmax(x != 0, axis=1)
    return x * np.sign(x[np.arange(len(x)), j])[:, None]


def source(source_dir):
    source_dir = Path(source_dir)
    gram_path = source_dir / 'p48p_gram.npy'
    gram_bytes = gram_path.read_bytes()
    family_bytes = (source_dir / 'family_refilled.npz').read_bytes()
    require(hashlib.sha256(gram_bytes).hexdigest() == GRAM_SHA256,
            'The Gram file is not the pinned published P48 Gram')
    gram = integer_array(np.load(io.BytesIO(gram_bytes), allow_pickle=False))
    with np.load(io.BytesIO(family_bytes), allow_pickle=False) as data:
        base = integer_array(data['base'])
        require(np.array_equal(gram, data['gram']), 'Source Gram disagreement')
    require(gram.shape == (48, 48) and np.array_equal(gram, gram.T), 'Gram shape or symmetry')
    require(base.shape == (7783, 48), 'Base shape')
    norms(base, gram)
    require(len(np.unique(canonical(base), axis=0)) == len(base), 'Duplicate base line')
    bg = product(base, gram)
    count = 0
    for start in range(0, len(base), 256):
        h = product(bg[start:start + 256], base.T)
        for offset, row in enumerate(h):
            i = start + offset
            require(np.all(np.abs(row[:i]) <= 2) and np.all(np.abs(row[i + 1:]) <= 2),
                    'Base pair exceeds absolute raw product two')
            count += i
    return base, gram, count


def certificate(path, base, gram, expected_classes=None):
    path = Path(path)
    raw = path.read_bytes()
    with np.load(io.BytesIO(raw), allow_pickle=False) as data:
        transforms = integer_array(data['transforms'])
        bits = data['retained_bits'].copy()
        extra = integer_array(data['extra_points'])
        offsets = integer_array(data['extra_offsets'])
        if 'base' in data:
            require(np.array_equal(data['base'], base), 'Certificate base disagreement')
        if 'gram' in data:
            require(np.array_equal(data['gram'], gram), 'Certificate Gram disagreement')
    m = len(transforms)
    require(m > 0 and transforms.shape == (m, 48, 48), 'Transform shape')
    require(expected_classes is None or m == expected_classes, 'Wrong number of tail classes')
    require(bits.dtype == np.uint8 and bits.shape == (m, (len(base) + 7) // 8), 'Mask shape')
    bits = np.unpackbits(bits, axis=1)
    require(not bits[:, len(base):].any(), 'Nonzero mask padding')
    bits = bits[:, :len(base)].astype(bool)
    require(extra.ndim == 2 and extra.shape[1] == 48, 'Extra head shape')
    require(offsets.shape == (m + 1,) and offsets[0] == 0 and offsets[-1] == len(extra)
            and np.all(np.diff(offsets) >= 0), 'Extra partition')
    norms(extra, gram)
    keys, checked = [], 0
    sizes = []
    for i, transform in enumerate(transforms):
        require(np.array_equal(product(product(transform, gram), transform.T), gram),
                'Transform does not preserve the Gram')
        retained = base[bits[i]]
        added = extra[offsets[i]:offsets[i + 1]]
        ag = product(added, gram)
        for start in range(0, len(added), 256):
            h = product(ag[start:start + 256], retained.T)
            require(np.all(np.abs(h) <= 2), 'Added/retained class conflict')
        checked += len(added) * len(retained)
        for start in range(0, len(added), 256):
            h = product(ag[start:start + 256], added.T)
            for offset, row in enumerate(h):
                j = start + offset
                require(row[j] == 6 and np.all(np.abs(row[:j]) <= 2)
                        and np.all(np.abs(row[j + 1:]) <= 2), 'Added/added class conflict')
                checked += j
        x = canonical(product(np.concatenate((retained, added)), transform))
        require(magnitude(x) < 2**15, 'Full-coordinate serialization range')
        norms(x, gram)
        x = np.ascontiguousarray(x, dtype=np.int16)
        keys.append(x.view(np.dtype((np.void, 96))).ravel())
        sizes.append(len(x))
    allkeys = np.concatenate(keys)
    require(len(np.unique(allkeys)) == len(allkeys), 'Duplicate global mother-line ownership')
    return dict(status='PASS', certificate_sha256=hashlib.sha256(raw).hexdigest(),
                classes=m, selected_lines=len(allkeys), extra_lines=len(extra),
                extra_pairs=checked, class_sizes=sizes,
                arithmetic='int64 with checked absolute product and partial-sum bounds',
                ownership='exact complete canonical 48-coordinate int16 records')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate')
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parent / 'source')
    args = parser.parse_args()
    base, gram, count = source(args.source)
    report = certificate(args.certificate, base, gram)
    report['base_pairs'] = count
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
