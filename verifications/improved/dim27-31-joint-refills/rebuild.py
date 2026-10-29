"""Rebuild a certificate from a pinned public input and the explicit delta."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from verify_geometry import canonical_hash, require

ROOT = Path(__file__).resolve().parent


def rebuild(dimension, source, output):
    records = json.loads((ROOT/'data/provenance.json').read_text())['records']
    record = next(r for r in records if r['dimension'] == dimension)
    require(hashlib.sha256(source.read_bytes()).hexdigest() == record['source_sha256'], 'Pinned input SHA-256')
    with np.load(source, allow_pickle=False) as z:
        denominator = int(z['coordinate_denominator'])
        require(denominator in (1, 2), 'Input scale')
        x = z['points'].astype(np.int16)*(2//denominator)
    require(canonical_hash(np.unique(x, axis=0)) == record['source_canonical_sha256'], 'Input coordinate hash')
    with np.load(ROOT/record['delta'], allow_pickle=False) as z:
        removed, added = z['removed'], z['added']
    points = {row.tobytes(): row for row in x}
    for row in removed:
        require(row.tobytes() in points, 'Removed-point membership')
        del points[row.tobytes()]
    for row in added:
        require(row.tobytes() not in points, 'Insertion distinctness')
        points[row.tobytes()] = row
    x = np.unique(np.array(list(points.values()), dtype=np.int16), axis=0)
    require(len(x) == record['total'] and canonical_hash(x) == record['canonical_point_sha256'], 'Rebuilt coordinates')
    np.savez_compressed(output, points=x, squared_norm=128, inner_product_limit=64, coordinate_denominator=2)
    print(json.dumps(dict(status='PASS', dimension=dimension, points=len(x), output=str(output))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dimension', type=int)
    parser.add_argument('source', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rebuild(args.dimension, args.source, args.output)
