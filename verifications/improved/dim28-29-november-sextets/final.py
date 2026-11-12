"""Derive registered counts from the supplied complete coordinate arrays."""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parent


def results():
    out = []
    for record in json.loads((ROOT/'data/provenance.json').read_text())['records']:
        with np.load(ROOT/record['certificate'], allow_pickle=False) as z:
            n, d = z['points'].shape
        if (n, d) != (record['total'], record['dimension']):
            raise ValueError('Coordinate count or dimension')
        if n <= record['packingstar']:
            raise ValueError('Pinned comparison not exceeded')
        out.append(dict(dimension=d, total=n))
    return out


if __name__ == '__main__':
    for record in results():
        print(json.dumps(record))
