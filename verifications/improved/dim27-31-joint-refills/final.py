"""Read each explicit coordinate count; verify.py establishes its geometry."""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parent


def results():
    records = json.loads((ROOT/'data/provenance.json').read_text())['records']
    out = []
    for record in records:
        with np.load(ROOT/record['certificate'], allow_pickle=False) as z:
            n, d = z['points'].shape
            assert d == record['dimension'] and n == record['total']
            assert n > record['packingstar']
        out.append(dict(dimension=d, total=n))
    return out


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
