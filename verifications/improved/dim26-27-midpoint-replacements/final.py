"""Derive both midpoint-replacement counts from their coordinate certificates."""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parent


def results():
    rows = []
    for dimension, total, removed, inserted, pure in (
            (26, 200918, 4352, 8704, 6), (27, 204636, 8064, 16128, 12)):
        record = json.loads((ROOT / ('source/record%d.json' % dimension)).read_text())
        with np.load(ROOT / 'source' / record['filename'], allow_pickle=False) as data:
            new = data['new_points']
            old = data['removed_leech']
            assert new.shape == (inserted, dimension) and old.shape == (removed, 24)
            assert np.issubdtype(new.dtype, np.integer) and np.issubdtype(old.dtype, np.integer)
            assert np.abs(new.astype(np.int64)).max() <= 4
            assert np.abs(old.astype(np.int64)).max() <= 4
            assert data['midpoint_owner_indices'].shape == (inserted, 2)
        assert record['dimension'] == dimension and record['points'] == total
        assert total == 196560 - removed + inserted + pure
        rows.append(dict(dimension=dimension, total=total))
    return rows


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
