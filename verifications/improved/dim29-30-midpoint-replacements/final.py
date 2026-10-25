"""Derive counts and validate coordinate ranges before the exact replay."""
from pathlib import Path
import json
import numpy as np
ROOT = Path(__file__).resolve().parent

def results():
    rows = []
    for d, removed, inserted, pure in ((29,19968,39936,40),(30,27648,55296,72)):
        total = 196560-removed+inserted+pure
        with np.load(ROOT/'source'/('certificate%d_%d.npz'%(d,total)), allow_pickle=False) as data:
            for key in ('new_points','removed_leech','midpoint_owner_indices','metric_weights'):
                a = data[key]
                assert np.issubdtype(a.dtype, np.integer)
                assert np.all(a >= -32768) and np.all(a <= 32767)
            assert data['new_points'].shape == (inserted,d)
            assert data['removed_leech'].shape == (removed,24)
            assert data['midpoint_owner_indices'].shape == (inserted,2)
            assert data['metric_weights'].shape == (d,)
            assert np.all(data['metric_weights'] == 1)
            assert np.max(np.abs(data['new_points'].astype(np.int64))) <= 4
            assert np.max(np.abs(data['removed_leech'].astype(np.int64))) <= 4
        rows.append(dict(dimension=d,total=total))
    return rows

if __name__ == '__main__':
    for row in results(): print(json.dumps(row))
