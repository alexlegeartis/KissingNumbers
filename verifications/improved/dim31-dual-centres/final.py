"""Derive the certificate count; run verify.py for the complete geometric proof."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def results():
    with np.load(ROOT/'source/certificate31_241742.npz',allow_pickle=False) as z:
        for k in z.files:
            a=z[k]; assert np.issubdtype(a.dtype,np.integer)
            assert np.all(a>=-65535) and np.all(a<=65535)
        for k in ('midpoints','centres','removed_leech','centre_owner_pairs','pure_tails'):
            assert np.max(np.abs(z[k].astype(np.int64)))<=4
        assert z['midpoints'].shape==(86528,31) and z['centres'].shape==(1792,31)
        assert z['removed_leech'].shape==(43264,24) and z['pure_tails'].shape==(126,7)
        assert np.array_equal(z['centre_scale_squared'],[8,3])
        total=196560-len(z['removed_leech'])+len(z['midpoints'])+len(z['centres'])+len(z['pure_tails'])
        assert total==241742
    return [dict(dimension=31,total=total)]

if __name__=='__main__':
    for row in results(): print(json.dumps(row))
