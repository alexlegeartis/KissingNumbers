"""Derive the certificate count; run verify.py for the complete geometric proof."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def results():
    with np.load(ROOT/'source/core.npz',allow_pickle=False) as c, np.load(ROOT/'source/auxiliary.npz',allow_pickle=False) as a:
        arrays=[c['coefficients'],a['norm36'],a['norm30']]
        assert [x.shape for x in arrays]==[(1084590,40),(40590,40),(1134,40)]
        assert all(x.dtype==np.int8 for x in arrays)
        total=sum(len(x) for x in arrays)
        assert total==1126314
    return [dict(dimension=40,total=total)]

if __name__=='__main__':
    for row in results(): print(json.dumps(row))
