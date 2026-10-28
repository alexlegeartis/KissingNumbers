"""Derive certificate counts; verify.py proves all geometric inequalities."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent


def results():
    meta=json.loads((ROOT/'data/provenance.json').read_text());rows=[]
    for record in meta['records']:
        path=ROOT/record['certificate'];assert hashlib.sha256(path.read_bytes()).hexdigest()==record['certificate_sha256']
        with np.load(path,allow_pickle=False) as z:
            a=z['rational'];b=z['radical'];assert a.dtype.kind in 'iu' and b.dtype.kind in 'iu'
            assert a.shape==b.shape==(record['total'],record['dimension']) and not np.any(b)
            assert np.all(a>=-8) and np.all(a<=8)
            assert np.all((a.astype(np.int32)**2).sum(1)==32)
            assert int(z['radicand'])==2 and int(z['squared_norm'])==32 and int(z['inner_product_limit'])==16
        rows.append(dict(dimension=record['dimension'],total=len(a)))
    return rows


if __name__=='__main__':
    for row in results():print(json.dumps(row))
