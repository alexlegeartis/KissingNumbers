"""Derive the packaged count; verify.py independently proves every pair bound."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parent


def results():
    manifest=json.loads((ROOT/'data/provenance.json').read_text(encoding='utf-8'))
    path=ROOT/'data/coordinates31_251110.npz'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==manifest['certificate_sha256']
    with np.load(path,allow_pickle=False) as z:
        x=z['points'];assert x.dtype.kind in 'iu' and x.shape[1]==31
        assert np.all(x>=-8) and np.all(x<=8)
        assert np.all((x.astype(np.int32)**2).sum(1)==32)
        assert int(z['squared_norm'])==32 and int(z['inner_product_limit'])==16
        assert int(z['coordinate_denominator'])==1
        count=len(x);assert count==251110
    return [dict(dimension=31,total=count)]


if __name__=='__main__':
    for result in results():print(json.dumps(result))
