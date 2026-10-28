"""Read the certificate count; verify.py independently proves every pair bound."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent


def results():
    manifest=json.loads((ROOT/'data/provenance.json').read_text())
    path=ROOT/'data/coordinates27_209600.npz'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==manifest['certificate_sha256']
    with np.load(path,allow_pickle=False) as z:
        a=z['rational'];b=z['radical']
        assert a.shape==b.shape==(209600,27) and a.dtype.kind in 'iu' and b.dtype.kind in 'iu'
        assert np.max(abs(a))<=8 and np.max(abs(b))<=8
        a=a.astype(np.int32);b=b.astype(np.int32)
        assert np.all((a*a+2*b*b).sum(1)==32) and np.all((a*b).sum(1)==0)
        assert int(z['radicand'])==2 and int(z['squared_norm'])==32 and int(z['inner_product_limit'])==16
    return [dict(dimension=27,total=len(a))]


if __name__=='__main__':
    for row in results():print(json.dumps(row))
