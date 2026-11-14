"""Pinned source exchange, separate int64 delta checks, complete final geometry."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
import numpy as np
if not __debug__:raise RuntimeError('Run without -O: assertions are required.')
ROOT=Path(__file__).resolve().parent

def main():
    for name,digest in json.loads((ROOT/'manifest.json').read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    z=np.load(ROOT/'data/coordinates30_248212.npz',allow_pickle=False)
    A,B=z['rational'],z['radical'];assert A.shape==B.shape==(248212,30) and not np.any(B)
    s=np.load(ROOT/'data/source30_248180.npz',allow_pickle=False);O=s['points']
    assert O.shape==(248180,30) and int(s['squared_norm'])==128 and int(s['inner_product_limit'])==64
    e=np.load(ROOT/'data/exchange30.npz',allow_pickle=False);C,removed=e['added'],e['removed_indices']
    assert C.shape[1]==30 and len(C)-len(removed)==32
    assert len(np.unique(removed))==len(removed) and np.all((removed>=0)&(removed<len(O)))
    keep=np.ones(len(O),bool);keep[removed]=False
    keys=[v.tobytes() for v in A];assert len(set(keys))==len(A)
    assert {v.tobytes() for v in np.r_[O[keep],C]}==set(keys)
    donors=[np.load(ROOT/f'data/donor30_{n}.npz',allow_pickle=False)['points'] for n in (248020,247972)]
    donor_keys={v.tobytes() for D in donors for v in D}
    assert all(v.tobytes() in donor_keys for v in C)
    AI=A.astype(np.int64);index={v:i for i,v in enumerate(keys)}
    for v in C:
        dot=np.zeros(len(A),np.int64)
        for k in np.flatnonzero(v):dot+=int(v[k])*AI[:,k]
        dot[index[v.tobytes()]]=0
        assert dot.max()<=64
    print('EXACT_EXCHANGE',len(O),'+',len(C),'-',len(removed),'=',len(A),flush=True)
    with tempfile.TemporaryDirectory() as td:
        subprocess.run([sys.executable,str(ROOT/'verify_geometry.py'),str(ROOT/'data/coordinates30_248212.npz'),'--out',str(Path(td)/'verification.json')],check=True)
    print('PASS: 30D 248212; all final pairs exact, without a source-validity premise.',flush=True)
if __name__=='__main__':main()
