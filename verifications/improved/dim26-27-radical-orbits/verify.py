"""Pinned exchange reconstruction, separate integer delta proof, full geometry."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
import numpy as np
if not __debug__:raise RuntimeError('Run without -O: assertions are required.')
ROOT=Path(__file__).resolve().parent
def main():
    manifest=json.loads((ROOT/'manifest.json').read_text())
    for name,digest in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    s=np.load(ROOT/'data/source26_201180.npz');z=np.load(ROOT/'data/coordinates26_201182.npz')
    e=np.load(ROOT/'data/exchange26.npz');A=z['rational'];B=z['radical']
    O=2*np.r_[s['rational_points'],s['quadratic_A']]
    OB=2*np.r_[np.zeros_like(s['rational_points']),s['quadratic_B']]
    assert A.shape==B.shape==(201182,26)
    addA,addB,removed=e['added_A'],e['added_B'],e['removed_indices']
    assert addA.shape==addB.shape==(68,26) and removed.shape==(66,)
    assert len(np.unique(removed))==66 and np.all((removed>=0)&(removed<len(O)))
    keep=np.ones(len(O),bool);keep[removed]=False
    rebuilt={np.r_[v,b].tobytes() for v,b in zip(np.r_[O[keep],addA],np.r_[OB[keep],addB])}
    keys=[np.r_[v,b].tobytes() for v,b in zip(A,B)]
    assert len(set(keys))==len(A) and rebuilt==set(keys)
    # Separate coordinate-at-a-time integer dot products and algebraic signs.
    AI,BI=A.astype(np.int64),B.astype(np.int64)
    keyindex={k:i for i,k in enumerate(keys)}
    for ca,cb in zip(addA,addB):
        p=np.zeros(len(A),np.int64);q=p.copy()
        for k in np.flatnonzero(ca):
            p+=int(ca[k])*AI[:,k];q+=int(ca[k])*BI[:,k]
        for k in np.flatnonzero(cb):
            p+=3*int(cb[k])*BI[:,k];q+=int(cb[k])*AI[:,k]
        t=64-p
        bad=(q==0)&(p>64)
        positive=q>0;bad|=positive&((t<0)|(3*q*q>t*t))
        negative=q<0;bad|=negative&(t<0)&(t*t>3*q*q)
        bad[keyindex[np.r_[ca,cb].tobytes()]]=False
        assert not bad.any()
    print('SOURCE_AND_EXCHANGE_EXACT',len(O),'+',len(addA),'-',len(removed),'=',len(A),flush=True)
    with tempfile.TemporaryDirectory() as td:
        subprocess.run([sys.executable,str(ROOT/'verify_geometry.py'),str(ROOT/'data/coordinates26_201182.npz'),
                        '--out',str(Path(td)/'verification.json')],check=True)
    print('PASS: 26D 201182, full exact final geometry; no source-validity premise.',flush=True)
if __name__=='__main__':main()
