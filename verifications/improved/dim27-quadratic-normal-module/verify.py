"""Standalone exact all-pairs verifier for x=A+B sqrt(r), r=2 or 3.

Only NumPy is needed. Integer-valued dot products use exact bounded float64
BLAS; irrational comparisons use integer sign cases and squaring.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np


def leq(a,b,limit,r):
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64);q=limit-a
    return np.where(b>=0,(q>=0)&(r*b*b<=q*q),(q>=0)|(r*b*b>=q*q))


def verify(path,block_size=128):
    start=time.time()
    with np.load(path,allow_pickle=False) as z:
        a=z['rational'];b=z['radical'];r=int(z['radicand'])
        assert int(z['squared_norm'])==32 and int(z['inner_product_limit'])==16
        assert r in (2,3) and a.shape==b.shape and a.ndim==2
        assert a.dtype.kind in 'iu' and b.dtype.kind in 'iu'
        assert 26<=a.shape[1]<=30 and np.max(abs(a))<=8 and np.max(abs(b))<=8
    a=a.astype(np.int64);b=b.astype(np.int64)
    assert np.all((a*a+r*b*b).sum(1)==32) and np.all((a*b).sum(1)==0)
    d=a.shape[1];allc=np.c_[a,b];u=np.unique(allc,axis=0)
    assert len(u)==len(allc),'duplicate points'
    positive=u[np.arange(len(u)),np.argmax(u!=0,axis=1)]>0
    half=u[positive];assert np.array_equal(half,np.unique(-u[~positive],axis=0))
    assert len(half)*2==len(u)
    integer=np.all(half[:,d:]==0,axis=1);x=half[integer,:d].astype(float)
    ra=half[~integer,:d].astype(float);rb=half[~integer,d:].astype(float)
    pairs=0;maximum=0
    for at in range(0,len(x),block_size):
        block=x[at:at+block_size];g=block@x[at:].T
        i,j=np.tril_indices(len(block));g[i,j]=0
        assert np.max(abs(g))<=16,('integer pair',at)
        maximum=max(maximum,int(np.max(abs(g))))
        pairs+=len(block)*(len(x)-at)-len(block)*(len(block)+1)//2
        if at%(64*block_size)==0:print('integer lines',at,'/',len(x),'seconds',round(time.time()-start,1),flush=True)
    for at in range(0,len(x),4096):
        aa=x[at:at+4096]@ra.T;bb=x[at:at+4096]@rb.T
        assert np.all(leq(aa,bb,16,r)&leq(-aa,-bb,16,r)),('mixed pair',at)
        pairs+=len(aa)*len(ra)
    aa=ra@ra.T+r*rb@rb.T;bb=ra@rb.T+rb@ra.T
    valid=leq(aa,bb,16,r)&leq(-aa,-bb,16,r);np.fill_diagonal(valid,True)
    assert np.all(valid),'radical pair';pairs+=len(ra)*(len(ra)-1)//2
    assert pairs==len(half)*(len(half)-1)//2
    return dict(status='PASS',dimension=d,points=len(a),radicand=r,squared_norm=32,inner_product_limit=16,antipodal_lines=len(half),irrational_lines=len(ra),distinct_line_pairs_checked=pairs,integer_max_absolute_inner_product=maximum,coefficient_sha256=hashlib.sha256(allc.astype('<i2').tobytes()).hexdigest(),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.time()-start,arithmetic='Exact bounded integer dot products; exact sign-sensitive quadratic comparisons, no tolerance.',scope='All norms, distinctness, antipodal partners, and every distinct line-pair inequality; no external baseline assumptions.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('coordinates',type=Path,nargs='?',default=Path(__file__).resolve().parent/'data/coordinates27_209600.npz');p.add_argument('--report',type=Path);p.add_argument('--block-size',type=int,default=128);args=p.parse_args()
    result=verify(args.coordinates,args.block_size)
    if args.report:args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
