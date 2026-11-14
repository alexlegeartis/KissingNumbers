"""Independent full-pair checker for integral A+B*sqrt(3) kissing codes.

Rational pairs use CPU float32 matrix multiplication exactly: every coordinate
is integral, and Cauchy bounds the sum of absolute products by 128, far below
the 2^24 exact-integer limit. Radical pairs use integer rational intervals for
sqrt(3), with the interval obtained by an integer square root. No discovery
or optimization code is imported.
"""
import argparse,hashlib,json,math,time
from pathlib import Path
import numpy as np
if not __debug__:raise RuntimeError('Run without -O: assertions are required.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--block',type=int,default=384)
    a=ap.parse_args();start=time.monotonic();z=np.load(a.candidate)
    assert np.issubdtype(z['rational'].dtype,np.integer) and np.issubdtype(z['radical'].dtype,np.integer)
    assert int(z['coordinate_denominator'])==2
    A=z['rational'].astype(np.int64);B=z['radical'].astype(np.int64)
    assert A.ndim==2 and B.shape==A.shape and int(z['radicand'])==3
    assert int(z['squared_norm'])==128 and int(z['inner_product_limit'])==64
    assert np.all(np.sum(A*A+3*B*B,axis=1)==128)
    assert not np.any(np.sum(A*B,axis=1))
    rational=np.flatnonzero(~np.any(B,axis=1));radical=np.flatnonzero(np.any(B,axis=1))
    F=A[rational].astype(np.float32);rational_max=-128
    for at in range(0,len(F),a.block):
        end=min(len(F),at+a.block);D=F[at:end]@F[at:].T
        D[np.arange(end-at),np.arange(end-at)]=-128
        maximum=int(D.max()) if D.size else -128;rational_max=max(rational_max,maximum)
        if maximum>64:
            i,j=np.argwhere(D>64)[0];raise AssertionError(('rational pair',int(rational[at+i]),int(rational[at+j]),int(D[i,j])))
        if at%(a.block*24)==0:print('RATIONAL_PAIRS',at,len(F),'max',rational_max,'seconds',round(time.monotonic()-start,2),flush=True)
    scale=1<<40;lo=math.isqrt(3*scale*scale);assert lo*lo<3*scale*scale<(lo+1)*(lo+1)
    for j in radical:
        p=A@A[j]+3*(B@B[j]);q=A@B[j]+B@A[j]
        lower=p*scale+q*np.where(q>=0,lo,lo+1)
        upper=p*scale+q*np.where(q>=0,lo+1,lo)
        lower[j]=upper[j]=0
        if np.any(lower>64*scale):raise AssertionError(('radical pair',int(j),np.flatnonzero(lower>64*scale)[:8].tolist()))
        assert not np.any(upper>64*scale),('interval inconclusive',int(j))
    report=dict(verified=True,all_pairs_checked=True,dimension=A.shape[1],count=len(A),
                unordered_pairs=len(A)*(len(A)-1)//2,rational_points=len(rational),radical_points=len(radical),
                squared_norm=128,inner_product_limit=64,rational_pair_maximum=rational_max,
                coordinate_denominator=int(z['coordinate_denominator']),
                sha256=hashlib.sha256(a.candidate.read_bytes()).hexdigest(),seconds=time.monotonic()-start)
    a.out.write_text(json.dumps(report,indent=2)+'\n');print('VERIFIED_FULL_CODE',json.dumps(report),flush=True)

if __name__=='__main__':main()
