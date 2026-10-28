"""Independent exact verifier for possibly asymmetric quadratic kissing codes.

Only NumPy is required. Existing opposite pairs are exploited for speed, but
antipodality is neither assumed nor required. Every unordered point pair is
covered by a signed or absolute exact quadratic comparison.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np


def leq(a,b,limit):
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64);q=limit-a
    return np.where(b>=0,(q>=0)&(2*b*b<=q*q),(q>=0)|(2*b*b>=q*q))


def check(a,b,c,e,absolute):
    g=a@c.T+2*b@e.T;q=a@e.T+b@c.T
    ok=leq(g,q,16)
    return ok&leq(-g,-q,16) if absolute else ok


def verify(path,block_size=128):
    start=time.time();raw=path.read_bytes()
    with np.load(path,allow_pickle=False) as z:
        a=z['rational'];b=z['radical']
        assert int(z['radicand'])==2 and int(z['squared_norm'])==32 and int(z['inner_product_limit'])==16
        assert a.shape==b.shape and a.ndim==2 and 26<=a.shape[1]<=31
        assert a.dtype.kind in 'iu' and b.dtype.kind in 'iu'
        assert np.all(a>=-8)&np.all(a<=8)&np.all(b>=-8)&np.all(b<=8)
    a=a.astype(np.int64);b=b.astype(np.int64);n,d=a.shape
    assert np.all((a*a+2*b*b).sum(1)==32) and np.all((a*b).sum(1)==0)
    coeff=np.c_[a,b];assert len(np.unique(coeff,axis=0))==n,'duplicate points'
    canonical=coeff.copy();sign=np.where(canonical[np.arange(n),np.argmax(canonical!=0,axis=1)]<0,-1,1);canonical*=sign[:,None]
    representatives,inv,counts=np.unique(canonical,axis=0,return_inverse=True,return_counts=True)
    assert np.all((counts==1)|(counts==2))
    paired=representatives[counts==2];unpaired=coeff[counts[inv]==1]
    assert 2*len(paired)+len(unpaired)==n
    integer=np.all(paired[:,d:]==0,axis=1);x=paired[integer,:d].astype(float);pa=paired[~integer,:d].astype(float);pb=paired[~integer,d:].astype(float)
    ua=unpaired[:,:d].astype(float);ub=unpaired[:,d:].astype(float);line_pairs=0
    for at in range(0,len(x),block_size):
        block=x[at:at+block_size];g=block@x[at:].T;i,j=np.tril_indices(len(block));g[i,j]=0
        assert np.max(abs(g))<=16,('paired integer inequality',at)
        line_pairs+=len(block)*(len(x)-at)-len(block)*(len(block)+1)//2
        if at%(64*block_size)==0:print('integer line rows',at,'/',len(x),'seconds',round(time.time()-start,1),flush=True)
    for at in range(0,len(x),4096):
        block=x[at:at+4096]
        assert np.all(leq(block@pa.T,block@pb.T,16)&leq(-block@pa.T,-block@pb.T,16)),('mixed paired inequality',at)
        line_pairs+=len(block)*len(pa)
    ok=check(pa,pb,pa,pb,True);np.fill_diagonal(ok,True);assert ok.all(),'quadratic paired inequality'
    line_pairs+=len(pa)*(len(pa)-1)//2;assert line_pairs==len(paired)*(len(paired)-1)//2
    for at in range(0,len(x),4096):
        block=x[at:at+4096];g=block@ua.T;q=block@ub.T
        assert np.all(leq(g,q,16)&leq(-g,-q,16)),('unpaired versus integer pair',at)
    assert check(ua,ub,pa,pb,True).all(),'unpaired versus quadratic pair'
    ok=check(ua,ub,ua,ub,False);np.fill_diagonal(ok,True);assert ok.all(),'unpaired inequality'
    cross=len(paired)*len(unpaired);singles=len(unpaired)*(len(unpaired)-1)//2
    covered=4*line_pairs+len(paired)+2*cross+singles;assert covered==n*(n-1)//2
    return dict(status='PASS',dimension=d,points=n,antipodal_pairs=len(paired),unpaired_points=len(unpaired),quadratic_paired_lines=len(pa),distinct_line_pairs_checked=line_pairs,paired_to_unpaired_checks=cross,unpaired_pairs_checked=singles,unordered_point_pairs_covered=covered,squared_norm=32,inner_product_limit=16,coefficient_sha256=hashlib.sha256(coeff.astype('<i2').tobytes()).hexdigest(),file_sha256=hashlib.sha256(raw).hexdigest(),seconds=time.time()-start,arithmetic='Bounded integer dot products exactly represented in float64; sign-sensitive integer comparisons in Q(sqrt(2)), without tolerance.',scope='All norms, distinctness and all unordered point-pair inequalities; no external baseline assumption and no antipodality requirement.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('coordinates',type=Path);p.add_argument('--report',type=Path);q=p.parse_args();result=verify(q.coordinates)
    if q.report:q.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
