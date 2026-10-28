"""Independently verify the exact optimum of a supplied finite point pool."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np


def leq(a,b,limit):
    q=limit-a
    return np.where(b>=0,(q>=0)&(2*b*b<=q*q),(q>=0)|(2*b*b>=q*q))


def verify(path):
    z=np.load(path,allow_pickle=False);a=z['rational'].astype(np.int64);b=z['radical'].astype(np.int64);col=z['clique'];ind=z['independent']
    assert int(z['radicand'])==2 and np.max(abs(a))<=8 and np.max(abs(b))<=8 and a.shape[1]<=15
    assert a.shape==b.shape and len(col)==len(a) and np.all((a*a+2*b*b).sum(1)==32) and np.all((a*b).sum(1)==0)
    assert len(np.unique(np.c_[a,b],axis=0))==len(a)
    assert len(np.unique(ind))==len(ind) and np.all(ind>=0) and np.all(ind<len(a))
    for k in np.unique(col):
        ids=np.flatnonzero(col==k);g=a[ids]@a[ids].T+2*b[ids]@b[ids].T;q=a[ids]@b[ids].T+b[ids]@a[ids].T
        assert not leq(g,q,16).any(),'A claimed conflict clique is not complete'
    g=a[ind]@a[ind].T+2*b[ind]@b[ind].T;q=a[ind]@b[ind].T+b[ind]@a[ind].T
    ok=leq(g,q,16);np.fill_diagonal(ok,True);assert np.all(ok),'Lower-bound witness is not independent'
    assert len(ind)==len(np.unique(col))
    return dict(status='PASS',pool_points=len(a),exact_finite_optimum=len(ind),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),scope='Only the supplied finite pool: a matching feasible subset and complete conflict-clique cover. Integer arithmetic; no optimizer needed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate',type=Path);p.add_argument('--report',type=Path);q=p.parse_args();r=verify(q.certificate)
    if q.report:q.report.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
