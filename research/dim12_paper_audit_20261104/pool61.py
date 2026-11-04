"""Test a new exact integer-ray pool for 61 points in the frozen-bridge region.

A 61-point result copied into both R6 factors, plus the 720 bridges, gives 842.
Failure of this finite MILP says nothing about unrestricted configurations.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import itertools,json,sys,time,warnings
from pathlib import Path
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack
from audit import FACTORS

ROOT=Path(__file__).resolve().parent


def main():
    start=time.monotonic();rows=[]
    for i in range(6):
        for s in [-1,1]:x=np.zeros(6,dtype=np.int64);x[i]=s;rows.append(x)
    for support in itertools.combinations(range(6),4):
        for signs in itertools.product([-1,1],repeat=4):
            x=np.zeros(6,dtype=np.int64);x[list(support)]=signs;rows.append(x)
    for i in range(6):
        for signs in itertools.product([-1,1],repeat=6):
            x=np.array(signs,dtype=np.int64);x[i]*=2;rows.append(x)
    for signs in itertools.product([-1,1],repeat=6):rows.append(np.array(signs,dtype=np.int64))
    raw=np.array(rows);n=len(raw);assert n==700
    norms=np.sum(raw*raw,axis=1);top=np.sort(abs(raw),axis=1)[:,-2:].sum(axis=1)
    assert np.all(top*top<=norms)
    dot=raw@raw.T;bad=(dot>0)&(4*dot*dot>norms[:,None]*norms[None,:]);np.fill_diagonal(bad,False)
    ii,jj=np.nonzero(np.triu(bad,1));m=len(ii)
    supports=[{0,1,2,3},{0,1,4,5},{2,3,4,5}]
    known=[i for i,row in enumerate(raw) if np.count_nonzero(row)==1
           or (norms[i]==4 and set(np.flatnonzero(row)) in supports)]
    assert len(known)==60 and not bad[np.ix_(known,known)].any()
    a=coo_matrix((np.ones(2*m),(np.r_[np.arange(m),np.arange(m)],np.r_[ii,jj])),shape=(m,n)).tocsc()
    a=vstack([a,coo_matrix(np.ones((1,n)))]).tocsc()
    lo=np.r_[np.full(m,-np.inf),61.];hi=np.r_[np.ones(m),np.inf]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected')
        opt=milp(np.zeros(n),integrality=np.ones(n),bounds=Bounds(0,1),
                 constraints=LinearConstraint(a,lo,hi),
                 options=dict(time_limit=45.,presolve=True,threads=1))
    result=dict(pool=n,conflicting_pairs=m,status=int(opt.status),message=str(opt.message),
                elapsed_seconds=time.monotonic()-start,found61=False,
                exact_known60_positive_control=True)
    if opt.x is not None:
        ids=np.flatnonzero(opt.x>.5)
        assert len(ids)>=61 and not bad[np.ix_(ids,ids)].any()
        head=raw[ids[:61]];bridges=[]
        for matching in FACTORS:
            for (i,j),(k,l) in itertools.product(matching,repeat=2):
                for signs in itertools.product([-1,1],repeat=4):
                    x=np.zeros(12);x[[i,j,k+6,l+6]]=signs;bridges.append(x)
        x=np.vstack([np.pad(head,((0,0),(0,6))),np.pad(head,((0,0),(6,0))),np.array(bridges)])
        np.save(ROOT/'integer842_candidate.npy',x)
        sys.path.insert(0,str(ROOT.parent/'dim12_quaternionic_20261104'))
        from verify import verify
        certificate=verify(ROOT/'integer842_candidate.npy');assert certificate['passed']
        result.update(found61=True,selected=ids[:61].tolist(),certificate=certificate)
    (ROOT/'pool61_report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
