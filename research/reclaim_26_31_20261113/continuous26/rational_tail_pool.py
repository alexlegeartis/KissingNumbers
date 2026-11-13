"""Enumerate non-grid rational norm16 tails against the exact 26D exterior.

New points have integer E8 deep-hole heads (squared norm16) and tails
4(a,b)/D where a^2+b^2=D^2. The old sixteen Q(sqrt3) points are retained as
independent vertices, with exact sign-sensitive conflict comparisons.
"""
from itertools import combinations,product
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix


def heads():
    h=[]
    for i in range(8):
        for s in [-4,4]:
            r=np.zeros(8,np.int16);r[i]=s;h.append(r)
    for I in combinations(range(8),4):
        for s in product([-2,2],repeat=4):
            r=np.zeros(8,np.int16);r[list(I)]=s;h.append(r)
    for i in range(8):
        for s in product([-1,1],repeat=8):
            if sum(v<0 for v in s)%2!=1:continue
            r=np.asarray(s,np.int16);r[i]*=3;h.append(r)
    h=np.asarray(h);assert len(h)==2160 and np.all((h.astype(int)**2).sum(1)==16)
    return h


def greater3(p,q,limit):
    p=np.asarray(p,np.int64);q=np.asarray(q,np.int64);r=limit-p
    return np.where(q>=0,(r<0)|(3*q*q>r*r),(r<0)&(r*r>3*q*q))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--den',type=int,default=65)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=int,default=120)
    ap.add_argument('--include-axes',action='store_true')
    a=ap.parse_args();start=time.monotonic();z=np.load(Path(__file__).with_name('normal26.npz'))
    b=np.rint(z['B']).astype(np.int16);b=b[(b.astype(int)**2).sum(1)==16]
    H=heads();D=a.den;tails=[(x,y) for x in range(-D,D+1) for y in range(-D,D+1) if x*x+y*y==D*D and (a.include_axes or x*y!=0)]
    parts=[];rows=[]
    for t in tails:
        p=np.zeros((len(H),10),np.int32);p[:,:8]=D*H;p[:,8:]=4*np.asarray(t)
        good=(p.astype(np.float32)@b.astype(np.float32).T<=16*D).all(1)
        parts.append(p[good]);rows.append({'tail_numerators':list(t),'admissible':int(good.sum())})
    P=np.unique(np.vstack(parts),axis=0)
    old=np.rint(z['Y'][:364]).astype(np.int32); oldset={tuple(r) for r in D*old}
    P=np.asarray([p for p in P if tuple(p) not in oldset],np.int32).reshape((-1,10))
    N=np.vstack([D*old,P]);n=len(N)
    edges=np.argwhere(np.triu(N.astype(np.float64)@N.astype(np.float64).T>16*D*D,1))
    src=Path(__file__).resolve().parents[4]/'PackingStar/25D-31D/new_bounds_20261111/D26_201180_coordinates.npz'
    f=np.load(src);I=z['indices'];A=f['quadratic_A'][:,I].astype(int);B=f['quadratic_B'][:,I].astype(int)
    rr,cc=np.where(greater3(N@A.T,N@B.T,16*D));e2=np.c_[rr,n+cc]
    qr,qc=np.where(np.triu(greater3(A@A.T+3*(B@B.T),A@B.T+B@A.T,16),1));e3=np.c_[n+qr,n+qc]
    edges=np.vstack([edges,e2,e3]);total=n+len(A)
    M=coo_matrix((np.ones(2*len(edges)),(np.repeat(np.arange(len(edges)),2),edges.ravel())),shape=(len(edges),total)).tocsr()
    report={'denominator':D,'tail_count':len(tails),'tails':rows,'new_rational_points':len(P),'pool_points':total,'conflicts':len(edges),'baseline_sector':380,'scope':'Numerical MILP discovery; final changes require independent exact verification'}
    print(json.dumps(report),flush=True)
    np.savez_compressed(a.out.with_suffix('.npz'),rational=N,quadratic_A=A,quadratic_B=B,denominator=D,edges=edges,indices=I)
    res=milp(-np.ones(total),integrality=np.ones(total),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,1),options={'time_limit':a.seconds,'mip_rel_gap':0,'threads':2})
    report.update(status=int(res.status),message=res.message,seconds=time.monotonic()-start)
    if res.x is not None:
        ids=np.flatnonzero(res.x>0.5);assert not np.any(np.isin(edges[:,0],ids)&np.isin(edges[:,1],ids))
        report.update(selected_count=len(ids),selected_indices=ids.tolist(),mip_gap=float(res.mip_gap),mip_dual_bound=float(res.mip_dual_bound))
    a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['selected_indices','tails']}),flush=True)


if __name__=='__main__':main()
