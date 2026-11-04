"""Replace SDP cones by exact integer quadratic-form cuts, then certify an LP.

Each PSD cut uses an explicitly stored integer vector v and the integer
Fourier blocks: sum_j (v' B_j v) x_j >= 0. No eigenvalue tolerance enters
the validity of a cut. Floating point is used only to find useful vectors.
"""
import argparse,json,time
import numpy as np
from scipy.optimize import linprog
from labels_sdp import prepare,four_cliques
from model import DATA

def linear(z,four=True):
    nv=len(z['reps']);I=np.eye(nv,dtype=np.int64)
    E=I[[int(z['orbits'][0,0,0]),*map(int,np.flatnonzero(z['edges']))]]
    e=np.r_[1,np.zeros(len(E)-1)].astype(np.int64)
    A=np.r_[I,-I,-z['cliques'].astype(np.int64)]
    b=np.r_[np.zeros(nv),-np.ones(nv),-np.ones(len(z['cliques']))].astype(np.int64)
    if four:
        R=four_cliques(z).astype(np.int64);A=np.r_[A,-R];b=np.r_[b,-np.ones(len(R),dtype=np.int64)]
    return E,e,A,b

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--exclude',type=int,default=972);p.add_argument('--iterations',type=int,default=300);a=p.parse_args()
    start=time.monotonic();z=prepare();E,e,A,b=linear(z);blocks=np.rint(z['blocks']).astype(np.int64)
    cuts=[];seen=set();history=[]
    for it in range(a.iterations):
        norms=np.maximum(1,np.max(abs(A),axis=1))
        r=linprog(-z['row'].astype(float),A_ub=-A/norms[:,None],b_ub=-b/norms,
                  A_eq=E,b_eq=e,bounds=[(None,None)]*A.shape[1],method='highs')
        assert r.success,r.message
        bound=-r.fun;history.append(float(bound))
        if it%5==0:print(json.dumps(dict(event='bound',iteration=it,numerical_bound=bound,cuts=len(cuts),seconds=time.monotonic()-start)),flush=True)
        if bound<a.exclude-1e-5:break
        added=0
        for bi,B in enumerate(blocks):
            M=np.einsum('abi,i->ab',B,r.x);ev,Q=np.linalg.eigh(M)
            for j in np.flatnonzero(ev<-1e-7):
                v=np.rint(Q[:,j]*4096).astype(np.int64)
                c=np.einsum('a,abi,b->i',v,B,v)
                divisor=int(np.gcd.reduce(c));divisor=abs(divisor) or 1
                c=c//divisor;key=tuple(map(int,c))
                if key in seen:continue
                seen.add(key);A=np.r_[A,c[None,:]];b=np.r_[b,0]
                cuts.append(dict(block=bi,vector=list(map(int,v)),divisor=divisor));added+=1
        if not added:break
    report=dict(exclude=a.exclude,numerical_upper_bound=bound,cuts=cuts,history=history,
                iterations=len(history),seconds=time.monotonic()-start,certificate_checked=False)
    (DATA/'outer_cuts.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(DATA/'outer_lp.npz',E=E,e=e,A=A,b=b,row=z['row'])
    print(json.dumps({k:v for k,v in report.items() if k not in ('cuts','history')}),flush=True)
