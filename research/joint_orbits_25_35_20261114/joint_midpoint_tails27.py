"""Joint midpoint/root/rational-head exchanges with shared source deletions.

Four sign orbits are optimized simultaneously; the current head frames are
not fixed. Full profile multiplicities price every incumbent deletion.
A positive solution is a proposal pending independent full geometry replay.
"""
import argparse,hashlib,itertools,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reclaim_26_31_20261113'))
from search_rational_sphere_tails27 import heads,tails
sys.path.append('/workspace/work/deps')
from ortools.sat.python import cp_model
if not __debug__:raise RuntimeError('Run without -O.')
ROOT=Path(__file__).resolve().parents[2]
N=np.array([0,2,4,5,6,10,11,23]);I=np.r_[N,[24,25,26]]
def canon(v):
    v=v.copy()
    for sl in [slice(0,8),slice(8,11)]:
        a=v[sl];nz=np.flatnonzero(a)
        if len(nz) and a[nz[0]]<0:v[sl]*=-1
    return v

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--cost',type=int,default=16);ap.add_argument('--seconds',type=float,default=180);ap.add_argument('--seed',type=int,default=701);ap.add_argument('--K',type=int,default=2)
    ap.add_argument('--body-cost',type=int,default=0);ap.add_argument('--require-fractional',action='store_true');ap.add_argument('--local',action='store_true');a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic()
    src=ROOT/'verifications/improved/dim26-27-radical-orbits/data/coordinates27_209884.npz';z=np.load(src);assert not np.any(z['radical']) and not np.any(z['rational']%2)
    O=(z['rational']//2).astype(np.int16);proj=O[:,I];canonical=np.array([canon(v) for v in proj]);P,inv,w=np.unique(canonical,axis=0,return_inverse=True,return_counts=True);oldkeys={v.tobytes() for v in np.unique(proj,axis=0)}
    cat={}
    def add(v,d=1,kind='head'):
        v=canon(np.asarray(v,np.int64));g=int(np.gcd.reduce(np.r_[np.abs(v),d]));v//=g;d//=g
        assert int(v@v)==32*d*d
        if d==1 and v.astype(np.int16).tobytes() in oldkeys:return
        key=(v.tobytes(),d);cat.setdefault(key,(v,d,kind))
    modeldata=ROOT/'verifications/improved/dim26-27-radical-orbits/data/model_orbits27.npz';m=np.load(modeldata)
    for h,c in zip(m['H'],m['colors']):
        assert not np.any(h%2);add(np.r_[h[N]//2,m['tail_palette'][c]//2],kind='midpoint')
    L=np.load(ROOT/'research/26_27_continuation_20261114/data/leech.npy');J=np.setdiff1d(np.arange(24),N)
    for v in L[~np.any(L[:,J],axis=1)]:add(np.r_[v[N],np.zeros(3,int)],kind='root')
    H=heads();T,D=tails(a.K)
    for t,d in zip(T,D):
        for h in H:add(np.r_[h*d,t],int(d))
    vals=[];blocks=[];prices=[]
    PC,PT=P[:,:8].astype(np.float32),P[:,8:].astype(np.float32)
    allvals=list(cat.values())
    for at in range(0,len(allvals),96):
        part=allvals[at:at+96];V=np.array([v[0] for v in part],np.float32);den=np.array([v[1] for v in part]);bad=np.abs(V[:,:8]@PC.T)+np.abs(V[:,8:]@PT.T)>16*den[:,None];cost=bad@w
        limit=np.array([a.body_cost if v[2]!='head' and a.body_cost else a.cost for v in part])
        for i in np.flatnonzero(cost<=limit):
            vals.append(part[i]);blocks.append(np.flatnonzero(bad[i]));prices.append(int(cost[i]))
    A=np.array([v[0] for v in vals],np.int64);den=np.array([v[1] for v in vals],np.int64);weight=2+2*(np.any(A[:,:8],axis=1)&np.any(A[:,8:],axis=1));used=np.unique(np.concatenate(blocks));where={int(v):i for i,v in enumerate(used)}
    print('PRICED',json.dumps(dict(raw=len(cat),admitted=len(vals),profiles=len(P),kinds={k:sum(v[2]==k for v in vals) for k in ['root','midpoint','head']},cost=a.cost,seconds=time.monotonic()-start)),flush=True)
    model=cp_model.CpModel();x=[model.new_bool_var(f'x{i}') for i in range(len(A))];y=[model.new_bool_var(f'y{i}') for i in used]
    for i,bs in enumerate(blocks):
        for b in bs:model.add(x[i]<=y[where[int(b)]])
    edges=0;AF=A.astype(np.float32)
    for at in range(0,len(A),96):
        p=np.abs(AF[at:at+96,:8]@AF[:,:8].T)+np.abs(AF[at:at+96,8:]@AF[:,8:].T)
        bad=p>16*(den[at:at+96,None]*den[None,:])
        for i,j in zip(*np.nonzero(bad)):
            i+=at
            if i<j:model.add_at_most_one(x[int(i)],x[int(j)]);edges+=1
    gain=sum(int(c)*v for c,v in zip(weight,x))-sum(int(w[i])*v for i,v in zip(used,y));model.add(gain>=1)
    if a.require_fractional:model.add(sum(x[i] for i in np.flatnonzero(den>1))>=1)
    if a.local:model.minimize(sum(int(w[i])*v for i,v in zip(used,y)))
    else:model.maximize(gain)
    for v in x+y:model.add_hint(v,False)
    print('MODEL',len(A),'edges',edges,'seconds',time.monotonic()-start,flush=True)
    solver=cp_model.CpSolver();solver.parameters.num_search_workers=1;solver.parameters.max_time_in_seconds=a.seconds;solver.parameters.random_seed=a.seed;solver.parameters.log_search_progress=True
    st=solver.solve(model);report=dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),status=solver.status_name(st),candidates=len(A),edges=edges,cost=a.cost,body_cost=a.body_cost,require_fractional=a.require_fractional,local=a.local,K=a.K,seed=a.seed,seconds=time.monotonic()-start)
    if st in [cp_model.FEASIBLE,cp_model.OPTIMAL]:
        sel=np.flatnonzero([solver.boolean_value(v) for v in x]);dp=used[np.flatnonzero([solver.boolean_value(v) for v in y])];removed=np.flatnonzero(np.isin(inv,dp));common=int(np.lcm.reduce(den[sel]));added=[]
        for j in sel:
            v=A[j]*(common//den[j]);cs=(-1,1) if np.any(v[:8]) else (1,);ts=(-1,1) if np.any(v[8:]) else (1,)
            for s,t in itertools.product(cs,ts):
                point=np.zeros(27,np.int64);point[I]=np.r_[s*v[:8],t*v[8:]];added.append(point)
        added=np.array(added,np.int64);keep=np.ones(len(O),bool);keep[removed]=False;RA=np.r_[O[keep].astype(np.int64)*common,added]
        assert len(RA)>len(O) and len(RA)-len(O)==solver.value(gain)
        assert len(np.unique(RA,axis=0))==len(RA) and np.all(np.sum(RA*RA,axis=1)==32*common*common)
        for i,point in enumerate(added):
            dots=RA@point;dots[len(O)-len(removed)+i]=0;assert dots.max()<=16*common*common
        np.savez_compressed(a.out/'candidate.npz',points=RA,coordinate_denominator=common,squared_norm=32*common*common,inner_product_limit=16*common*common)
        report.update(count=len(RA),gain=len(RA)-len(O),added=len(added),removed=len(removed),denominator=common)
        print('POSITIVE_DELTA_CHECKED',json.dumps(report),flush=True)
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
