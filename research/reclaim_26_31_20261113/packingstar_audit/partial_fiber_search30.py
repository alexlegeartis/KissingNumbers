"""Shared-deletion MILPs on coupled tail neighborhoods of the 30D module.

All admitted points have four incumbent blockers. Neighborhoods comprise
four, six, or eight incumbent tail fibers; every point whose complete blocker
set lies in the neighborhood is considered up to the explicit pool limit.
The exact orthoplex capacity16 for each fixed tail strengthens the model.
Only independently checked positive witnesses are improvements. Numerical
solver bounds and interrupted searches are not exact no-gain certificates.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

NORMAL=np.array([0,2,4,5,6,10,11,23])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('source',type=Path)
    ap.add_argument('pool',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--models',type=int,default=18)
    ap.add_argument('--seconds',type=float,default=6)
    ap.add_argument('--limit',type=int,default=5000)
    ap.add_argument('--seed',type=int,default=20261113)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    z=np.load(args.source);X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:X//=2
    p=np.load(args.pool);C=p['candidates'].astype(np.int16);B=p['blockers'].reshape((-1,4))
    assert len(C)==len(B)
    victims=np.unique(B);oldtail,ti=np.unique(X[victims,24:],axis=0,return_inverse=True)
    assert len(victims)==960 and len(oldtail)==60
    assert np.all(np.bincount(ti)==16)
    outside=[j for j in range(24) if j not in NORMAL]
    assert not np.any(X[victims][:,outside])
    victim_index={int(v):i for i,v in enumerate(victims)}
    BT=np.array([[ti[victim_index[int(v)]] for v in b] for b in B],np.int16)
    patterns=defaultdict(list)
    for i,b in enumerate(BT):
        key=tuple(sorted(set(map(int,b))))
        if len(key)==4:patterns[key].append(i)
    ranked=sorted(patterns,key=lambda k:-len(patterns[k]))
    rng=np.random.default_rng(args.seed);neighborhoods=[]
    for key in ranked[:max(6,args.models//3)]:neighborhoods.append(key)
    for key in ranked:
        matches=[k for k in ranked if len(set(k)&set(key))==2]
        if matches:
            k=matches[int(rng.integers(len(matches)))];neighborhoods.append(tuple(sorted(set(k)|set(key))))
        if len(neighborhoods)>=args.models:break
    records=[];best=0;started=time.monotonic()
    for m,tailset in enumerate(neighborhoods[:args.models]):
        at=time.monotonic();idx=np.flatnonzero(np.all(np.isin(BT,tailset),axis=1))
        if len(idx)>args.limit:
            # Explicit random cap, preserving all four-tail bridge candidates first.
            bridge=idx[np.array([len(set(v))==4 for v in BT[idx]])]
            other=np.setdiff1d(idx,bridge)
            if len(bridge)>=args.limit:idx=np.sort(rng.choice(bridge,args.limit,replace=False))
            else:idx=np.sort(np.r_[bridge,rng.choice(other,args.limit-len(bridge),replace=False)])
        Y=C[idx];blocks=B[idx];old=np.unique(blocks);n=len(Y);oldmap={int(v):n+i for i,v in enumerate(old)}
        rows=[];cols=[];values=[];upper=[]
        for i,bs in enumerate(blocks):
            for v in bs:
                r=len(upper);rows.extend([r,r]);cols.extend([i,oldmap[int(v)]]);values.extend([1,-1]);upper.append(0)
        new_conflicts=0
        for start in range(0,n,256):
            D=Y[start:start+256].astype(np.float32)@Y.astype(np.float32).T
            for loc in range(len(D)):
                i=start+loc;js=np.flatnonzero(D[loc,i+1:]>16)+i+1
                new_conflicts+=len(js)
                for j in js:
                    r=len(upper);rows.extend([r,r]);cols.extend([i,int(j)]);values.extend([1,1]);upper.append(1)
        groups=defaultdict(list)
        for i,t in enumerate(Y[:,8:]):groups[tuple(map(int,t))].append(i)
        oldgroups=defaultdict(list)
        for v in old:oldgroups[tuple(map(int,X[v,24:]))].append(int(v))
        for key,inds in groups.items():
            r=len(upper);rows.extend([r]*len(inds));cols.extend(inds);values.extend([1]*len(inds))
            if key in oldgroups:
                vs=oldgroups[key];rows.extend([r]*len(vs));cols.extend([oldmap[v] for v in vs]);values.extend([-1]*len(vs))
                # All current tails have a complete16-point fiber; undeletable
                # points still count, so the new-minus-deleted upper bound is0.
                upper.append(0)
            else:upper.append(16)
        matrix=coo_matrix((values,(rows,cols)),shape=(len(upper),n+len(old))).tocsc()
        ans=milp(np.r_[-np.ones(n),np.ones(len(old))],integrality=np.ones(n+len(old)),bounds=Bounds(0,1),
                 constraints=LinearConstraint(matrix,-np.inf,upper),
                 options={'time_limit':args.seconds,'mip_rel_gap':0,'threads':2})
        add=np.array([],np.int64) if ans.x is None else np.flatnonzero(ans.x[:n]>.5)
        remove=np.array([],np.int64) if ans.x is None else old[np.flatnonzero(ans.x[n:]>.5)]
        gain=len(add)-len(remove)
        rec=dict(model=m,tail_indices=list(tailset),tail_vectors=oldtail[list(tailset)].tolist(),
                 candidates=n,old_variables=len(old),new_conflicts=new_conflicts,status=int(ans.status),
                 solver_message=str(ans.message),gain=gain,added=len(add),removed=len(remove),
                 numerical_gain_upper_bound=None if getattr(ans,'mip_dual_bound',None) is None else -float(ans.mip_dual_bound),
                 seconds=time.monotonic()-at)
        if gain>0:
            A=np.zeros((len(add),30),np.int16);A[:,NORMAL]=Y[add,:8];A[:,24:]=Y[add,8:]
            keep=np.ones(len(X),bool);keep[remove]=False;F=np.r_[X[keep],A]
            assert len(F)==len(X)+gain and len(np.unique(F,axis=0))==len(F)
            for j,a in enumerate(A):
                dots=np.zeros(len(F),np.int32)
                for k in np.flatnonzero(a):dots+=F[:,k].astype(np.int32)*int(a[k])
                dots[len(X)-len(remove)+j]=-32
                assert dots.max()<=16
            np.savez_compressed(args.out/f'exchange30_gain{gain}_model{m}.npz',removed=X[remove],added=A)
            np.savez_compressed(args.out/f'proposal30_{len(F)}.npz',points=F,squared_norm=32,inner_product_limit=16)
            rec['exact_changed_pairs']='PASS';best=max(best,gain)
        records.append(rec);print(json.dumps(rec),flush=True)
        report=dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),source_count=len(X),
                    pool_sha256=hashlib.sha256(args.pool.read_bytes()).hexdigest(),models_completed=len(records),
                    best_gain=best,seed=args.seed,seconds=time.monotonic()-started,records=records,
                    scope='Numerical bounded local MILPs; only direct-verified positive exchanges certify improvements.')
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        if best>0:break


if __name__=='__main__':main()
