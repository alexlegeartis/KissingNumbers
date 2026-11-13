#!/usr/bin/env python3
"""Exact partial-fiber search in a 27D E8-core x D3-tail candidate family.

Enumerates all 2,160 integral norm-16 E8 heads times 12 norm-16 tails
sqrt(2)*(±2,±2,0) and permutations. Source comparisons are exact in Q(sqrt2);
candidate-candidate comparisons are integral. Arbitrary partial head fibers
and shared incumbent deletions are allowed. This is a finite candidate family,
not a global kissing-number search.
"""
import argparse, hashlib, json, time
from itertools import combinations, product
from pathlib import Path
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

CORE=np.array([0,2,4,5,6,10,11,23],dtype=np.int64)
TAIL=np.array([24,25,26],dtype=np.int64)


def greater2(p,q,limit=16):
    """Exact test p + sqrt(2)q > limit for integer arrays."""
    p=np.asarray(p,dtype=np.int64); q=np.asarray(q,dtype=np.int64); r=limit-p
    return ((q>=0)&((r<0)|(2*q*q>r*r))) | ((q<0)&(r<0)&(r*r>2*q*q))


def heads2160():
    out=[]
    for i in range(8):
        for s in (-4,4):
            x=np.zeros(8,np.int16); x[i]=s; out.append(x)
    for I in combinations(range(8),4):
        for signs in product((-2,2),repeat=4):
            x=np.zeros(8,np.int16); x[list(I)]=signs; out.append(x)
    for i in range(8):
        for signs in product((-1,1),repeat=8):
            if sum(s<0 for s in signs)%2==0: continue
            x=np.asarray(signs,np.int16); x[i]*=3; out.append(x)
    H=np.asarray(out,np.int16)
    assert H.shape==(2160,8) and np.all((H.astype(np.int64)**2).sum(1)==16)
    assert len(np.unique(H,axis=0))==2160
    return H


def tails12():
    out=[]
    for i,j in combinations(range(3),2):
        for a,b in product((-2,2),repeat=2):
            t=np.zeros(3,np.int16); t[i]=a; t[j]=b; out.append(t)
    T=np.asarray(out,np.int16)
    assert T.shape==(12,3) and np.all((T.astype(np.int64)**2).sum(1)==8)
    return T


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',type=Path,required=True); ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--cutoff',type=int,default=256); ap.add_argument('--seconds',type=float,default=180)
    ap.add_argument('--max-candidates',type=int,default=12000); ap.add_argument('--seed',type=int,default=20261114)
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=False); t0=time.monotonic()
    z=np.load(args.source,allow_pickle=False); X=z['points'].astype(np.int64)
    if int(z['coordinate_denominator'])!=1 or not np.all((X*X).sum(1)==32):
        raise ValueError("expected exact integral norm-32 source")
    J=np.r_[CORE,TAIL]
    profiles,inverse,weights=np.unique(X[:,J],axis=0,return_inverse=True,return_counts=True)
    H=heads2160(); T=tails12(); P=H.astype(np.int64)@profiles[:,:8].astype(np.int64).T
    cand=[]; bl=[]; costs=[]
    for ti,t in enumerate(T):
        q=t.astype(np.int64)@profiles[:,8:].astype(np.int64).T
        bad=greater2(P,q,16); price=bad@weights.astype(np.int64)
        for hi in np.flatnonzero(price<=args.cutoff):
            cand.append((ti,int(hi))); bl.append(np.flatnonzero(bad[hi])); costs.append(int(price[hi]))
    if not cand: raise RuntimeError("no candidates under cutoff; raise --cutoff")
    cand=np.asarray(cand,np.int32); costs=np.asarray(costs,np.int32)
    if len(cand)>args.max_candidates:
        rng=np.random.default_rng(args.seed); jitter=rng.random(len(costs))*1e-6
        take=np.argsort(costs+jitter)[:args.max_candidates]; take.sort()
        cand=cand[take]; costs=costs[take]; bl=[bl[i] for i in take]
    hc=H[cand[:,1]].astype(np.int64); tc=T[cand[:,0]].astype(np.int64)
    assert np.all((hc*hc).sum(1)+2*(tc*tc).sum(1)==32)
    # Mandatory exact replay of every candidate's complete blocker profile set.
    for i in range(len(cand)):
        actual=np.flatnonzero(greater2(hc[i]@profiles[:,:8].astype(np.int64).T,
                                      tc[i]@profiles[:,8:].astype(np.int64).T,16))
        if not np.array_equal(actual,bl[i]): raise RuntimeError(f"blocker replay failed at {i}")
    # Exact new-new conflicts: h.h' + 2 t.t' > 16.
    ei=[]; ej=[]
    for lo in range(0,len(cand),128):
        D=hc[lo:lo+128]@hc.T+2*(tc[lo:lo+128]@tc.T)
        for r in range(len(D)):
            i=lo+r; js=np.flatnonzero(D[r,i+1:]>16)+i+1
            ei.extend([i]*len(js)); ej.extend(js.tolist())
    edges=np.c_[ei,ej].astype(np.int32).reshape((-1,2))
    used=np.unique(np.concatenate(bl)) if bl else np.empty(0,np.int64)
    vm={int(v):len(cand)+i for i,v in enumerate(used)}
    rr=[];cc=[];vv=[];ub=[]
    for i,bs in enumerate(bl):
        for b in bs:
            r=len(ub); rr.extend((r,r)); cc.extend((i,vm[int(b)])); vv.extend((1,-1)); ub.append(0)
    for i,j in edges:
        r=len(ub); rr.extend((r,r)); cc.extend((int(i),int(j))); vv.extend((1,1)); ub.append(1)
    M=coo_matrix((np.asarray(vv,np.float64),(rr,cc)),shape=(len(ub),len(cand)+len(used))).tocsc()
    obj=np.r_[-np.ones(len(cand)),weights[used].astype(np.float64)]
    report=dict(dimension=27,source_count=len(X),profiles=len(profiles),head_count=len(H),tail_count=len(T),
        raw_candidates=len(H)*len(T),admitted_candidates=len(cand),blocker_cutoff=args.cutoff,
        candidate_conflicts=len(edges),used_blocker_profiles=len(used),
        source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        candidate_family="2160 integral norm-16 E8 heads x 12 sqrt(2)-scaled D3 norm-16 tails",
        exact_norm_and_blocker_replay="PASS",
        scope="Finite algebraic candidate family; cutoff and pool cap are explicit")
    np.savez_compressed(args.out/'pool.npz',heads=hc.astype(np.int16),tails=tc.astype(np.int16),
        blocker_profile_ids=np.concatenate(bl),blocker_indptr=np.r_[0,np.cumsum([len(x) for x in bl])],
        blocker_costs=costs,edges=edges,profiles=profiles,profile_weights=weights)
    (args.out/'pool_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    res=milp(obj,integrality=np.ones(len(obj)),bounds=Bounds(0,1),
        constraints=LinearConstraint(M,-np.inf,np.asarray(ub)),
        options={'time_limit':args.seconds,'presolve':True,'mip_rel_gap':0.0})
    report.update(solver_status=int(res.status),solver_message=str(res.message),
        seconds=time.monotonic()-t0,dual_bound=None if res.mip_dual_bound is None else float(-res.mip_dual_bound))
    if res.x is not None:
        chosen=np.flatnonzero(res.x[:len(cand)]>.5)
        deleted_profiles=used[np.flatnonzero(res.x[len(cand):]>.5)]
        required=np.unique(np.concatenate([bl[i] for i in chosen])) if len(chosen) else np.empty(0,np.int64)
        if not set(map(int,required)).issubset(set(map(int,deleted_profiles))):
            raise RuntimeError("MILP omitted selected-candidate blockers")
        if len(chosen)>1:
            D=hc[chosen]@hc[chosen].T+2*(tc[chosen]@tc[chosen].T)
            if np.any(np.triu(D,1)>16): raise RuntimeError("selected candidates conflict")
        deleted=np.flatnonzero(np.isin(inverse,deleted_profiles)); gain=len(chosen)-len(deleted)
        # Exact replay against every retained projected source profile.
        retained_profiles=np.setdiff1d(np.arange(len(profiles)),deleted_profiles,assume_unique=False)
        for k in chosen:
            if np.any(greater2(hc[k]@profiles[retained_profiles,:8].astype(np.int64).T,
                               tc[k]@profiles[retained_profiles,8:].astype(np.int64).T,16)):
                raise RuntimeError("selected point conflicts with a retained source profile")
        report.update(selected=len(chosen),deleted=len(deleted),gain=int(gain),
            selected_indices=chosen.tolist(),deleted_profile_ids=deleted_profiles.tolist())
        if gain>0:
            A=np.zeros((len(chosen),27),np.int16); B=np.zeros_like(A)
            A[:,CORE]=hc[chosen].astype(np.int16); B[:,TAIL]=tc[chosen].astype(np.int16)
            keep=np.ones(len(X),bool); keep[deleted]=False
            np.savez_compressed(args.out/f"proposal27_{len(X)+gain}.npz",
                rational=np.vstack([X[keep].astype(np.int16),A]),
                radical=np.vstack([np.zeros((keep.sum(),27),np.int16),B]),
                radicand=2,squared_norm=32,inner_product_limit=16)
            report['final_geometry_replay']="REQUIRED: independently verify complete Q(sqrt(2)) configuration"
    (args.out/'solver_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('selected_indices','deleted_profile_ids')}),flush=True)

if __name__=='__main__': main()
