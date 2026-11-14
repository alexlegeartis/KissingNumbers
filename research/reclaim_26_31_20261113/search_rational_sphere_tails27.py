#!/usr/bin/env python3
"""27D search with a rationally dense family of radius-4 S^2 tails.

Each point is (h,t), ||h||^2=||t||^2=16. h ranges over all 2160
integral norm-16 E8 heads. t ranges over rational stereographic points
on the radius-4 sphere, not just the 12 D3 roots. All source blockers and
candidate conflicts are checked by integer arithmetic after clearing
per-point denominators.
"""
import argparse, hashlib, json, math, time
from itertools import combinations, product
from pathlib import Path
import numpy as np
if not __debug__:raise RuntimeError('Run without -O: exact input checks require assertions.')
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix
from search_sources import load_source27

CORE=np.array([0,2,4,5,6,10,11,23])
TAIL=np.array([24,25,26])

def heads():
    H=[]
    for i in range(8):
        for s in (-4,4):
            x=np.zeros(8,np.int16); x[i]=s; H.append(x)
    for I in combinations(range(8),4):
        for s in product((-2,2),repeat=4):
            x=np.zeros(8,np.int16); x[list(I)]=s; H.append(x)
    for i in range(8):
        for s in product((-1,1),repeat=8):
            if sum(v<0 for v in s)%2==0: continue
            x=np.asarray(s,np.int16); x[i]*=3; H.append(x)
    H=np.asarray(H,np.int64)
    assert H.shape==(2160,8) and np.all((H*H).sum(1)==16)
    return H

def tails(K):
    out={}
    for c in range(1,K+1):
        for a in range(-K,K+1):
            for b in range(-K,K+1):
                if math.gcd(math.gcd(abs(a),abs(b)),c)!=1: continue
                d=c*c+a*a+b*b
                v=(8*a*c,8*b*c,4*(c*c-a*a-b*b))
                g=math.gcd(math.gcd(math.gcd(abs(v[0]),abs(v[1])),abs(v[2])),d)
                out[(tuple(x//g for x in v),d//g)]=1
    T=sorted(out)
    V=np.asarray([x[0] for x in T],np.int64); D=np.asarray([x[1] for x in T],np.int64)
    assert np.all((V*V).sum(1)==16*D*D)
    return V,D

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',type=Path,required=True); ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--K',type=int,default=2); ap.add_argument('--blocker-cost',type=int,default=512)
    ap.add_argument('--max-candidates',type=int,default=10000); ap.add_argument('--seconds',type=float,default=180)
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=False); t0=time.monotonic()
    X=load_source27(a.source)
    J=np.r_[CORE,TAIL]; P,inv,w=np.unique(X[:,J],axis=0,return_inverse=True,return_counts=True)
    PC=P[:,:8]; PT=P[:,8:]; H=heads(); T,D=tails(a.K)
    hs=[]; ts=[]; ds=[]; bl=[]; prices=[]
    for t,d in zip(T,D):
        dot=H@PC.T*int(d)+t@PT.T
        bad=dot>16*int(d); cost=bad@w
        for hi in np.flatnonzero(cost<=a.blocker_cost):
            hs.append(H[hi]); ts.append(t); ds.append(d); bl.append(np.flatnonzero(bad[hi])); prices.append(int(cost[hi]))
    if not hs: raise RuntimeError('empty candidate pool; raise blocker cost or K')
    Hc=np.asarray(hs,np.int64); Tc=np.asarray(ts,np.int64); Dc=np.asarray(ds,np.int64); prices=np.asarray(prices)
    if len(Hc)>a.max_candidates:
        take=np.sort(np.lexsort((np.arange(len(prices)),prices))[:a.max_candidates])
        Hc,Tc,Dc,prices=Hc[take],Tc[take],Dc[take],prices[take]; bl=[bl[i] for i in take]
    # Replay exact blockers.
    for i in range(len(Hc)):
        actual=np.flatnonzero((Hc[i]@PC.T)*Dc[i]+Tc[i]@PT.T>16*Dc[i])
        if not np.array_equal(actual,bl[i]): raise RuntimeError('blocker replay mismatch')
    # Exact pair conflicts: (h_i.h_j-16)d_i d_j + t_i.t_j > 0.
    ei=[]; ej=[]
    for lo in range(0,len(Hc),96):
        HH=Hc[lo:lo+96]@Hc.T; TT=Tc[lo:lo+96]@Tc.T
        for r in range(len(HH)):
            i=lo+r; v=(HH[r]-16)*Dc[i]*Dc+TT[r]
            js=np.flatnonzero(v[i+1:]>0)+i+1; ei.extend([i]*len(js)); ej.extend(js.tolist())
    ei=np.asarray(ei,np.int32); ej=np.asarray(ej,np.int32)
    used=np.unique(np.concatenate(bl)); vm={int(v):len(Hc)+i for i,v in enumerate(used)}
    rr=[]; cc=[]; vv=[]; ub=[]
    for i,bs in enumerate(bl):
        for b in bs:
            r=len(ub); rr.extend((r,r)); cc.extend((i,vm[int(b)])); vv.extend((1,-1)); ub.append(0)
    for i,j in zip(ei,ej):
        r=len(ub); rr.extend((r,r)); cc.extend((int(i),int(j))); vv.extend((1,1)); ub.append(1)
    M=coo_matrix((np.asarray(vv,float),(rr,cc)),shape=(len(ub),len(Hc)+len(used))).tocsc()
    obj=np.r_[-np.ones(len(Hc)),w[used].astype(float)]
    report=dict(source_count=len(X),profiles=len(P),heads=len(H),tails=len(T),raw_pairs=len(H)*len(T),
        candidates=len(Hc),conflicts=len(ei),used_profiles=len(used),K=a.K,blocker_cost=a.blocker_cost,
        source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest(),
        exact_norms='PASS',exact_blocker_replay='PASS',
        tail_family='rational stereographic radius-4 S^2; core and tail squared norms 16',
        scope='Finite rational-tail catalogue, not a global exclusion')
    np.savez_compressed(a.out/'pool.npz',heads=Hc.astype(np.int16),tail_num=Tc,tail_den=Dc,
        blocker_ids=np.concatenate(bl),blocker_ptr=np.r_[0,np.cumsum([len(x) for x in bl])],
        conflicts=np.c_[ei,ej],prices=prices,profile_weights=w,profiles=P)
    (a.out/'pool_report.json').write_text(json.dumps(report,indent=2)+'\n')
    res=milp(obj,integrality=np.ones(len(obj)),bounds=Bounds(0,1),
        constraints=LinearConstraint(M,-np.inf,np.asarray(ub)),
        options={'time_limit':a.seconds,'presolve':True,'mip_rel_gap':0.0,'threads':1})
    report.update(status=int(res.status),message=str(res.message),seconds=time.monotonic()-t0,
        dual_bound=None if res.mip_dual_bound is None else float(-res.mip_dual_bound))
    if res.x is not None:
        chosen=np.flatnonzero(res.x[:len(Hc)]>.5); dels=used[np.flatnonzero(res.x[len(Hc):]>.5)]
        req=np.unique(np.concatenate([bl[i] for i in chosen])) if len(chosen) else np.empty(0,int)
        if not set(map(int,req)).issubset(set(map(int,dels))): raise RuntimeError('missing blocker')
        dels=req
        if len(chosen)<int(w[dels].sum()):
            chosen=np.empty(0,np.int64);dels=np.empty(0,np.int64)
        if len(chosen)>1:
            for pos,i in enumerate(chosen):
                v=(Hc[i]@Hc[chosen[pos+1:]].T-16)*Dc[i]*Dc[chosen[pos+1:]]+Tc[i]@Tc[chosen[pos+1:]].T
                if np.any(v>0): raise RuntimeError('selected pair conflict')
        removed=np.flatnonzero(np.isin(inv,dels)); gain=len(chosen)-len(removed)
        report.update(selected=len(chosen),removed=len(removed),gain=int(gain))
        if gain>0:
            np.savez_compressed(a.out/f'proposal27_{len(X)+gain}.npz',
                source_sha256=np.asarray(report['source_sha256']),added_heads=Hc[chosen].astype(np.int16),
                added_tail_num=Tc[chosen],added_tail_den=Dc[chosen],removed_source_indices=removed)
    (a.out/'solver_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('selected_indices',)}),flush=True)

if __name__=='__main__': main()
