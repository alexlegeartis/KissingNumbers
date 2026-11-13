#!/usr/bin/env python3
"""Component-wise exact MILP for the 26D 759-octad exchange pool.

A candidate is selected only if every incumbent blocker is deleted. Selected
candidates must be mutually compatible. Objective: number added - number deleted.
Every blocker is replayed exactly in Q(sqrt(3)); every candidate conflict is
rebuilt with int64 products. A timeout is not an exclusion proof.
"""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components


def above16(a, b):
    """Exact predicate a + b*sqrt(3) > 16 for integer arrays."""
    a=np.asarray(a,dtype=np.int64); b=np.asarray(b,dtype=np.int64); c=a-16
    return ((b==0)&(c>0)) | ((b>0)&((c>=0)|(3*b*b>c*c))) | ((b<0)&(c>0)&(c*c>3*b*b))


def replay_blockers(x, A, B):
    da=np.zeros(len(A),np.int64); db=np.zeros(len(A),np.int64)
    for k in np.flatnonzero(x):
        da += int(x[k])*A[:,k]; db += int(x[k])*B[:,k]
    return np.flatnonzero(above16(da,db))


def conflict_edges(C):
    ei=[]; ej=[]; C=C.astype(np.int64)
    for lo in range(0,len(C),128):
        D=C[lo:lo+128]@C.T
        for r in range(len(D)):
            i=lo+r; js=np.flatnonzero(D[r,i+1:]>16)+i+1
            ei.extend([i]*len(js)); ej.extend(js.tolist())
    return np.asarray(ei,np.int32),np.asarray(ej,np.int32)


def components(nc, blockers, ei, ej):
    vids=np.unique(np.concatenate(blockers)) if blockers else np.empty(0,np.int64)
    vm={int(v):nc+i for i,v in enumerate(vids)}
    rr=[];cc=[]
    for i,bs in enumerate(blockers):
        for v in bs:
            j=vm[int(v)]; rr.extend((i,j)); cc.extend((j,i))
    for i,j in zip(ei,ej):
        rr.extend((int(i),int(j))); cc.extend((int(j),int(i)))
    G=csr_matrix((np.ones(len(rr),np.int8),(rr,cc)),shape=(nc+len(vids),)*2)
    nlab,lab=connected_components(G,directed=False,return_labels=True)
    out=[]
    for k in range(nlab):
        cs=np.flatnonzero(lab[:nc]==k)
        if len(cs): out.append((cs,vids[np.flatnonzero(lab[nc:]==k)]))
    return out


def greedy(cs, blockers, edges, seed):
    """Multi-start feasible selection; shared blockers are naturally reused."""
    rng=np.random.default_rng(seed)
    sig={int(i):tuple(map(int,blockers[int(i)])) for i in cs}
    freq={}
    for i in cs:
        for v in sig[int(i)]: freq[v]=freq.get(v,0)+1
    reuse={int(i):sum(1/freq[v] for v in sig[int(i)]) for i in cs}
    orders=[sorted(cs,key=lambda i:(len(sig[int(i)]),-reuse[int(i)],int(i))),
            sorted(cs,key=lambda i:(-reuse[int(i)],len(sig[int(i)]),int(i)))]
    for _ in range(32):
        jitter=rng.random(len(cs)); orders.append(sorted(cs,key=lambda i:
            len(sig[int(i)])-.15*reuse[int(i)]+.25*jitter[np.searchsorted(cs,i)]))
    best=(0,[],[])
    for order in orders:
        chosen=[]; deleted=set(); used=set()
        for i0 in order:
            i=int(i0)
            if any((min(i,j),max(i,j)) in edges for j in used): continue
            chosen.append(i); used.add(i); deleted.update(sig[i])
        gain=len(chosen)-len(deleted)
        if gain>best[0]: best=(gain,chosen,sorted(deleted))
    return best


def solve_component(cs,vs,blockers,ei,ej,seconds):
    cm={int(x):i for i,x in enumerate(cs)}
    vm={int(x):len(cs)+i for i,x in enumerate(vs)}
    rr=[];cc=[];vv=[];ub=[]
    for c in cs:
        for v in blockers[int(c)]:
            r=len(ub);rr.extend((r,r));cc.extend((cm[int(c)],vm[int(v)]));vv.extend((1,-1));ub.append(0)
    for a,b in zip(ei,ej):
        r=len(ub);rr.extend((r,r));cc.extend((cm[int(a)],cm[int(b)]));vv.extend((1,1));ub.append(1)
    M=coo_matrix((vv,(rr,cc)),shape=(len(ub),len(cs)+len(vs))).tocsc()
    obj=np.r_[-np.ones(len(cs)),np.ones(len(vs))]
    res=milp(obj,integrality=np.ones(len(obj)),bounds=Bounds(0,1),
        constraints=LinearConstraint(M,-np.inf,np.asarray(ub)),
        options={'time_limit':seconds,'presolve':True,'mip_rel_gap':0.0})
    if res.x is None:
        return dict(status=int(res.status),message=str(res.message),gain=None,selected=None,deleted=None,dual_bound=None)
    sel=cs[np.flatnonzero(res.x[:len(cs)]>.5)].astype(int)
    dele=vs[np.flatnonzero(res.x[len(cs):]>.5)].astype(int)
    req=np.unique(np.concatenate([blockers[int(i)] for i in sel])) if len(sel) else np.empty(0,int)
    if not set(map(int,req)).issubset(set(map(int,dele))): raise RuntimeError("MILP omitted a blocker")
    return dict(status=int(res.status),message=str(res.message),gain=int(len(sel)-len(dele)),
        selected=sel.tolist(),deleted=dele.tolist(),
        dual_bound=None if res.mip_dual_bound is None else float(-res.mip_dual_bound))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',type=Path,required=True); ap.add_argument('--pool',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True); ap.add_argument('--seconds',type=float,default=180)
    ap.add_argument('--seed',type=int,default=20261114); args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False); t0=time.monotonic()
    z=np.load(args.source,allow_pickle=False)
    assert int(z['rational_denominator'])==1
    assert np.all(z['quadratic_denominators']==1) and np.all(z['quadratic_radicands']==3)
    A=np.r_[z['rational_points'],z['quadratic_A']].astype(np.int64)
    B=np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']].astype(np.int64)
    p=np.load(args.pool,allow_pickle=False); C=p['candidate_A'].astype(np.int64)
    ptr,idx=p['blocker_indptr'],p['blocker_indices']
    supplied=[idx[ptr[i]:ptr[i+1]].astype(np.int64) for i in range(len(C))]
    if len(C)!=len(np.unique(C,axis=0)) or np.any((C*C).sum(1)!=32): raise ValueError("invalid pool")
    blockers=[]
    for i,x in enumerate(C):
        actual=replay_blockers(x,A,B)
        if not np.array_equal(np.sort(actual),np.sort(supplied[i])): raise ValueError(f"blocker mismatch {i}")
        blockers.append(actual)
    ei,ej=conflict_edges(C); edge_set={(int(a),int(b)) for a,b in zip(ei,ej)}
    comps=components(len(C),blockers,ei,ej)
    greedy_rows=[]
    for k,(cs,vs) in enumerate(comps):
        g=greedy(cs,blockers,edge_set,args.seed+k)
        greedy_rows.append(dict(candidates=len(cs),victims=len(vs),gain=g[0],selected=g[1],deleted=g[2]))
    rows=[]; solved={}; remaining=args.seconds
    for k,(cs,vs) in enumerate(comps):
        if remaining<=0: break
        mask=np.isin(ei,cs)&np.isin(ej,cs); t=time.monotonic()
        budget=max(.25,remaining*len(cs)/max(1,len(C)))
        r=solve_component(cs,vs,blockers,ei[mask],ej[mask],budget)
        r.update(component=k,candidates=len(cs),victims=len(vs),seconds=time.monotonic()-t)
        rows.append(r); solved[k]=r; remaining=max(0,args.seconds-(time.monotonic()-t0))
        print(json.dumps(r),flush=True)
    chosen=[]; deleted=[]
    for k,(cs,vs) in enumerate(comps):
        r=solved.get(k); g=greedy_rows[k]
        if r and r['selected'] is not None and r['gain']>=g['gain']:
            chosen.extend(r['selected']); deleted.extend(r['deleted'])
        else: chosen.extend(g['selected']); deleted.extend(g['deleted'])
    chosen=np.asarray(sorted(set(chosen)),np.int64); deleted=np.asarray(sorted(set(deleted)),np.int64)
    required=np.unique(np.concatenate([blockers[int(i)] for i in chosen])) if len(chosen) else np.empty(0,np.int64)
    if not set(map(int,required)).issubset(set(map(int,deleted))): raise RuntimeError("assembled exchange omits blockers")
    if len(chosen)>1 and np.any(np.triu(C[chosen]@C[chosen].T,1)>16): raise RuntimeError("candidate conflict")
    # Materialize and independently replay every changed pair. Unchanged
    # source-source pairs inherit the frozen source certificate.
    keep=np.ones(len(A),dtype=bool); keep[deleted]=False
    FA=np.vstack([A[keep],C[chosen]])
    FB=np.vstack([B[keep],np.zeros((len(chosen),A.shape[1]),dtype=np.int64)])
    if len(FA)!=len(A)+len(chosen)-len(deleted): raise RuntimeError("bad exchange cardinality")
    if not (np.all((FA*FA+3*FB*FB).sum(1)==32) and np.all((FA*FB).sum(1)==0)):\n        raise RuntimeError("proposal exact Q(sqrt(3)) norm failure")
    for k,x in enumerate(C[chosen]):
        da=np.zeros(len(FA),np.int64); db=np.zeros(len(FA),np.int64)
        for j in np.flatnonzero(x):
            da += int(x[j])*FA[:,j]; db += int(x[j])*FB[:,j]
        # The new point itself has product 32; exclude its own appended row.
        da[len(A)-len(deleted)+k]=-32; db[len(A)-len(deleted)+k]=0
        if np.any(above16(da,db)): raise RuntimeError("changed-pair replay failed")
    np.savez_compressed(args.out/'proposal_coordinates.npz',rational=FA,radical=FB,
        radicand=3,squared_norm=32,inner_product_limit=16)
    report=dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        pool_sha256=hashlib.sha256(args.pool.read_bytes()).hexdigest(),source_count=len(A),
        candidates=len(C),candidate_conflicts=len(ei),components=len(comps),
        exact_blocker_replay="PASS",greedy_components=greedy_rows,solved_components=rows,
        added=len(chosen),removed=len(deleted),gain=int(len(chosen)-len(deleted)),
        scope="Finite 2,004-candidate pool only. Positive exchange requires full final exact replay; timeout is not exclusion.",
        seconds=time.monotonic()-t0)
    np.savez_compressed(args.out/'proposal_indices.npz',selected_candidates=chosen,deleted_source_indices=deleted)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('greedy_components','solved_components')}),flush=True)

if __name__=='__main__': main()
