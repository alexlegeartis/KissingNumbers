"""Joint shared-deletion selection for an explicitly priced candidate pool."""
from pathlib import Path
import argparse, json, time
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path)
    ap.add_argument('--seconds',type=float,default=180);args=ap.parse_args();p=args.directory;start=time.monotonic()
    meta=json.loads((p/'report.json').read_text());source=Path(meta['source'])
    z=np.load(source);X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:X//=2
    z=np.load(p/'pool.npz');C=z['candidates'].astype(np.int16);ptr=z['blocker_indptr'];old=z['blocker_indices']
    unique,inverse=np.unique(old,return_inverse=True);n=len(C);m=len(unique)
    if not n:raise ValueError('Empty pool')
    rr=[];cc=[];vv=[];upper=[]
    # All old-point deletion costs are shared, irrespective of the old layer.
    e=len(old);rows=np.arange(e,dtype=np.int32);candidate=np.repeat(np.arange(n,dtype=np.int32),np.diff(ptr))
    rr.extend([rows,rows]);cc.extend([candidate,n+inverse]);vv.extend([np.ones(e,np.float64),-np.ones(e,np.float64)]);upper.extend(np.zeros(e))
    D=C.astype(np.float32)@C.astype(np.float32).T
    conflict=D>16;np.fill_diagonal(conflict,False)
    i,j=np.nonzero(np.triu(conflict,1));e2=len(i);rows=np.arange(e,e+e2,dtype=np.int32)
    rr.extend([rows,rows]);cc.extend([i,j]);vv.extend([np.ones(e2),np.ones(e2)]);upper.extend(np.ones(e2))
    for vectors in (C[:,:24],C[:,24:]):
        B,bi=np.unique(vectors,axis=0,return_inverse=True)
        for b in range(len(B)):
            norm=int(np.sum(B[b].astype(np.int64)**2))
            if norm<=16:continue
            cap=16//(norm-16)
            ids=np.flatnonzero(bi==b)
            if len(ids)>cap:
                rows=np.full(len(ids),len(upper),np.int32);rr.append(rows);cc.append(ids);vv.append(np.ones(len(ids)));upper.append(cap)
    # Valid mutual-conflict cliques strengthen the otherwise weak edge LP.
    rng=np.random.default_rng(20261112);seen=set();cliques=0
    for anchor in rng.permutation(n):
        ids=np.flatnonzero(conflict[anchor]);rng.shuffle(ids);members=[int(anchor)]
        for j in ids:
            if np.all(conflict[j,members]):members.append(int(j))
            if len(members)>=48:break
        if len(members)<3:continue
        key=tuple(sorted(members))
        if key in seen:continue
        seen.add(key);rr.append(np.full(len(members),len(upper),np.int32));cc.append(np.asarray(members));vv.append(np.ones(len(members)));upper.append(1);cliques+=1
    matrix=coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(len(upper),n+m)).tocsc()
    print(json.dumps({'candidates':n,'relevant_old_points':m,'candidate_conflicts':e2,'cliques':cliques,'constraints':len(upper),'build_seconds':time.monotonic()-start}),flush=True)
    result=milp(np.r_[-np.ones(n),np.ones(m)],integrality=np.ones(n+m),bounds=Bounds(0,1),
        constraints=LinearConstraint(matrix,-np.inf,np.asarray(upper)),options={'time_limit':args.seconds,'mip_rel_gap':0,'disp':True})
    selected=np.zeros(0,dtype=int) if result.x is None else np.flatnonzero(result.x[:n]>.5)
    # Ignore potentially unnecessary deletion variables; exact blocker union is
    # the deletion set used for both cardinality and output verification.
    removed=np.unique(np.concatenate([old[ptr[i]:ptr[i+1]] for i in selected])) if len(selected) else np.zeros(0,dtype=int)
    if len(selected) and np.any(np.triu(D[np.ix_(selected,selected)]>16,1)):raise ValueError('Candidate conflict')
    gain=len(selected)-len(removed);report={'solver_status':int(result.status),'message':str(result.message),'gain':gain,
        'added':len(selected),'removed':len(removed),'selection':selected.tolist(),'deletions':removed.tolist(),
        'numerical_gain_upper_bound':None if getattr(result,'mip_dual_bound',None) is None else -float(result.mip_dual_bound),
        'seconds':time.monotonic()-start,'scope':'Numerical selection; an improving output gets exact changed-pair replay, not a source proof.'}
    if gain>0:
        keep=np.ones(len(X),bool);keep[removed]=False;F=np.r_[X[keep],C[selected]]
        if len(np.unique(F,axis=0))!=len(F):raise ValueError('Duplicate output')
        for i,row in enumerate(C[selected]):
            products=np.sum(F.astype(np.int64)*row.astype(np.int64),axis=1);products[len(X)-len(removed)+i]=-32
            if products.max(initial=-32)>16:raise ValueError('Exact output pair failure')
        np.savez_compressed(p/'improved_proposal.npz',points=F,squared_norm=32,inner_product_limit=16,coordinate_denominator=1)
        report['exact_changed_pairs']='PASS'
    (p/'selection.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
