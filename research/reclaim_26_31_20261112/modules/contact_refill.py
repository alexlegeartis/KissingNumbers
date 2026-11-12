"""Negative-contact closure from explicit replacement anchors on a new seed."""
from pathlib import Path
import argparse, hashlib, json, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import milp, Bounds, LinearConstraint

def select(C,blockers,seconds):
    n=len(C);old=np.unique([j for b in blockers for j in b]);mapping={int(x):n+i for i,x in enumerate(old)}
    r=[];c=[];v=[];upper=[]
    for i,bs in enumerate(blockers):
        for j in bs:
            k=len(upper);r.extend([k,k]);c.extend([i,mapping[j]]);v.extend([1,-1]);upper.append(0)
    edges=0;F=C.astype(np.float32)
    for at in range(0,n,128):
        dots=F[at:at+128]@F.T
        for local in range(len(dots)):
            i=at+local
            for j0 in np.flatnonzero(dots[local,i+1:]>64)+i+1:
                j=int(j0);k=len(upper);r.extend([k,k]);c.extend([i,j]);v.extend([1,1]);upper.append(1);edges+=1
    matrix=coo_matrix((v,(r,c)),shape=(len(upper),n+len(old))).tocsc()
    result=milp(np.r_[-np.ones(n),np.ones(len(old))],integrality=np.ones(n+len(old)),bounds=Bounds(0,1),
        constraints=LinearConstraint(matrix,-np.inf,upper),options={'time_limit':seconds,'mip_rel_gap':0,'threads':2})
    add=[] if result.x is None else np.flatnonzero(result.x[:n]>.5).tolist()
    remove=sorted(set(j for i in add for j in blockers[i]))
    return {'gain':len(add)-len(remove),'added':add,'removed':remove,'status':int(result.status),'message':str(result.message),
        'candidate_conflicts':edges,'numerical_gain_bound':None if getattr(result,'mip_dual_bound',None) is None else -float(result.mip_dual_bound),
        'scope':'Numerical finite-pool selection; only exact changed-pair checks establish a positive output.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('anchors',type=Path)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--cutoff',type=int,default=2)
    ap.add_argument('--seconds',type=float,default=60);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    X=np.load(args.source)['points'].astype(np.int16);A=np.load(args.anchors)['added'].astype(np.int16)
    if not np.all(np.sum(X.astype(np.int64)**2,axis=1)==128):raise ValueError('Source normalization')
    if not np.all(np.sum(A.astype(np.int64)**2,axis=1)==128):raise ValueError('Anchor normalization')
    lookup={row.tobytes():i for i,row in enumerate(X)};anchors=np.asarray([lookup[row.tobytes()] for row in A])
    F=X.astype(np.float32);order=np.random.default_rng(20261112).permutation(len(X))
    known=set(lookup);seen=set();pool=[];blockers=[];records=[];incidences=checks=0
    if X.shape[1]*2*int(np.abs(X).max())**2>=2**24:raise ValueError('Exact binary32 bound')
    for at in range(0,len(anchors),32):
        ids=anchors[at:at+32];dot=F[ids]@F.T;ii,jj=np.nonzero(dot==-64);Y=X[ids[ii]]+X[jj];incidences+=len(Y)
        fresh=[]
        for row in Y:
            key=row.tobytes()
            if key not in known and key not in seen:seen.add(key);fresh.append(row)
        Y=np.asarray(fresh,np.int16).reshape((-1,X.shape[1]));active=np.arange(len(Y));counts=np.zeros(len(Y),np.int16)
        if not np.all(np.sum(Y.astype(np.int64)**2,axis=1)==128):raise ValueError('Contact-sum norm')
        for off in range(0,len(order),256):
            if not len(active):break
            bad=Y[active].astype(np.float32)@F[order[off:off+256]].T>64;checks+=bad.size
            counts[active]+=bad.sum(axis=1).astype(np.int16);active=active[counts[active]<=args.cutoff]
        bs=[]
        for row in Y[active]:
            values=np.zeros(len(X),np.int32)
            for j in np.flatnonzero(row):values+=X[:,j].astype(np.int32)*int(row[j])
            b=np.flatnonzero(values>64).tolist()
            if len(b)>args.cutoff:raise ValueError('Exact blocker replay')
            bs.append(b)
        pool.extend(Y[active]);blockers.extend(bs)
        rec={'anchor_start':at,'anchors':len(ids),'sum_incidences':len(ii),'new_sums':len(Y),'admitted':len(active),'free':sum(not b for b in bs),'seconds':time.monotonic()-start}
        records.append(rec);print(json.dumps(rec),flush=True);(args.out/'partial.json').write_text(json.dumps(records,indent=2)+'\n')
    C=np.asarray(pool,np.int16).reshape((-1,X.shape[1]));ptr=np.r_[0,np.cumsum(list(map(len,blockers)))];edges=np.asarray([j for b in blockers for j in b],np.int64)
    np.savez_compressed(args.out/'pool.npz',candidates=C,blocker_indptr=ptr,blocker_indices=edges)
    selection=select(C,blockers,args.seconds) if len(C) else {'gain':0,'added':[],'removed':[]}
    if selection['gain']>0:
        keep=np.ones(len(X),bool);keep[selection['removed']]=False;added=C[selection['added']];out=np.r_[X[keep],added]
        if len(np.unique(out,axis=0))!=len(out):raise ValueError('Output distinctness')
        for i,row in enumerate(added):
            values=np.zeros(len(out),np.int32)
            for j in np.flatnonzero(row):values+=out[:,j].astype(np.int32)*int(row[j])
            values[len(X)-len(selection['removed'])+i]=-128
            if values.max()>64:raise ValueError('Independent exact changed pair')
        selection['exact_changed_pairs']='PASS';selection['final_count']=len(out)
        np.savez_compressed(args.out/f'proposal31_{len(out)}.npz',points=out,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
        np.savez_compressed(args.out/'exchange.npz',removed=X[selection['removed']],added=added,coordinate_denominator=2,squared_norm=128,inner_product_limit=64)
    report={'source':str(args.source.resolve()),'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),'anchors_file':str(args.anchors.resolve()),
        'anchors_sha256':hashlib.sha256(args.anchors.read_bytes()).hexdigest(),'baseline':len(X),'anchors':len(anchors),
        'sum_incidences':incidences,'unique_new_sums':len(seen),'admitted':len(C),'cutoff':args.cutoff,'pricing_products':checks,
        'optimization':selection,'seconds':time.monotonic()-start,
        'scope':'Complete contact closure on the explicit added-head anchors, with all source points; only the stated blocker-cutoff selection pool.'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
