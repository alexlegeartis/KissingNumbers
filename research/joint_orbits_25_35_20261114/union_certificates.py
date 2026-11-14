"""Exact optimal union of two verified codes via bipartite matching.

This is discovery, with exact delta checks; every positive needs full final
verification. Matching proves optimality only within this two-code union.
"""
import argparse,collections,json,time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
if not __debug__:raise RuntimeError('Assertions are required.')
ROOT=Path(__file__).resolve().parents[2]
CACHE={}

def read(path):
    z=np.load(path,allow_pickle=False)
    A=z['points'] if 'points' in z else z['rational']
    B=z['radical'] if 'radical' in z else np.zeros_like(A)
    assert np.issubdtype(A.dtype,np.integer) and np.issubdtype(B.dtype,np.integer)
    assert A.shape==B.shape and A.ndim==2
    assert int(z['inner_product_limit'])*2==int(z['squared_norm'])
    scale={32:2,128:1}[int(z['squared_norm'])]
    return (scale*A).astype(np.int16),(scale*B).astype(np.int16),int(z['radicand']) if 'radicand' in z else 3

def above(p,q,k):
    p=p.astype(np.int64);q=q.astype(np.int64);t=64-p
    return ((q>=0)&((t<0)|(k*q*q>t*t)))|((q<0)&(t<0)&(t*t>k*q*q))

def keys(A,B):return [v.tobytes() for v in np.c_[A,B]]

def union(A,B,C,D,k):
    old=set(keys(A,B));keep=np.array([v not in old for v in keys(C,D)])
    C,D=C[keep],D[keep]
    if not len(C):return A,B,dict(new=0,gain=0)
    AF,BF=A.astype(np.float32),B.astype(np.float32);CF,DF=C.astype(np.float32),D.astype(np.float32)
    rows=[];cols=[];ckeys=[(key,k if np.any(d) else 0) for key,d in zip(keys(C,D),D)];cache=CACHE.setdefault((len(A),k if np.any(B) else 0),{})
    missing=np.array([i for i,key in enumerate(ckeys) if key not in cache],int)
    print('UNION_START',len(C),'uncached',len(missing),flush=True)
    for offset in range(0,len(missing),64):
        ids=missing[offset:offset+64];p=CF[ids]@AF.T
        if np.any(B) or np.any(DF[ids]):
            q=np.zeros_like(p)
            if np.any(B):q+=CF[ids]@BF.T
            if np.any(DF[ids]):
                q+=DF[ids]@AF.T
                if np.any(B):p+=k*(DF[ids]@BF.T)
            bad=above(p,q,k)
        else:bad=p>64
        for i,row in zip(ids,bad):cache[ckeys[i]]=np.flatnonzero(row).astype(np.int32)
    for i,key in enumerate(ckeys):
        js=cache[key];rows.extend([i]*len(js));cols.extend(js.tolist())
    print('EXACT_GRAPH',len(C),'edges',len(rows),flush=True)
    G=csr_matrix((np.ones(len(rows),np.int8),(rows,cols)),shape=(len(C),len(A)))
    match=maximum_bipartite_matching(G,perm_type='column');rev=np.full(len(A),-1,int)
    for i,j in enumerate(match):
        if j>=0:rev[j]=i
    cseen=match<0;oseen=np.zeros(len(A),bool);todo=collections.deque(np.flatnonzero(cseen))
    while todo:
        i=todo.popleft()
        for j in G.indices[G.indptr[i]:G.indptr[i+1]]:
            if j==match[i] or oseen[j]:continue
            oseen[j]=True;m=rev[j]
            assert m>=0,'An augmenting path contradicts maximum matching.'
            if not cseen[m]:cseen[m]=True;todo.append(int(m))
    assert int(cseen.sum()-oseen.sum())==int((match<0).sum())
    rec=dict(new=len(C),edges=G.nnz,added=int(cseen.sum()),removed=int(oseen.sum()),gain=int((match<0).sum()))
    if not rec['gain']:return A,B,rec
    print('POSITIVE_MATCHING',json.dumps(rec),flush=True)
    RA,RB=np.r_[A[~oseen],C[cseen]],np.r_[B[~oseen],D[cseen]]
    assert len(set(keys(RA,RB)))==len(RA)
    assert np.all(np.sum(RA.astype(np.int64)**2+k*RB.astype(np.int64)**2,axis=1)==128)
    assert not np.any(np.sum(RA.astype(np.int64)*RB.astype(np.int64),axis=1))
    for ca,cb in zip(C[cseen],D[cseen]):
        p=RA.astype(np.int64)@ca.astype(np.int64)+k*(RB.astype(np.int64)@cb.astype(np.int64))
        q=RA.astype(np.int64)@cb.astype(np.int64)+RB.astype(np.int64)@ca.astype(np.int64)
        bad=above(p,q,k);own=np.all(RA==ca,axis=1)&np.all(RB==cb,axis=1);bad[own]=False;assert not bad.any()
    return RA,RB,rec

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--start',type=Path);a=ap.parse_args();a.out.mkdir(exist_ok=False)
    files=sorted((ROOT/'verifications/improved').glob(f'**/coordinates{a.dimension}_*.npz'))
    loaded=[(f,*read(f)) for f in files]
    if a.start:loaded.append((a.start,*read(a.start)))
    loaded.sort(key=lambda v:len(v[1]),reverse=True)
    source,A,B,k=loaded[0];rows=[];start=time.monotonic();baseline=len(A)
    for f,C,D,kk in loaded[1:]:
        if np.any(B) and np.any(D) and kk!=k:continue
        usek=k if np.any(B) else kk
        A,B,row=union(A,B,C,D,usek);k=usek
        row.update(donor=str(f.relative_to(ROOT)),count=len(A));rows.append(row);print(json.dumps(row),flush=True)
        (a.out/'progress.json').write_text(json.dumps(dict(dimension=a.dimension,baseline=baseline,count=len(A),rows=rows),indent=2)+'\n')
        if row.get('gain',0)>0:
            np.savez_compressed(a.out/f'candidate{a.dimension}_{len(A)}.npz',rational=A,radical=B,radicand=k if np.any(B) else 3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
    report=dict(dimension=a.dimension,baseline=baseline,count=len(A),source=str(source),rows=rows,seconds=time.monotonic()-start)
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
