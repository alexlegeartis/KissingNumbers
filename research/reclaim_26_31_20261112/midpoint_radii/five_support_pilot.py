"""Price a new exact norm-twelve tail alphabet against the current witness.

This is a bounded research experiment, not a kissing-number certificate.  Every
old point is eligible for deletion.  Histograms factor out repeated old tails,
so a complete signed normal-octad body catalogue is cheap to price.  Norm32 is
used throughout; core bodies have shape (2,2,2,2,2), and tails (2,2,2), with arbitrary signs.
"""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
import argparse, hashlib, json, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("source",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--body-samples",type=int,default=256)
    ap.add_argument("--max-blockers",type=int,default=64)
    ap.add_argument("--pool-limit",type=int,default=6000)
    args=ap.parse_args();start=time.monotonic();args.out.mkdir(parents=True,exist_ok=False)
    z=np.load(args.source,allow_pickle=False);X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:
        if np.any(X%2):raise ValueError("Nonintegral norm32 source")
        X=X//2
    n,d=X.shape;k=d-24
    if not np.all(np.sum(X.astype(np.int64)**2,axis=1)==32):raise ValueError("Source norm")
    if k<3:raise ValueError("Tail shape needs three coordinates")
    T=[]
    for support in combinations(range(k),3):
        for signs in product((-2,2),repeat=3):
            t=np.zeros(k,np.int16);t[list(support)]=signs;T.append(t)
    T=np.asarray(T,np.int16)
    if not np.all(np.sum(T*T,axis=1)==12):raise ValueError("Tail norm")
    normal=[0,2,4,5,6,10,11,23];B=[]
    for support in combinations(normal,5):
        for signs in product((-2,2),repeat=5):
            b=np.zeros(24,np.int16);b[list(support)]=signs;B.append(b)
    q=np.sum(X[:,24:].astype(np.int64)**2,axis=1)
    incumbent_bodies=np.unique(X[q==12,:24],axis=0)
    rng=np.random.default_rng(20261112)
    extra=incumbent_bodies[rng.choice(len(incumbent_bodies),min(args.body_samples,len(incumbent_bodies)),replace=False)]
    B=np.unique(np.vstack([B,extra]),axis=0)
    S,si=np.unique(X[:,24:],axis=0,return_inverse=True)
    sq=np.sum(S.astype(np.int64)**2,axis=1)
    dotT=T.astype(np.int32)@S.astype(np.int32).T
    # Every source core dot lies within [-40,40] by Cauchy--Schwarz.
    limit=16-dotT;lookup=np.clip(limit+40,0,80)
    if np.any((limit<=-40)|(limit>=40)):raise ValueError("Histogram endpoint")
    Sidx=np.arange(len(S))[None,:]
    allcosts=np.empty((len(B),len(T)),np.int32)
    min_by_body=[];witnesses=[]
    Xcore=X[:,:24].astype(np.float32)
    bound=24*int(np.abs(B).max())*int(np.abs(X).max())
    if bound>=2**24:raise ValueError("Exact integer accumulation bound")
    for at in range(0,len(B),16):
        bd=B[at:at+16].astype(np.float32)@Xcore.T
        if np.any(np.abs(bd)>40):raise ValueError("Histogram range")
        for local,delta in enumerate(bd):
            i=at+local
            H=np.bincount(si*81+delta.astype(np.int32)+40,minlength=len(S)*81).reshape(len(S),81)
            greater=np.cumsum(H[:,::-1],axis=1)[:,::-1]-H
            counts=greater[Sidx,lookup]
            allcosts[i]=counts.sum(axis=1)
            best=int(np.argmin(allcosts[i]));min_by_body.append(int(allcosts[i,best]))
            if i%64==0:
                witnesses.append({'body_index':i,'tail_index':best,'blockers':int(allcosts[i,best]),
                    'by_old_tail_norm':{str(v):int(counts[best,sq==v].sum()) for v in (0,8,16,32)}})
        if at%256==0:print(json.dumps({'priced_bodies':min(at+16,len(B)),'of':len(B),'seconds':time.monotonic()-start}),flush=True)
    np.savez_compressed(args.out/'prices.npz',bodies=B,tails=T,blocker_counts=allcosts)
    ids=np.flatnonzero(allcosts.ravel()<=args.max_blockers)
    # Rank by exact deletion cost, then deterministic shuffled order for ties.
    rng.shuffle(ids);ids=ids[np.argsort(allcosts.ravel()[ids],kind='stable')]
    eligible=len(ids);ids=ids[:args.pool_limit]
    C=np.concatenate([B[ids//len(T)],T[ids%len(T)]],axis=1)
    oldset={row.tobytes() for row in X};keep=[i for i,row in enumerate(C) if row.tobytes() not in oldset]
    C=C[keep];ids=ids[keep]
    ptr=[0];edges=[]
    # Deliberately independent direct integer sums reconstruct every blocker.
    for row,flatid in zip(C,ids):
        dots=np.zeros(n,np.int32)
        for j in np.flatnonzero(row):dots+=X[:,j].astype(np.int32)*int(row[j])
        bad=np.flatnonzero(dots>16)
        if len(bad)!=allcosts.ravel()[flatid]:raise ValueError("Factored/direct blocker disagreement")
        edges.extend(bad.tolist());ptr.append(len(edges))
    ptr=np.asarray(ptr,np.int64);edges=np.asarray(edges,np.int64)
    graph=csr_matrix((np.ones(len(edges),np.int8),edges,ptr),shape=(len(C),n))
    owners=maximum_bipartite_matching(graph,perm_type='column') if len(C) else np.zeros(0,np.int64)
    matched=int(np.sum(owners>=0))
    if matched==len(C) and len(C):
        if len(set(owners.tolist()))!=len(C):raise ValueError("Owner injection")
        products=np.sum(C.astype(np.int64)*X[owners].astype(np.int64),axis=1)
        if np.any(products<=16):raise ValueError("Owner conflict")
    np.savez_compressed(args.out/'pool.npz',candidates=C,blocker_indptr=ptr,blocker_indices=edges,owners=owners)
    report={'dimension':d,'source':str(args.source.resolve()),'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),
        'source_count':n,'body_catalogue':'all signed five-subsets of normal octad plus deterministic random incumbent midpoint bodies',
        'normal_octad':normal,'body_samples':args.body_samples,'bodies':len(B),'tails':len(T),'distinct_old_tails':len(S),
        'priced_candidates':int(allcosts.size),'all_blocker_counts':dict(sorted(Counter(map(int,allcosts.ravel())).items())),
        'body_minimum_blocker_counts':dict(sorted(Counter(min_by_body).items())),
        'layer_examples':witnesses,'cutoff':args.max_blockers,'eligible_at_cutoff':eligible,'selected_pool':len(C),
        'pool_limit':args.pool_limit,'pool_edges':len(edges),'injectively_matched':matched,
        'minimum_cost':int(allcosts.min()),'zero_blocker_candidates':int(np.sum(allcosts==0)),
        'status':'bounded exact pricing; no geometric lower-bound claim',
        'pool_conclusion':'Every subset of selected finite pool has nonpositive gain by independent checked forbidden-owner injection' if matched==len(C) else 'Matching leaves candidates unmatched; selection remains open',
        'seconds':time.monotonic()-start,'partial':False}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('all_blocker_counts','body_minimum_blocker_counts','layer_examples')},indent=2),flush=True)

if __name__=='__main__':main()
