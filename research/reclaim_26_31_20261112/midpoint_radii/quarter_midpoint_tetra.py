"""Research pricing of tetrahedral lifts of quarter-angle Leech edges.

For ||u||^2=||v||^2=32, u.v=8, h=(u+v)/2 has norm20.
Four norm12 tetrahedral tails give four norm32 points over h.  Every
other Leech point clears h; only the two parents can be core blockers.
This script prices complete four-point bundles against a supplied source.
It does not assert any bundles can coexist or any new kissing bound.
"""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
import argparse, hashlib, json, time
import numpy as np


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--leech',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--bodies',type=int,default=2048)
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic(); rng=np.random.default_rng(26111208)
    z=np.load(args.source,allow_pickle=False); X=z['points'].astype(np.int16)
    norm=int(z['squared_norm'])
    if norm==32:X*=2
    elif norm!=128:raise ValueError('source normalization')
    if not np.all(np.sum(X.astype(np.int64)**2,axis=1)==128):raise ValueError('source norm')
    L=np.load(args.leech,allow_pickle=False).astype(np.int16)
    if L.shape!=(196560,24) or not np.all(np.sum(L.astype(np.int64)**2,axis=1)==32):raise ValueError('Leech data')
    q=np.sum(X[:,24:].astype(np.int64)**2,axis=1)
    eq={r.tobytes() for r in X[q==0,:24]}
    missing=np.asarray([i for i,r in enumerate(L) if (2*r).tobytes() not in eq],np.int32)
    bodies=[]; parents=[]; origins=[]; seen=set()
    for label,ids,target in [('missing',missing,args.bodies//2),('full',np.arange(len(L)),args.bodies-args.bodies//2)]:
        attained=0; LF=L[ids].astype(np.float32)
        while attained<target:
            chosen=rng.choice(ids,min(32,target-attained),replace=False)
            dots=L[chosen].astype(np.float32)@LF.T
            for i,row in zip(chosen,dots):
                js=ids[np.flatnonzero(row==8)]
                if len(js)==0:continue
                j=int(rng.choice(js)); h=L[i]+L[j]
                if h.tobytes() in seen:continue
                if int(np.dot(h.astype(np.int32),h.astype(np.int32)))!=80:raise ValueError('midpoint norm')
                seen.add(h.tobytes()); bodies.append(h); parents.append([int(i),j]); origins.append(label); attained+=1
    B=np.asarray(bodies,np.int16); k=X.shape[1]-24
    tetra=[]; keys=[]
    for support in combinations(range(k),3):
        for parity in (0,1):
            tails=[]
            for signs in product((-4,4),repeat=3):
                if sum(v<0 for v in signs)%2!=parity:continue
                t=np.zeros(k,np.int16);t[list(support)]=signs;tails.append(t)
            tails=np.asarray(tails,np.int16)
            G=tails.astype(np.int32)@tails.astype(np.int32).T
            if not np.array_equal(G,64*np.eye(4,dtype=np.int32)-16):raise ValueError('tetrahedron')
            tetra.append(tails);keys.append([list(support),parity])
    T=np.asarray(tetra,np.int16)
    S,si=np.unique(X[:,24:],axis=0,return_inverse=True)
    sq=np.sum(S.astype(np.int32)**2,axis=1)
    maxima=(T.reshape(-1,k).astype(np.int32)@S.astype(np.int32).T).reshape(len(T),4,len(S)).max(axis=1)
    thresholds=64-maxima; offset=128;width=257
    if np.any(np.abs(thresholds)>=offset):raise ValueError('histogram endpoint')
    lookup=thresholds+offset; Sidx=np.arange(len(S))[None,:]
    prices=np.empty((len(B),len(T)),np.int32)
    examples=[]; XF=X[:,:24].astype(np.float32)
    for at in range(0,len(B),16):
        dots=B[at:at+16].astype(np.float32)@XF.T
        for loc,row in enumerate(dots):
            i=at+loc
            H=np.bincount(si*width+row.astype(np.int32)+offset,minlength=len(S)*width).reshape(len(S),width)
            greater=np.cumsum(H[:,::-1],axis=1)[:,::-1]-H
            counts=greater[Sidx,lookup]; prices[i]=counts.sum(axis=1)
            if i%32==0:
                best=int(np.argmin(prices[i]));examples.append({'body':i,'origin':origins[i],'tetrahedron':best,'blockers':int(prices[i,best]),'by_tail_norm':{str(t):int(counts[best,sq==t].sum()) for t in np.unique(sq)}})
        if at%256==0:print(json.dumps({'priced_bodies':min(at+16,len(B)),'of':len(B),'seconds':time.monotonic()-started}),flush=True)
    np.savez_compressed(args.out/'prices.npz',bodies=B,tetrahedra=T,parents=np.asarray(parents,np.int32),blocker_counts=prices)
    bestids=np.argsort(prices.ravel(),kind='stable')[:128]
    replay=[]
    for flat in bestids:
        bi,ti=divmod(int(flat),len(T)); C=np.concatenate([np.repeat(B[bi:bi+1],4,axis=0),T[ti]],axis=1)
        if not np.all(np.sum(C.astype(np.int64)**2,axis=1)==128):raise ValueError('new norm')
        union=np.zeros(len(X),bool)
        for c in C:
            dot=np.zeros(len(X),np.int32)
            for j in np.flatnonzero(c):dot+=X[:,j].astype(np.int32)*int(c[j])
            union|=dot>64
        if int(union.sum())!=int(prices[bi,ti]):raise ValueError('independent blocker mismatch')
        replay.append({'body':bi,'tetrahedron':ti,'blockers':np.flatnonzero(union).tolist()})
    report={'source':str(args.source.resolve()),'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),'dimension':X.shape[1],'baseline':len(X),'body_count':len(B),'bundle_count':int(prices.size),'points_per_bundle':4,'tetrahedra':keys,'missing_core_points':len(missing),'minimum_blockers':int(prices.min()),'price_quantiles':np.quantile(prices,[0,.01,.1,.5,.9,1]).tolist(),'counts':dict(Counter(map(int,prices.ravel()))),'examples':examples,'independent_replays':replay,'seconds':time.monotonic()-started,'scope':'Sampled quarter-angle midpoint bodies; all coordinate tetrahedra; exact prices with independent direct-integer checks of the cheapest128bundles. Cross-bundle compatibility not asserted.'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['counts','examples','tetrahedra','independent_replays']},indent=2),flush=True)


if __name__=='__main__':main()
