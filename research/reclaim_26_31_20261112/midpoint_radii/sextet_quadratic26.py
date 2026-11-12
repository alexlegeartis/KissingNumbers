"""Exact finite sextet unions over Q(sqrt(3)), using the latest 26D source."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
from scipy.sparse import csr_matrix, vstack
from scipy.sparse.csgraph import maximum_bipartite_matching

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'verifications/improved/dim27-31-joint-refills'))
from search_sextets import frames


def greater(p,q):
    p=p.astype(np.int64);q=q.astype(np.int64);z=64-p
    return np.where(q>=0,(z<0)|(3*q*q>z*z),(z<0)&(z*z>3*q*q))


def products(a,b):
    A=a[:,:26].astype(np.float32);B=a[:,26:].astype(np.float32)
    C=b[:,:26].astype(np.float32);D=b[:,26:].astype(np.float32)
    return greater(A@C.T+3*(B@D.T),A@D.T+B@C.T)


def transform(x,blocks):
    y=x.copy()
    for shift in (0,26):
        for block in blocks:
            b=np.asarray(block)+shift;s=x[:,b].sum(axis=1)
            if np.any(s%2):raise ValueError('fractional transform')
            y[:,b]=s[:,None]//2-x[:,b]
        y[:,np.asarray(blocks[0])+shift]*=-1
    A=y[:,:26].astype(np.int32);B=y[:,26:].astype(np.int32)
    if not np.all(np.sum(A*A+3*B*B,axis=1)==128) or np.any(np.sum(A*B,axis=1)):
        raise ValueError('transformed norm')
    return y


def refill(x,y,cutoff):
    old={r.tobytes():i for i,r in enumerate(x)};new={r.tobytes():i for i,r in enumerate(y)}
    common=old.keys()&new.keys()
    ai=np.asarray([i for k,i in old.items() if k not in common],np.int32)
    bi=np.asarray([i for k,i in new.items() if k not in common],np.int32)
    a=x[ai];b=y[bi]
    if not len(b):return x,{'gain':0,'exclusive':0,'candidates':0}
    mats=[];selected=[]
    for at in range(0,len(b),128):
        G=products(b[at:at+128],a);good=np.flatnonzero(G.sum(axis=1)<=cutoff)
        if len(good):mats.append(csr_matrix(G[good]));selected.extend((good+at).tolist())
    if not selected:return x,{'gain':0,'exclusive':len(b),'candidates':0}
    exclusive=len(b);b=b[selected];G=vstack(mats,format='csr')
    matching=maximum_bipartite_matching(G,perm_type='column');gain=len(b)-int(np.count_nonzero(matching>=0))
    if gain==0:return x,{'gain':0,'exclusive':exclusive,'candidates':len(b),'edges':int(G.nnz)}
    reverse=np.full(len(a),-1,np.int32)
    for i,j in enumerate(matching):
        if j>=0:reverse[j]=i
    left=set(map(int,np.flatnonzero(matching<0)));right=set();queue=list(left)
    for i in queue:
        for jj in G.indices[G.indptr[i]:G.indptr[i+1]]:
            j=int(jj)
            if matching[i]==j or j in right:continue
            right.add(j)
            if reverse[j]>=0 and int(reverse[j]) not in left:
                left.add(int(reverse[j]));queue.append(int(reverse[j]))
    keep=np.ones(len(x),bool);keep[ai[sorted(right)]]=False
    result=np.unique(np.vstack([x[keep],b[sorted(left)]]),axis=0)
    if len(result)!=len(x)+gain:raise ValueError('exchange count')
    return result,{'gain':gain,'exclusive':exclusive,'candidates':len(b),'edges':int(G.nnz),'removed':len(right),'added':len(left)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--frames',type=int,default=8);ap.add_argument('--cutoff',type=int,default=32);ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    z=np.load(args.source,allow_pickle=False)
    if int(z['rational_denominator'])!=1 or not np.all(z['quadratic_denominators']==1) or not np.all(z['quadratic_radicands']==3):raise ValueError('format')
    A=np.vstack([z['rational_points'],z['quadratic_A']])*2
    B=np.vstack([np.zeros_like(z['rational_points']),z['quadratic_B']])*2
    X=np.hstack([A,B]).astype(np.int16);t0=time.monotonic();records=[]
    for i,blocks in enumerate(frames()[:args.frames]):
        X,r=refill(X,transform(X,blocks),args.cutoff);r.update(frame=i+1,total=len(X),seconds=time.monotonic()-t0);records.append(r)
        print(json.dumps(r),flush=True)
        if len(X)>201180:
            np.savez_compressed(args.out/f'sextet26_{len(X)}.npz',rational=X[:,:26],radical=X[:,26:],radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
        (args.out/'sextet26.json').write_text(json.dumps({'baseline':201180,'current':len(X),'cutoff':args.cutoff,'frames':records,'scope':'Exact finite union search; positive witnesses pending complete geometric replay.'},indent=2)+'\n')


if __name__=='__main__':main()
