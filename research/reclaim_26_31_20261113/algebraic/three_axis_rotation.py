"""Exact three-axis deformations of unsaturated 31D antipodal fibers.

For an old fiber b +/- p, b +/- q and a vacant tail axis r (all three
residual axes have squared norm 16), use the orthogonal matrix

 [[2*n*n, m*m, 2*m*n],
  [m*m, 2*n*n, -2*m*n],
  [-2*m*n, 2*m*n, 2*n*n-m*m]] / (2*n*n+m*m).

If 4*m*n <= 2*n*n+m*m, the projection perpendicular to r of each new
point lies in the convex hull of the four old projected points.  Thus
only those four parents and old points using the vacant tail coordinate
can obstruct the six new points.  The axes p and q may themselves have
both core and tail coordinates.  All decisions use integers.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np

NORMAL = np.array([0,2,4,5,6,10,11,23])
COMPLEMENT = np.array([i for i in range(24) if i not in NORMAL])
PARAMETERS = [(1,16),(1,8),(1,6),(1,4),(1,3),(1,2),(4,7),(7,12),(24,41),(41,70)]


def load_source(path):
    z = np.load(path, allow_pickle=False)
    assert int(z['squared_norm']) == 128
    assert np.all(z['points'] % 2 == 0)
    x = z['points'].astype(np.int16)//2
    assert len(x) == 251254 and x.shape[1] == 31
    assert np.all(np.sum(x.astype(np.int32)**2, axis=1) == 32)
    return x


def fibers(x):
    norm = np.sum(x[:, COMPLEMENT].astype(np.int32)**2, axis=1)
    ids = np.flatnonzero(norm == 16)
    bodies, inv, counts = np.unique(x[ids][:, COMPLEMENT], axis=0,
                                    return_inverse=True, return_counts=True)
    order = np.argsort(inv, kind='stable')
    cumulative = np.r_[0,np.cumsum(counts)]
    for index in np.flatnonzero(counts == 28):
        members = ids[order[cumulative[index]:cumulative[index+1]]]
        y = x[members]
        empty = np.flatnonzero(~np.any(y[:,24:],axis=0))
        equator = np.flatnonzero(~np.any(y[:,24:],axis=1))
        representatives = []
        for local in range(len(y)):
            p = y[local,np.r_[NORMAL,np.arange(24,31)]]
            if p[np.flatnonzero(p)[0]] > 0:
                representatives.append(local)
        assert len(representatives)*2 == len(y)
        # Every residual fiber really is an orthogonal antipodal 14-frame.
        residual = y[:,np.r_[NORMAL,np.arange(24,31)]].astype(np.int32)
        gram = residual @ residual.T
        assert np.all(np.isin(gram,[-16,0,16]))
        assert np.sum(gram == -16) == 28
        assert np.sum(gram == 16) == 28
        # Sixteen other rank-14 fibers have a four-support tail null vector;
        # this entry point intentionally treats the coordinate-axis fibers.
        if len(empty) != 1:
            continue
        yield dict(index=int(index),body=bodies[index],members=members,
                   equator=equator,representatives=representatives,
                   tail_axis=int(empty[0])+24)


def construct(b,p,q,axis,m,n):
    den=2*n*n+m*m
    assert m > 0 and 4*m*n <= den
    r=np.zeros(31,dtype=np.int64);r[axis]=4
    axes=np.array([2*n*n*p+m*m*q-2*m*n*r,
                   m*m*p+2*n*n*q+2*m*n*r,
                   2*m*n*p-2*m*n*q+(2*n*n-m*m)*r])
    assert np.array_equal(axes@axes.T,16*den*den*np.eye(3,dtype=np.int64))
    added=np.array([den*b+sign*v for v in axes for sign in [-1,1]])
    assert np.all(np.sum(added*added,axis=1)==32*den*den)
    return added,den


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--max-fibers',type=int,default=4004)
    ap.add_argument('--pairs',type=int,default=3,help='Maximum residual axis pairs per fiber; 0 means all.')
    ap.add_argument('--save-cutoff',type=int,default=10,help='Save bundles with at most this many exterior blockers.')
    ap.add_argument('--seed',type=int,default=20261113)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();x=load_source(args.source)
    source_hash=hashlib.sha256(args.source.read_bytes()).hexdigest()
    fs=list(fibers(x));available=len(fs);rng=np.random.default_rng(args.seed)
    rng.shuffle(fs);fs=fs[:args.max_fibers]
    tailids={j:np.flatnonzero(x[:,j]) for j in range(24,31)}
    float_x={j:x[ids].astype(np.float32) for j,ids in tailids.items()}
    records=[];saved=[];hist=Counter();best=10**9;count=0;unique_states=set()
    for fi,f in enumerate(fs):
        axis=f['tail_axis'];oldids=tailids[axis];yf=float_x[axis]
        reps=f['representatives'];pairs=list(itertools.combinations(reps,2))
        if args.pairs:
            pure=set(f['equator'].tolist())
            groups=[[pair for pair in pairs if sum(local in pure for local in pair)==k] for k in [2,1,0]]
            for group in groups:rng.shuffle(group)
            pairs=[]
            while len(pairs)<args.pairs and any(groups):
                for group in groups:
                    if group and len(pairs)<args.pairs:pairs.append(group.pop())
        b=np.zeros(31,dtype=np.int64);b[COMPLEMENT]=f['body']
        beta=(yf@b.astype(np.float32)).astype(np.int32)
        t=4*x[oldids,axis].astype(np.int32)
        for pair in pairs:
            p=x[f['members'][pair[0]]].astype(np.int64)-b
            q=x[f['members'][pair[1]]].astype(np.int64)-b
            assert int(p@p)==int(q@q)==16 and int(p@q)==0
            key=(tuple(b),tuple(p),tuple(q),axis)
            assert key not in unique_states;unique_states.add(key)
            parents=[]
            for local in range(len(f['members'])):
                row=x[f['members'][local]]
                if np.array_equal(row,b+p) or np.array_equal(row,b-p) or np.array_equal(row,b+q) or np.array_equal(row,b-q):
                    parents.append(int(f['members'][local]))
            assert len(parents)==4
            P=(yf@p.astype(np.float32)).astype(np.int32)
            Q=(yf@q.astype(np.float32)).astype(np.int32)
            c=16-beta
            assert np.all(np.abs(P)<=c) and np.all(np.abs(Q)<=c)
            pairbest=None
            for signq,signr in itertools.product([-1,1],repeat=2):
                PP=P;QQ=signq*Q;TT=signr*t
                S=PP+QQ;D=PP-QQ
                # An exact Cauchy envelope removes points safe for every parameter.
                K=D.astype(np.int64)**2+2*TT.astype(np.int64)**2
                critical=(K>(2*c-np.abs(S)).astype(np.int64)**2)|(K>2*c.astype(np.int64)**2)
                pp,qq,tt,cc=PP[critical].astype(np.int64),QQ[critical].astype(np.int64),TT[critical].astype(np.int64),c[critical].astype(np.int64)
                critical_ids=oldids[critical]
                for m,n in PARAMETERS:
                    den=2*n*n+m*m
                    a1=2*n*n*pp+m*m*qq-2*m*n*tt
                    a2=m*m*pp+2*n*n*qq+2*m*n*tt
                    a3=2*m*n*pp-2*m*n*qq+(2*n*n-m*m)*tt
                    blocked=(np.maximum.reduce([np.abs(a1),np.abs(a2),np.abs(a3)])>den*cc)
                    blockers=critical_ids[blocked]
                    cost=len(blockers);hist[cost]+=1;count+=1
                    rec=dict(fiber=f['index'],axis=axis,parents=parents,
                             parent_axis_tail_norms=[int(p[24:]@p[24:]),int(q[24:]@q[24:])],parameter=[m,n],
                             sign_q=signq,sign_r=signr,exterior_blockers=cost,
                             net_gain=2-cost)
                    if pairbest is None or cost<pairbest['exterior_blockers']:pairbest=rec
                    if cost<best:
                        best=cost
                        print(json.dumps(dict(event='best',fibers_completed=fi+1,**rec)),flush=True)
                    if cost<=args.save_cutoff:
                        add,den=construct(b,p,signq*q,axis,m,n)
                        if signr==-1:
                            add[:,axis]*=-1
                        # Independently replay every admitted exterior blocker.
                        dots=add.astype(np.int64)@x[oldids].T.astype(np.int64)
                        direct=oldids[np.any(dots>16*den,axis=0)]
                        assert np.array_equal(np.sort(direct),np.sort(blockers))
                        saved.append(dict(**rec,denominator=den,added_numerators=add.tolist(),
                                          blockers=blockers.tolist()))
            records.append(pairbest)
        if (fi+1)%100==0:
            print(json.dumps(dict(event='progress',fibers_completed=fi+1,states=count,
                                 best_exterior_cost=best,saved=len(saved),seconds=time.monotonic()-started)),flush=True)
    report=dict(source=str(args.source),source_sha256=source_hash,baseline=len(x),dimension=31,
                fibers_available=available,fibers_tested=len(fs),axis_pairs_tested=len(records),
                parameter_states=count,parameters=[list(p) for p in PARAMETERS],
                best_exterior_cost=best,best_net_gain=2-best,
                exterior_cost_histogram={str(k):v for k,v in sorted(hist.items())},
                saved_bundles=len(saved),records=records,seconds=time.monotonic()-started,
                scope='Finite exact three-axis rotations of old size28 fibers. Source validity is inherited; no global optimality is claimed.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (args.out/'bundles.json').write_text(json.dumps(saved,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['records','exterior_cost_histogram']},indent=2),flush=True)


if __name__=='__main__':main()
