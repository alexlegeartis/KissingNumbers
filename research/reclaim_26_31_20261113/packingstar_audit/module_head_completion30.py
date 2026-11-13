"""Price every integer norm16+norm16 head in the distinguished 30D module.

The body catalogue contains all 2160 norm16 E8 deep-hole vectors in the
normal coordinate model. The 252 tails have shapes (4,0^5) or (2^4,0^2).
Unchanged source validity remains an inherited premise. This scans a finite
catalogue, not arbitrary real heads, radii, or global configurations.
"""
import argparse
from collections import Counter
import hashlib
import json
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

NORMAL = np.array([0,2,4,5,6,10,11,23])


def catalogues():
    body, tail = [], []
    for dim, out in ((8, body), (6, tail)):
        for i in range(dim):
            for sign in (-4,4):
                row = np.zeros(dim, np.int16)
                row[i] = sign
                out.append(row)
        for support in combinations(range(dim),4):
            for signs in product((-2,2),repeat=4):
                row = np.zeros(dim,np.int16)
                row[list(support)] = signs
                out.append(row)
    for i in range(8):
        for signs in product((-1,1),repeat=8):
            if sum(s < 0 for s in signs) % 2 == 0:
                continue
            row = np.array(signs,np.int16)
            row[i] *= 3
            body.append(row)
    body, tail = np.array(body), np.array(tail)
    assert body.shape == (2160,8) and tail.shape == (252,6)
    assert np.all(np.sum(body.astype(np.int64)**2,axis=1)==16)
    assert np.all(np.sum(tail.astype(np.int64)**2,axis=1)==16)
    return body, tail


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('source',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--cutoff',type=int,default=4)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    z=np.load(args.source)
    X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:
        assert not np.any(X%2)
        X//=2
    assert X.shape[1]==30 and np.all(np.sum(X.astype(np.int64)**2,axis=1)==32)
    body,tail=catalogues()
    q=np.sum(X[:,24:].astype(np.int64)**2,axis=1)
    h=np.sum(X[:,NORMAL].astype(np.int64)**2,axis=1)
    relevant=((q==8)&(h>8))|((q==16)&(h>0))
    ids=np.flatnonzero(relevant)
    projection=X[:,np.r_[NORMAL,np.arange(24,30)]]
    # Explicit exact-integer guards cover every point omitted from the graph.
    assert np.max(body.astype(np.float32)@np.unique(projection[q==0,:8],axis=0).astype(np.float32).T)<=16
    assert np.max(body.astype(np.float32)@np.unique(projection[(q==8)&~relevant,:8],axis=0).astype(np.float32).T)<=8
    assert np.max(tail.astype(np.float32)@np.unique(projection[q==8,8:],axis=0).astype(np.float32).T)<=8
    assert np.max(tail.astype(np.float32)@np.unique(projection[q==32,8:],axis=0).astype(np.float32).T)<=16
    oldset={row.tobytes() for row in projection}
    core_products=body.astype(np.float32)@projection[ids,:8].astype(np.float32).T
    points=[]
    blockers=[]
    hist=Counter()
    for t in tail:
        mask=core_products+(projection[ids,8:].astype(np.int32)@t.astype(np.int32))[None,:]>16
        costs=mask.sum(axis=1)
        hist.update(map(int,costs))
        for i in np.flatnonzero(costs<=args.cutoff):
            row=np.r_[body[i],t]
            if row.tobytes() in oldset:
                continue
            points.append(row)
            blockers.append(ids[np.flatnonzero(mask[i])])
    C=np.array(points,np.int16).reshape((-1,14))
    offsets=np.r_[0,np.cumsum([len(v) for v in blockers])]
    old=np.concatenate(blockers) if len(C) else np.array([],np.int64)
    matrix=csr_matrix((np.ones(len(old),np.int8),old,offsets),shape=(len(C),len(X)))
    matching=maximum_bipartite_matching(matrix,perm_type='column')
    free=[i for i,b in enumerate(blockers) if len(b)==0]
    report=dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),source_count=len(X),
                candidate_universe=len(body)*len(tail),relevant_old=len(ids),cutoff=args.cutoff,
                admitted_new=len(C),conflict_incidences=len(old),matched_new=int(np.sum(matching>=0)),
                free=len(free),price_histogram=dict(sorted(hist.items())),
                scope='Full integer2160x252 head catalogue priced; selection limited by cutoff.')
    np.savez_compressed(args.out/'pool.npz',candidates=C,offsets=offsets,blockers=old,matching=matching)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    if free:
        F=C[free]
        for row in F:
            a=np.zeros(30,np.int16)
            a[NORMAL],a[24:]=row[:8],row[8:]
            dots=sum((X[:,j].astype(np.int32)*int(a[j]) for j in np.flatnonzero(a)),start=np.zeros(len(X),np.int32))
            assert dots.max()<=16
        np.savez_compressed(args.out/'free_candidates.npz',heads=F)


if __name__=='__main__':
    main()
