#!/usr/bin/env python3
"""Rebuild and verify the 232,236-point witness from the public 232,230 stock."""
import argparse
import hashlib
import json
from collections import deque
from pathlib import Path
import numpy as np

SOURCE_SHA256 = "66a515c03aff33122f4eb06920a9b7f6d13728aefb4e3bead12be80c7cf91c71"
SOURCE_POINTS_SHA256 = "3675c465c219c0e542f881356774f85dc9d67cb7ad0f83cfb805d0e808138233"
TARGET = (4,-2,0,-2,0,0,4,0,0,0,0,4,0,0,2,-2,2,2,-2,0,0,2,0,4,4,0,0,0,4)
A = np.array((0,4,0,-4,0),np.int32)
B = np.array((4,0,0,0,-4),np.int32)
EXTRA_BODIES = (
    ((-4, 0, 0, -2, 0, 0, 4, 0, 0, 0, 0, 4, 2, -2, 0, 0, 0, 2, 2, 2, 0, -2, -2, -4), 'B'),
    ((-4, 0, 0, 2, 0, 0, 4, 0, 0, 0, 0, 4, -2, 2, 0, 0, 0, -2, -2, -2, 0, 2, 2, -4), 'B'),
    ((-2, 0, 2, 0, 2, -6, -2, 2, 2, 2, 2, -2, -2, -2, 0, 0, 0, 0, 0, -2, 2, 0, -2, -2), 'A'),
    ((0, 0, -4, 0, -4, 4, 0, 2, -2, 2, 4, 0, 2, -2, 0, 0, 0, 0, 0, -2, -2, 0, 2, 0), 'B'),
    ((2, 0, 6, 0, -2, -2, 2, 2, 2, 2, -2, 2, -2, -2, 0, 0, 0, 0, 0, -2, 2, 0, -2, 2), 'A'),
)

def point_hash(x):
    order=np.lexsort(x.T[::-1])
    return hashlib.sha256(x[order].astype("<i2").tobytes()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("source",type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--full",action="store_true",help="also check every retained-retained pair")
    opt=p.parse_args()
    assert hashlib.sha256(opt.source.read_bytes()).hexdigest()==SOURCE_SHA256
    with np.load(opt.source,allow_pickle=False) as f:
        source=np.asarray(f["points"],dtype=np.int32)
    assert source.shape==(232230,29)
    assert point_hash(source)==SOURCE_POINTS_SHA256
    H,T=source[:,:24],source[:,24:]
    target=np.flatnonzero((source==TARGET).all(axis=1))
    assert len(target)==1
    h=H[target[0]]
    assert (H==h).all(axis=1).sum()==1
    frame=[]
    for t in (A,B):
        sides=[]
        for z in (t,-t):
            sides.append({tuple(H[i]):int(i) for i in np.flatnonzero((T==z).all(axis=1))})
        frame.append([(x,sides[0][x],sides[1][x]) for x in sorted(sides[0].keys()&sides[1].keys())])
    X,Y=[np.array([r[0] for r in group],np.int32) for group in frame]
    assert not (X@h+abs(A@B)>64).any()
    seed=np.flatnonzero(Y@h+int(B@B)>64)
    assert len(seed)>0
    graph=X@Y.T>32
    J=(set(),set(map(int,seed)))
    pending=deque((1,int(j)) for j in seed)
    while pending:
        side,j=pending.popleft()
        near=np.flatnonzero(graph[:,j] if side else graph[j])
        for k in near:
            k=int(k)
            if k not in J[1-side]:
                J[1-side].add(k)
                pending.append((1-side,k))
    assert tuple(map(len,J))==(1896,1896)
    keep=np.ones(len(source),bool);keep[target]=False
    new=[]
    for side,ids in enumerate(J):
        t=B if side==0 else A
        for j in sorted(ids):
            body,i,k=frame[side][j]
            keep[[i,k]]=False
            new.extend((np.r_[body,t],np.r_[body,-t]))
    new.extend((np.r_[h,B],np.r_[h,-B]))
    for body,colour in EXTRA_BODIES:
        locus=np.flatnonzero((H==body).all(axis=1))
        assert len(locus)==1 and keep[locus[0]]
        keep[locus]=False
        t=A if colour=="A" else B
        new.extend((np.r_[body,t],np.r_[body,-t]))
    added=np.array(new,np.int32)
    old=source[keep]
    assert len(old)==224640 and len(added)==7596
    assert (np.einsum("ij,ij->i",added,added)==128).all()
    # Each partial sum of each dot product is an integer of magnitude at most
    # 128 by Cauchy-Schwarz. Binary32 matmul is therefore exact on these inputs.
    oldT=old.astype(np.float32).T.copy()
    anew=added.astype(np.float32)
    newT=anew.T.copy()
    max_old,max_new=-128,-128
    for start in range(0,len(added),256):
        z=anew[start:start+256]
        max_old=max(max_old,int(np.max(z@oldT)))
        g=z@newT
        r=np.arange(len(z));g[r,r+start]=-128
        max_new=max(max_new,int(np.max(g)))
    assert max_old<=64 and max_new<=64
    if opt.full:
        assert (np.einsum("ij,ij->i",old,old)==128).all()
        for start in range(0,len(old),256):
            z=oldT[:,start:start+256].T.copy()
            g=z@oldT[:,start:]
            r=np.arange(len(z));g[r,r]=-128
            assert g.max()<=64
    final=np.vstack((old,added))
    dg=point_hash(final)
    print(json.dumps(dict(points=len(final),old_removed=int((~keep).sum()),new_rows=len(added),
                          switched_pair_bodies=sum(map(len,J)),new_vs_retained_max=max_old,
                          new_vs_new_max=max_new,sorted_int16_rows_sha256=dg),indent=2))
    if opt.output:
        np.savez_compressed(opt.output,points=final.astype(np.int16),
                            squared_norm=np.int64(128),inner_product_limit=np.int64(64))

if __name__=="__main__":main()
