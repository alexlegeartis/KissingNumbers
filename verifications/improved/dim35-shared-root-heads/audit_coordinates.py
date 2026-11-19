#!/usr/bin/env python3
"""Independent coordinate expansion check of all interactions with changed points."""
from pathlib import Path
import itertools,json
import numpy as np
if not __debug__:
    raise SystemExit("Run without -O: verification assertions must be enabled.")
P=Path(__file__).resolve().parent
src=json.loads((P/'code_data.json').read_text())
parities=json.loads((P/'full_parities.json').read_text())
assert {m for m,p in parities}=={m for m,p in src['full_parity_supports']} and len(parities)==2202
z=np.load(P/'shared_positive.npz',allow_pickle=False)
stars=z['new_stars'].astype(np.int64); heads=z['heads'].astype(np.int64);D=z['dense_positions'].astype(np.int64)
span=np.array([src['dense_coset_shift']],dtype=np.uint32)
for row in src['dense_generator']:span=np.r_[span,span^np.uint32(row)]
dense=np.zeros((len(span),35),np.int8)
dense[:,D]=(1-2*((span[:,None]>>np.arange(32))&1)).astype(np.int8)
signs=np.array(list(itertools.product([2,-2],repeat=8)),dtype=np.int8); odd=(signs<0).sum(1)%2
middle=np.zeros((128*len(parities),35),np.int8)
for j,(support,par) in enumerate(parities):
 S=[i for i in range(35) if support>>i&1]; assert len(S)==8 and par in (0,1)
 middle[j*128:(j+1)*128,S]=signs[odd==par]
cut={tuple(r)for r in z['removed_roots']}
roots=[]
for i,j in itertools.combinations(range(35),2):
 for p,q in itertools.product([4,-4],repeat=2):
  x=np.zeros(35,np.int8);x[[i,j]]=[p,q]
  if tuple(x) not in cut:roots.append(x)
anchors=np.vstack([dense,middle,roots]).astype(np.int8)
N=np.vstack([stars,heads]);nz=np.sum(N*N,axis=1); assert len(N)==162
maxima=np.zeros(len(N),dtype=np.int64)
# Every partial sum is far below the float32 exact integer range.
for u in range(0,len(anchors),2048):
 A=anchors[u:u+2048].astype(np.float32)
 assert np.all(np.sum(A*A,axis=1)==32)
 m=(N.astype(np.float32)@A.T).astype(np.int64).max(1)
 maxima=np.maximum(maxima,m)
assert np.all(4*maxima**2<=32*nz)
H=N@N.T;np.fill_diagonal(H,0);assert np.all(4*np.maximum(H,0)**2<=nz[:,None]*nz)
assert len(anchors)+len(N)==415336
print('PASS 415336 points; 67258188 changed/anchor pair products verified' if len(N)*len(anchors)==67258188 else f'PASS 415336 points; {len(N)*len(anchors)} changed/anchor pair products verified')
