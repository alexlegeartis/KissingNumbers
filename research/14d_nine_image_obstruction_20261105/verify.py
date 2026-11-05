"""Verify a restricted 14D orbit-union obstruction, using exact integer arithmetic.

Run from repository root: python research/14d_nine_image_obstruction_20261105/verify.py
Requires numpy and scipy. The certificate makes no assertion outside this union.
"""
from collections import deque
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.sparse import csr_matrix,eye
root=Path(__file__).resolve().parents[2]; here=Path(__file__).resolve().parent
raw=(here/'baseline.bin').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='25fd242c2c5b52afb8f065fc4651d466d19f7b4016809eb847e17f1df5099e6b'
X=np.frombuffer(raw,dtype=np.int8).reshape(-1,14).astype(np.int64); assert len(X)==1932 and len({tuple(v) for v in X})==1932
assert np.all(np.sum(X*X,axis=1)==8)
D=X@X.T;np.fill_diagonal(D,0);assert D.max()<=4
recipe=json.loads((root/'discovery/recipes/group14_strict504.json').read_text()); gens=recipe['nodes'][0]['params']['generators']
assert len(gens)==2 and all(g['denominator']==2 for g in gens); generators=[np.array(g['matrix'],dtype=np.int64) for g in gens]
I=2*np.eye(14,dtype=np.int64)
for g in generators:assert np.array_equal(g@g.T,4*np.eye(14,dtype=np.int64))
queue=deque([I]); group=[I]; seen={tuple(I.ravel())}
while queue:
 g=queue.popleft()
 for h in generators:
  z=g@h;assert np.all(z%2==0);z=z//2;k=tuple(z.ravel())
  if k not in seen:
   assert len(group)<504;seen.add(k);group.append(z);queue.append(z)
assert len(group)==504
images={frozenset(map(tuple,X@g.T)) for g in group};assert len(images)==9 and all(len(s)==1932 for s in images)
core=set.intersection(*(set(s) for s in images));union=set.union(*(set(s) for s in images));assert len(core)==1764 and len(union)==3276
C=np.array(sorted(core),dtype=np.int64);Y=np.array(sorted(union-core),dtype=np.int64);assert len(Y)==1512
assert np.all(np.sum(Y*Y,axis=1)==32) and (Y@C.T).max()<=16
Dcore=C@C.T;np.fill_diagonal(Dcore,0);assert Dcore.max()<=16
D=Y@Y.T;A=csr_matrix((D>16).astype(np.int64));A.setdiag(0);A.eliminate_zeros();assert (A!=A.T).nnz==0 and np.all(np.asarray(A.sum(1)).ravel()==32)
# Exact polynomial certificate. Every intermediate is checked before the next
# multiplication; each entry of A is 0/1 and each row has 32 ones, so the
# following bound rules out int64 overflow at the next operation.
Id=eye(len(Y),dtype=np.int64,format='csr');P=Id.copy()
for r in [32,24,16,12,6,4,2,0,-4]:
 m=int(np.max(np.abs(P.data),initial=0));assert (32+abs(r))*m<2**62
 P=A@P-r*P;P.eliminate_zeros()
assert P.nnz==0
print('Group order: 504; 9 images; core: 1764; variable vertices: 1512.')
print('Exact p(A)=0 for roots -4,0,2,4,6,12,16,24,32; graph is symmetric and 32-regular.')
print('Hoffman bound: alpha <= 1512*4/(32+4) = 168. The entire union contains at most 1932 compatible points.')

# Extend the same calculation to an outside orbit with an order-nine stabilizer.
# It lies in the same 14D space but has a different norm, so all cross tests
# compare squared positive dot products rather than rounded normalizations.
seed=np.array([0,40,0,40,0,-40,0,0,9,0,9,0,-9,0],dtype=np.int64)
Z=np.unique(np.array([g@seed for g in group]),axis=0);assert len(Z)==56
nz=int(Z[0]@Z[0]);assert np.all(np.sum(Z*Z,axis=1)==nz)
d=Z@Z.T;np.fill_diagonal(d,0);assert np.all(4*d[d>0]**2<=nz*nz)
d=Z@C.T;assert np.all(4*d[d>0]**2<=nz*32)
assert not np.any((d>0)&(d*d==nz*32))
d=Z@Y.T;assert not np.any((d>0)&(d*d==nz*32))
B=(d>0)&(4*d*d>nz*32);assert np.all(B.sum(axis=1)==27)
counts=B.sum(axis=0);assert np.count_nonzero(counts==2)==756 and np.count_nonzero(counts==0)==756
hit=counts==2;other=~hit
assert (A[hit][:,other]).nnz==0 and np.all(np.asarray(A[other][:,other].sum(axis=1)).ravel()==32)
for mask in B:
 H=A[mask][:,mask].toarray();assert H.shape==(27,27) and np.all(H.sum(axis=1)==16)
 assert np.array_equal(H@H,8*np.eye(27,dtype=np.int64)+2*H+8*np.ones((27,27),dtype=np.int64))
# H has spectrum contained in {16,4,-2}, so Hoffman gives alpha(H)<=3.
# The 56 blocks cover the hit vertices twice. If k external points and m hit
# vertices are chosen, summing block inequalities gives 2m+3k<=168; hence
# m+k<=84. The other 756-vertex 32-regular component contributes at most 84,
# by the same -4 eigenvalue bound as above.
print('Extended certificate: 56 extra orbit points; 56 Schlaefli blocks cover one 756-vertex component twice.')
print('For m selected old vertices and k extra points, 2m+3k<=168. The extended union still has maximum 1932.')

# Broader candidate universe: allow arbitrary re-pairings, not just 9 images.
left=np.unique(C[:,:7],axis=0);right=np.unique(C[:,7:],axis=0)
assert left.shape==(126,7) and right.shape==(126,7)
assert np.all(np.sum(left*left,axis=1)==16) and np.all(np.sum(right*right,axis=1)==16)
assert all((np.all(v[:7]==0) and v[7:]@v[7:]==32) or
           (np.all(v[7:]==0) and v[:7]@v[:7]==32) for v in Y)
# For each fixed left root u, compatible product points (u,v) require
# v dot v' <= 0. The orthoplex/Rankin bound in R^7 gives at most 14 such v.
# Thus all 126*126 possible product points contribute at most 126*14=1764.
# The independent axis pool contributes at most 168, even when the 56 extra
# orbit points are added, by the preceding certificate. Cross-conflicts only
# lower the total. The original baseline shows equality is attained.
print('Broader scope: all 126*126 root-product points, 1512 axis points, and the 56 outside-orbit points have optimum 1932.')
