"""Exact obstruction for one 13D projection of a 14D E7-frame configuration."""
import pathlib,numpy as np
z=np.load(pathlib.Path(__file__).with_name('certificate.npz'));X=z['source'].astype(np.int64);M=z['matching'].astype(np.int64);u=np.array([1,-1,-1]+[0]*11,np.int64);assert X.shape==(1932,14) and len(set(map(tuple,X.tolist())))==len(X) and np.all((X*X).sum(1)==8)
for st in range(0,len(X),96):
 D=X[st:st+96]@X.T
 for k in range(len(D)):D[k,:st+k+1]=-99
 assert D.max(initial=-99)<=4
Y=3*X-(X@u)[:,None]*u;N=(Y*Y).sum(1);assert np.all(N>0) and not np.any(Y@u)
equ=np.flatnonzero(X@u==0);assert len(equ)==1122
for st in range(0,len(equ),96):
 ids=equ[st:st+96];D=Y[ids]@Y[equ].T;B=(D>0)&(4*D*D>N[ids,None]*N[equ][None,:])
 for k in range(len(B)):B[k,:st+k+1]=False
 assert not np.any(B)
assert M.shape==(810,2) and len(set(map(int,M.ravel())))==1620 and M.min()>=0 and M.max()<len(X)
for i,j in M:
 d=int(Y[i]@Y[j]);assert d>0 and 4*d*d>int(N[i]*N[j])
print('PASS: 810 pairwise disjoint conflicts give at most 1932-810=1122; the equator attains 1122.')
