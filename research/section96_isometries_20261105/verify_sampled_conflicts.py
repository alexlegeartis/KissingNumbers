"""Replay exact sampled edges and a matching; not a full graph proof."""
import json,pathlib
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
z=json.loads(pathlib.Path(__file__).with_name('sample92_conflicts.json').read_text());m=z['m'];mask=z['mask'];M=z['mapping'];assert m==4 and mask==1008 and len(M)==24
def gray(q):return sum(((((q>>2)&1)^((q&1)*x)^(((q>>1)&1)*y))<<(2*x+y)) for x in (0,1) for y in (0,1))
for row in M:
 assert len(row)==8 and len(set(row))==8 and all(type(v) is int and 0<=v<16 for v in row)
 for a in range(8):
  for b in range(8):assert (row[a]^row[b]).bit_count()==(gray(a)^gray(b)).bit_count()
for row in M[:m]:assert all(x.bit_count()%2==0 for x in row)
g=[7,1,4,4,2,7,3,3,4,7,2,1];G=[[0]*24 for _ in range(12)]
for i in range(12):
 for j,x in enumerate(g):G[i][i+j]=x
 G[i][23]=3
L={};R={};edges=set()
for astr,bstr,declared in z['witnesses']:
 assert len(astr)==len(bstr)==12 and set(astr+bstr)<=set('01234567');a=list(map(int,astr));b=list(map(int,bstr));assert all(a[j]%2==0 for j in range(12) if mask>>j&1)
 x=[sum(a[i]*G[i][j] for i in range(12))%8 for j in range(24)];y=[sum(b[i]*G[i][j] for i in range(12))%8 for j in range(24)]
 assert all(M[j][x[j]].bit_count()==2 and gray(y[j]).bit_count()==2 for j in range(m))
 d=sum((M[j][x[j]]^gray(y[j])).bit_count() for j in range(24));assert d==declared and d<24
 u=tuple(q%4 for q in a);v=tuple(q%4 for q in b);edges.add((L.setdefault(u,len(L)),R.setdefault(v,len(R))))
i,j=zip(*edges);A=coo_matrix((np.ones(len(edges),np.uint8),(i,j)),shape=(len(L),len(R))).tocsr();size=int((maximum_bipartite_matching(A,perm_type='column')>=0).sum());assert len(L)==497 and len(R)==499 and size==497
print('verified witnesses',len(z['witnesses']),'distinct edges',len(edges),'matching',size)
