from pathlib import Path
import json,hashlib
from fractions import Fraction
import numpy as np
P=Path(__file__).resolve().parent;S=np.load(P/'source_supports.npz')
# F4=F2[x]/(x^2+x+1), represented0,1,x,x+1.
def mul(a,b):
 r=0
 while b:
  if b&1:r^=a
  a<<=1
  if a&4:a^=7
  b>>=1
 return r
M=np.array([[mul(a,b) for b in range(4)] for a in range(4)],dtype=np.uint8)
Q=S['quaternary_generator'].copy();assert Q.shape==(9,24)
pivot=int(np.flatnonzero(Q[:,0])[0]);Q[[0,pivot]]=Q[[pivot,0]];Q[0]=M[[0,1,3,2][Q[0,0]],Q[0]]
for i in range(1,9):Q[i]^=M[Q[i,0],Q[0]]
assert Q[0,0]==1 and not Q[1:,0].any();G4=Q[1:,2:]
def span(G,q):
 V=np.zeros((1,G.shape[1]),np.uint8)
 for row in G:V=np.concatenate([V^M[a,row] for a in range(q)])
 return V
C=span(G4,4);assert C.shape==(65536,22) and len(np.unique(C,axis=0))==65536
w4=np.bincount(np.count_nonzero(C,axis=1),minlength=23);assert w4[0]==1 and not w4[1:11].any()
B=S['sign_generator'].copy();assert B.shape==(14,24)
for j in range(2):
 pivot=j+int(np.flatnonzero(B[j:,j])[0]);B[[j,pivot]]=B[[pivot,j]]
 for i in range(14):
  if i!=j and B[i,j]:B[i]^=B[j]
assert np.array_equal(B[:2,:2],np.eye(2,dtype=np.uint8)) and not B[2:,:2].any();G2=B[2:,2:]
D=span(G2,2);assert D.shape==(4096,22) and len(np.unique(D,axis=0))==4096
w2=np.bincount(np.count_nonzero(D,axis=1),minlength=23);assert w2[0]==1 and not w2[1:6].any()
# All counts and squared inequalities below exact.
a2=Fraction(48,11);assert 22*a2==96 and 22**2*a2<48**2 and 11*a2==48 and 4*a2*6>=96
count=65536*4096;assert count==268435456
np.save(P/'quaternary22_generator.npy',G4);np.save(P/'binary22_generator.npy',G2)
base={92:25769803776,93:38654705664,94:55834574848,95:77309411328}
obj=dict(status='PASS exact finite code enumeration and sparse geometry',quaternary_weight_histogram=w4.tolist(),binary_weight_histogram=w2.tolist(),supports=65536,signs=4096,sparse_count=count,amplitude_squared=str(a2),totals={d:n+count for d,n in base.items()},source_sha256=hashlib.sha256((P/'source_supports.npz').read_bytes()).hexdigest())
(P/'sparse_certificate.json').write_text(json.dumps(obj,indent=2)+'\n');print(json.dumps(obj,indent=2))
