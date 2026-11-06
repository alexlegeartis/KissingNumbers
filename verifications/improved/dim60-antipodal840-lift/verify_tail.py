"""Reconstruct and exactly verify Takhanov's antipodal 840-point tail."""
from pathlib import Path
import itertools,json,hashlib
import numpy as np
if not __debug__:raise SystemExit('Run without -O: exact verification uses assertions.')
D=Path(__file__).parent

def mul(a,b):
 return np.r_[a[0]*b[0]-a[1:]@b[1:],a[0]*b[1:]+b[0]*a[1:]+np.cross(a[1:],b[1:])].astype(np.int64)
C=[];F=[]
for p in np.array([[2,0,0,0],[1,1,1,1],[-1,1,1,1]],np.int64):
 c=np.array([mul(p,e)//2 for e in 2*np.eye(4,dtype=np.int64)])
 f=np.array([mul(np.array([1,1,0,0]),q)//2 for q in c]);C.append(c);F.append(f)
C=np.array(C);F=np.array(F)
assert all(np.array_equal(c@c.T,4*np.eye(4,dtype=int)) for c in C)
assert all(np.array_equal(f@f.T,2*np.eye(4,dtype=int)) for f in F)
rows=[]
# Coordinates are (A+B sqrt(2))/4.
for block,r,j,s in itertools.product(range(3),range(3),range(4),[-1,1]):
 A=np.zeros(12,np.int64);B=A.copy();A[4*block:4*block+4]=s*2*C[r,j];rows.append((A,B))
cells=[]
for r in range(3):cells.extend([[0,3+r,r,r,0],[1,r,3+r,3+r,[2,1,3][r]]])
for h,f,g,k,delta in cells:
 for a,b in itertools.product(range(4),repeat=2):
  labels=[a,b,a^b^delta]
  for signs in itertools.product([-1,1],repeat=3):
   A=np.zeros(12,np.int64);B=A.copy()
   for j,(frame,label,s) in enumerate(zip([f,g,k],labels,signs)):
    if j==h:
     assert frame>=3;A[4*j:4*j+4]=2*s*F[frame-3,label]
    elif frame<3:A[4*j:4*j+4]=s*C[frame,label]
    else:B[4*j:4*j+4]=s*F[frame-3,label]
   rows.append((A,B))
A=np.array([a for a,b in rows]);B=np.array([b for a,b in rows]);assert A.shape==B.shape==(840,12)
keys={tuple(np.r_[a,b]) for a,b in rows};assert len(keys)==840 and all(tuple(-np.r_[a,b]) in keys for a,b in rows)
assert np.all((A*A+2*B*B).sum(1)==16) and np.all((2*A*B).sum(1)==0)
U=A@A.T+2*(B@B.T);V=A@B.T+B@A.T

def sign(a,b):
 a=int(a);b=int(b)
 if b==0:return (a>0)-(a<0)
 if a==0:return (b>0)-(b<0)
 if a>0 and b>0:return 1
 if a<0 and b<0:return -1
 t=a*a-2*b*b
 return ((t>0)-(t<0))*(1 if a>0 else -1)
contacts=0
for i in range(840):
 for j in range(i):
  s=sign(U[i,j]-8,V[i,j]);assert s<=0,(i,j,int(U[i,j]),int(V[i,j]));contacts+=s==0
np.savez_compressed(D/'tail840_exact.npz',A=A,B=B)
report={'points':840,'dimension':12,'antipodal_pairs':420,'pairs_checked':840*839//2,'contacts':contacts,'violations':0,'field':'(A+B sqrt(2))/4','cells':cells,'coordinate_sha256':hashlib.sha256(A.tobytes()+B.tobytes()).hexdigest()}
(D/'tail_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
