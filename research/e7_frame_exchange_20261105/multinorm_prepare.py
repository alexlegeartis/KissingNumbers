"""Enumerate two unequal-norm shells around the checked E7 seed."""
from pathlib import Path
import numpy as np,itertools,collections,time
p=Path(__file__).resolve().parent/'cache';z=np.load(p/'candidates14.npz');old=z['old'].astype(np.int16);X=np.r_[old,-old];start=time.time()
def collect3():
 rows=[];blocks=[];hist=collections.Counter()
 for S in itertools.combinations(range(14),3):
  for signs in itertools.product((-1,1),repeat=3):
   q=np.zeros(14,np.int16);q[list(S)]=signs;bad=np.flatnonzero(X@q>2);hist[len(bad)]+=1
   if len(bad)<=9:rows.append(q);blocks.append(bad)
 return rows,blocks,hist
def collect18():
 rows=[];blocks=[];hist=collections.Counter()
 for i in range(14):
  for j in range(14):
   if i==j:continue
   for k,l in itertools.combinations([x for x in range(14) if x not in (i,j)],2):
    C=np.zeros((16,14),np.int16)
    for s in range(16):C[s,[i,j,k,l]]=[3*(1-2*(s&1)),1*(1-2*((s>>1)&1)),2*(1-2*((s>>2)&1)),2*(1-2*((s>>3)&1))]
    D=C@X.T
    for q,d in zip(C,D):
     bad=np.flatnonzero(d>6);hist[len(bad)]+=1
     if len(bad)<=9:rows.append(q.copy());blocks.append(bad)
 return rows,blocks,hist
r3,b3,h3=collect3();r18,b18,h18=collect18();assert h3[3]==112 and h18[4]==3936 and h18[9]==3456
Q=np.array(r3+r18,np.int8);norm=np.r_[np.full(len(r3),3,np.int16),np.full(len(r18),18,np.int16)];B=np.full((len(Q),9),65535,np.uint16)
for i,b in enumerate(b3+b18):B[i,:len(b)]=b
np.savez_compressed(p/'multinorm9.npz',Q=Q,B=B,norm=norm,X=X.astype(np.int8));print('PASS norm3',dict(h3),'norm18 low',[(k,h18[k]) for k in sorted(h18) if k<=9],'pool',len(Q),'seconds',round(time.time()-start,1))
