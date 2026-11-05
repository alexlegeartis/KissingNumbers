import numpy as np,itertools,collections,time
from pathlib import Path
z=np.load(Path(__file__).with_name('source.npz'));X=z['X'].astype(np.int16);g=z['g'].astype(np.int16);Q=[];B=[];hist=collections.Counter()
for i in range(14):
 for j,k in itertools.combinations([a for a in range(14) if a!=i],2):
  V=np.zeros((16,15),np.int16)
  for s in range(16):V[s,[i,j,k,14]]=[8*(1-2*(s&1)),2*(1-2*((s>>1)&1)),2*(1-2*((s>>2)&1)),6*(1-2*((s>>3)&1))]
  D=(V*g)@X.T
  for v,d in zip(V,D):
   ix=np.flatnonzero(d>24);hist[len(ix)]+=1
   if len(ix)<=10:Q.append(v*4);B.append(np.r_[ix,np.full(10-len(ix),65535)].astype(np.uint16))
print('hist',sorted(hist.items()),'pool',len(Q),flush=True);np.savez_compressed('generated_den3_pool.npz',Q=np.array(Q,np.int16),B=np.array(B,np.uint16))
