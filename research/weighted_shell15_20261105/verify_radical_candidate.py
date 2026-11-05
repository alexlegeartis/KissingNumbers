"""Independently verify a mixed-norm spherical code in metric diag(1^14,2)."""
import numpy as np,sys
z=np.load(sys.argv[1]); A=z['normal'].astype(np.int64);B=z['radical'].astype(np.int64);g=z['g'].astype(np.int64);assert len(g)==15 and tuple(g)==(1,)*14+(2,)
assert np.all(((A*A)*g).sum(1)==2304) and np.all(((B*B)*g).sum(1)==50)
V=np.r_[A,B];N=np.r_[np.full(len(A),2304),np.full(len(B),50)];assert len(V)>2564
for s in range(0,len(V),128):
 D=(V[s:s+128]*g)@V.T;mask=(D>0)&(4*D*D>N[s:s+128,None]*N[None,:]);ix=np.arange(s,min(s+128,len(V)));mask[np.arange(len(ix)),ix]=0
 if np.any(mask):raise ValueError(('incompatible',s,np.argwhere(mask)[0]))
print('PASS',len(V),'quadratic-field vectors',len(B),'all exact pair comparisons')
