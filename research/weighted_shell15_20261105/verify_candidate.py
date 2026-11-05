import numpy as np,sys
z=np.load(sys.argv[1]);X=z['X'].astype(np.int64);g=z['g'].astype(np.int64);d=int(z['den']);assert len(X)>2564 and np.all((X*X*g).sum(1)==16*d*d)
for k in range(0,len(X),128):
 G=(X[k:k+128]*g)@X.T;ii=np.arange(k,k+len(G))[:,None];jj=np.arange(len(X))[None,:];assert not np.any((G>8*d*d)&(ii<jj))
print('CERTIFIED 15D >=',len(X),flush=True)
