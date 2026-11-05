import sys,numpy as np
z=np.load(sys.argv[1]);X=z['X'].astype(np.int64);g=z['g'].astype(np.int64);norm=int(z['norm']);assert X.shape[1]==17 and tuple(g)==(1,)*16+(2,) and norm==288 and len(X)>5730;assert np.all(((X*X)*g).sum(1)==norm) and len(np.unique(X,axis=0))==len(X)
for st in range(0,len(X),96):
 D=(X[st:st+96]*g)@X.T;ix=np.arange(st,min(st+96,len(X)));D[np.arange(len(ix)),ix]=0;assert D.max()<=norm//2
print('EXACT 17D',len(X),'points; literature comparison required')
