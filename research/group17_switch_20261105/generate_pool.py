import numpy as np
X=np.load('coordinates.npz')['source'].astype(np.int32);g=np.array([1]*16+[2],np.int32);u=np.zeros(17,np.int32);v=u.copy();w=u.copy();u[[0,16]]=1;v[[1,16]]=1;w[[2,3,4,5]]=1;gen=[u,v,w]
def refl(Q,z):
 n=int((z*g)@z);P=2*np.outer(Q@(z*g),z);assert not np.any(P%n);return Q-P//n
A=X;imgs=[]
for word in [[],[0],[1],[2],[0,1],[0,2],[1,2],[0,1,2]]:
 Z=A.copy()
 for k in word:Z=refl(Z,gen[k])
 imgs.append(Z)
C=np.unique(np.concatenate(imgs),axis=0);norm=93312;assert np.all(((C.astype(np.int64)**2)*g).sum(1)==norm);S=set(map(tuple,A));o=np.array([tuple(x) in S for x in C],np.uint8);assert o.sum()==5730;E=[]
for st in range(0,len(C),96):
 dot=(C[st:st+96]*g)@C.T;i,j=np.where(dot>norm//2);i+=st;t=i<j;E.append(np.c_[i[t],j[t]].astype(np.int32))
E=np.concatenate(E);np.savez_compressed('search17w_multiimage_graph.npz',C=C,metric=g,norm=norm,origin=o,edges=E);o.tofile('search17w_multiimage_origin.bin');E.tofile('search17w_multiimage_edges.bin');print('17 weighted mixed-coordinate images',len(C),'edges',len(E),'source',o.sum(),flush=True)
