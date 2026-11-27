from tail26 import *
X,Y=np.where(pcol==0)[0],np.where(pcol==1)[0]
OLD=P[X]@P[Y].T
assert set(np.unique(OLD[OLD>32]))=={48,64}
e=csr_matrix(OLD==64)
n,h=connected_components(bmat([[None,e],[e.T,None]],format="csr"),directed=False)
g=np.empty(len(P),np.int32);g[np.r_[X,Y]]=h
group=2*g+pcol;order=np.argsort(group,kind="stable")
grp,st=np.unique(group[order],return_index=True)
stat=np.full((2*n,len(FB)),-10000000,np.int64)
stat[grp]=np.maximum.reduceat(P[order]@FB.T,st,axis=0)
