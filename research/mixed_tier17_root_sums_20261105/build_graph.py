"""Build joint orbit graph with new equator AND new cap families."""
import numpy as np
z=np.load('research/mixed_tier17_equator_exchange_20261105/cache/pool.npz');X=z['X'].astype(np.int16);ng=len(z['weights']);point_group=z['point_group'];Q=z['points'].reshape(-1,4,17).astype(np.int16);eqids=np.flatnonzero(z['cost']<=56);Q=Q[eqids];T=np.load('research/mixed_tier17_joint_layers_20261105/cache/cap_pool.npz');capids=np.flatnonzero(T['cost']<=8);N=15872;V=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);H=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_group.bin',np.uint16);F=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_flat_reps.bin',np.uint16);K=[int(h) for h in H if not int(h)&6];h0=next(int(h) for h in H if (int(h)&6)==2);h1=next(int(h) for h in H if (int(h)&6)==4);reps=[0,h0,h1,h0^h1];sg=np.array([[1-2*(h>>j&1) for j in range(16)]+[1] for h in K],np.int32);g=np.r_[np.ones(16,dtype=np.int32),2];cap=[]
for idx in capids:
 side=1 if idx<N else -1;i=int(idx%N)
 if i<512:
  word=int(F[i//4])^reps[i%4];v=[-4 if word>>j&1 else 4 for j in range(16)]+[4*side]
 else:
  q=i-512;row=V[q//4];word=reps[q%4];v=[3*int(row[j])*(-1 if word>>j&1 else 1) for j in range(16)]+[6*side]
 cap.append(v)
cap=np.asarray(cap,np.int16);Z=np.load('research/mixed_tier17_root_sums_20261105/cache/root_pool.npz');rootids=np.flatnonzero(Z['cost']<=36);offs=np.r_[0,np.cumsum(Z['sizes'])];rootrep=Z['points'][offs[rootids]];R=np.r_[Q[:,0,:],cap,rootrep].astype(np.int32);M=len(R);assert len(eqids)==2560 and len(capids)==4032 and np.array_equal(Z['X'],z['X']) and np.array_equal(Z['point_group'],point_group) and np.all(((R*R)*g).sum(1)==288);print('source groups',ng,'eq candidates',len(Q),'cap candidates',len(cap),'root-sum candidates',len(rootids),'total',ng+M,flush=True)
f=open('research/mixed_tier17_root_sums_20261105/cache/edges.bin','wb');count=0
# Complete candidate-to-source blockers for each full candidate orbit. Since
# source is invariant under K, a fixed representative against all source
# points covers every source/candidate orbit conflict.
for st in range(0,M,64):
 D=(R[st:st+64]*g)@X.astype(np.int32).T
 for q in range(len(D)):
  bad=np.flatnonzero(D[q]>144);groups=np.unique(point_group[bad]);E=np.c_[groups,np.full(len(groups),ng+st+q,np.int32)].astype(np.int32);E.tofile(f);count+=len(E)
 if st%1024==0:print('source pairs',st,'edges',count,flush=True)
# Candidate-candidate compatibility by exact relative subgroup products.
for st in range(0,M,64):
 U=R[st:st+64];D=np.zeros((len(U),M),bool)
 for s in sg:D|=((U*s)*g)@R.T>144
 i,j=np.where(D);i+=st;take=i<j;E=np.c_[i[take]+ng,j[take]+ng].astype(np.int32);E.tofile(f);count+=len(E)
 if st%1024==0:print('candidate pairs',st,'edges',count,flush=True)
f.close();weights=np.r_[z['weights'].astype(np.uint8),4*np.ones(len(eqids)+len(capids),np.uint8),Z['sizes'][rootids].astype(np.uint8)];weights.tofile('research/mixed_tier17_root_sums_20261105/cache/weights.bin');origin=np.r_[np.ones(ng,np.uint8),np.zeros(M,np.uint8)];origin.tofile('research/mixed_tier17_root_sums_20261105/cache/origin.bin');np.savez_compressed('research/mixed_tier17_root_sums_20261105/cache/nodes.npz',eqids=eqids,capids=capids,caprep=cap,rootids=rootids,ng=ng);print('COMPLETE edges',count,'source weight',int(weights[:ng].sum()),flush=True)
