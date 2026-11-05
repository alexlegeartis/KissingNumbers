"""Exact finite conflict graph for low-blocker equator orbits and source groups."""
import numpy as np
z=np.load('research/mixed_tier17_equator_exchange_20261105/cache/pool.npz');Q=z['points'].reshape(-1,4,17).astype(np.int16);cost=z['cost'];cand=np.flatnonzero(cost<=56);Q=Q[cand];V=Q[:,0,:16].astype(np.int32);weights=z['weights'].astype(np.int8);ng=len(weights);indptr=z['indptr'];indices=z['indices'];H=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_group.bin',np.uint16);K=[int(h) for h in H if not int(h)&6];sign=np.array([[1-2*(h>>j&1) for j in range(16)] for h in K],np.int32);assert len(K)==4
f=open('research/mixed_tier17_equator_exchange_20261105/cache/full_edges.bin','wb');count=0
for i,c in enumerate(cand):
 old=indices[indptr[c]:indptr[c+1]];E=np.c_[old,np.full(len(old),ng+i,np.int32)].astype(np.int32);E.tofile(f);count+=len(E)
for st in range(0,len(V),128):
 U=V[st:st+128];D=np.zeros((len(U),len(V)),bool)
 for sg in sign:D|=(U*sg)@V.T>144
 i,j=np.where(D);i+=st;sel=i<j;E=np.c_[i[sel]+ng,j[sel]+ng].astype(np.int32);E.tofile(f);count+=len(E)
f.close();W=np.r_[weights,4*np.ones(len(cand),np.int8)].astype(np.int8);W.tofile('research/mixed_tier17_equator_exchange_20261105/cache/full_weights.bin');origin=np.r_[np.ones(ng,np.uint8),np.zeros(len(cand),np.uint8)];origin.tofile('research/mixed_tier17_equator_exchange_20261105/cache/full_origin.bin');np.save('research/mixed_tier17_equator_exchange_20261105/cache/full_candidate_ids.npy',cand);print('vertices',len(W),'edges',count,'source groups',ng,'source weight',int(weights.sum()),flush=True)
