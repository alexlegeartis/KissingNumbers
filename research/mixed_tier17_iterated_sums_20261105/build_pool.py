from pathlib import Path
Path('research/mixed_tier17_iterated_sums_20261105/cache/').mkdir(parents=True,exist_ok=True)
"""Generate new 17D equator orbits by exact 120-degree vector sums."""
import numpy as np,collections
z=np.load('research/mixed_tier17_equator_exchange_20261105/cache/pool.npz');X=z['X'].astype(np.int16);E=X[X[:,16]==0,:16].astype(np.int16);R=np.load('research/mixed_tier17_root_sums_20261105/cache/root_pool.npz');offs=np.r_[0,np.cumsum(R['sizes'])];eligible=np.flatnonzero(R['cost']<=56);Q=R['points'][offs[eligible],:16].astype(np.int16);weights=z['weights'];point_group=z['point_group'];H=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_group.bin',np.uint16);K=[int(h) for h in H if not int(h)&6];rng=np.random.default_rng(171121);subset=rng.choice(len(Q),min(3000,len(Q)),replace=False);W=[]
for st in range(0,len(subset),64):
 U=Q[subset[st:st+64]].astype(np.int32);D=U@E.astype(np.int32).T
 for j in range(len(U)):
  ids=np.flatnonzero(D[j]==-144)
  if len(ids)>32:ids=rng.choice(ids,32,replace=False)
  if len(ids):W.append(U[j]+E[ids])
 if st%512==0:print('pairs processed',st,'candidate batches',len(W),flush=True)
W=np.unique(np.concatenate(W).astype(np.int16),axis=0);assert np.all((W.astype(np.int32)**2).sum(1)==288);old=set(map(tuple,X[:,:16].tolist()));seen=set();orbits=[]
for w in W:
 t=tuple(map(int,w));orb={tuple((-1 if h>>i&1 else 1)*t[i] for i in range(16)) for h in K};key=min(orb)
 if key in seen or any(v in old for v in orb):continue
 seen.add(key);orbits.append(np.array([list(x)+[0] for x in sorted(orb)],np.int16))
print('new sum orbits',len(orbits),'sizes',collections.Counter(map(len,orbits)),flush=True);hist=collections.Counter();cost=[];lists=[]
for i,orb in enumerate(orbits):
 d=orb[:,:16].astype(np.int32)@X[:,:16].astype(np.int32).T;bad=np.flatnonzero(np.any(d>144,axis=0));groups=np.unique(point_group[bad]);assert int(weights[groups].sum())==len(bad);lists.append(groups);cost.append(len(bad));hist[len(bad)]+=1
 if i%2000==0:print('checked',i,'best blocker points',min(cost),flush=True)
indptr=np.r_[0,np.cumsum([len(a) for a in lists])];np.savez_compressed('research/mixed_tier17_iterated_sums_20261105/cache/root_pool2.npz',X=X,points=np.concatenate(orbits),sizes=np.array(list(map(len,orbits))),weights=weights,point_group=point_group,indptr=indptr,indices=np.concatenate(lists),cost=np.array(cost));print('hist',sorted(hist.items())[:25],flush=True)
