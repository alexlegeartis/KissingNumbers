"""Build shared-subgroup equator candidates from blockwise coordinate maps."""
import numpy as np,collections,runpy,sys
from pathlib import Path
Path('research/mixed_tier17_equator_exchange_20261105/cache').mkdir(parents=True,exist_ok=True)
saved=sys.argv;sys.argv=['verify.py'];seed=runpy.run_path('research/mixed_tier17_order4_20261105/verify.py');sys.argv=saved;X=seed['Y'].astype(np.int16);E=X[X[:,16]==0];H=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_group.bin',np.uint16);K=[int(h) for h in H if not int(h)&6];assert len(K)==4
classes={}
for j in range(16):classes.setdefault(tuple((h>>j)&1 for h in K),[]).append(j)
rng=np.random.default_rng(170425);images=[]
for t in range(24):
 perm=np.arange(17)
 for block in classes.values():perm[block]=rng.permutation(block)
 Y=E[:,perm];images.append(Y)
Q=np.unique(np.concatenate(images),axis=0);source=set(map(tuple,X.tolist()));Q=np.array([x for x in Q if tuple(x) not in source],np.int16);print('distinct new equator points',len(Q),flush=True)
remaining=set(map(tuple,Q.tolist()));orbits=[]
while remaining:
 x=next(iter(remaining));orb={tuple((-1 if h>>i&1 else 1)*x[i] for i in range(16))+(0,) for h in K};assert orb<=remaining;remaining-=orb;orbits.append(np.array(sorted(orb),np.int16))
print('candidate orbits',len(orbits),'size histogram',collections.Counter(map(len,orbits)),flush=True)
# Decompose the complete 17D source, including the two caps, into K-orbits.
lookup={tuple(x):i for i,x in enumerate(X.tolist())};unseen=set(lookup);groups=[];point_group=np.empty(len(X),np.int32)
while unseen:
 x=next(iter(unseen));orb={tuple((-1 if h>>i&1 else 1)*x[i] for i in range(16))+(x[16],) for h in K};assert orb<=unseen;unseen-=orb;ids=[lookup[y] for y in orb];point_group[ids]=len(groups);groups.append(ids)
weights=np.array(list(map(len,groups)),np.int32);print('source groups',len(groups),'hist',collections.Counter(weights),flush=True)
# An orbit's blockers are the union over all its points; old points are checked exactly.
hist=collections.Counter();cost=[];lists=[]
for k,orb in enumerate(orbits):
 d=(orb[:,:16].astype(np.int32)@X[:,:16].astype(np.int32).T);bad=np.flatnonzero(np.any(d>144,axis=0));blockers=np.unique(point_group[bad]);assert weights[blockers].sum()==len(bad);cost.append(len(bad));lists.append(blockers);hist[len(bad)]+=1
 if k%2000==0:print('checked',k,'min',min(cost),flush=True)
print('blocker histogram',sorted(hist.items())[:25],flush=True)
indptr=np.r_[0,np.cumsum([len(a) for a in lists])];indices=np.concatenate(lists)
np.savez_compressed('research/mixed_tier17_equator_exchange_20261105/cache/pool.npz',X=X,points=np.concatenate(orbits),sizes=np.array(list(map(len,orbits))),cost=np.array(cost),weights=weights,point_group=point_group,indptr=indptr,indices=indices,permutations=len(images))
