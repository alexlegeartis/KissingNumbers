"""Expand proposed positive-cap orbits and check the complete 17D code exactly."""
import argparse
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--selection',help='Text file of orbit IDs. Default: baseline 44 IDs.')
a=p.parse_args()
r=Path('research/mixed_tier17_eval_subgroup_20261105/seed')
c=Path('research/mixed_tier17_eval_subgroup_20261105/cache')
ids=np.loadtxt(a.selection,dtype=int,ndmin=1) if a.selection else np.fromfile(c/'baseline_seed.bin',np.int32)
assert len(set(ids.tolist()))==len(ids) and np.all((ids>=0)&(ids<3968))
H=np.fromfile(c/'mixed17_eval4_group.bin',np.uint16)
V=np.fromfile(c/'mixed17_eval4_orbits.bin',np.int8).reshape(-1,16).astype(np.int16)
F=np.fromfile(c/'mixed17_eval4_flat_reps.bin',np.uint16)
H0=np.fromfile(r/'mixed17_rank4_group.bin',np.uint16)
S0=np.array([[(-1 if int(h)>>i&1 else 1) for i in range(16)] for h in H0],np.int16)
S=np.array([[(-1 if int(h)>>i&1 else 1) for i in range(16)] for h in H],np.int16)
E=np.fromfile(r/'mixed17_eq_orbits.bin',np.int8).reshape(-1,16).astype(np.int16)
T=np.fromfile(r/'mixed17_seed_tier.bin',np.int8).reshape(-1,16).astype(np.int16)
rows=[]
for e in E:
 for y in {tuple(e*s) for s in S0}: rows.append((*[3*x for x in y],0))
assert len(rows)==4320
for t in T:
 for y in {tuple(t*s) for s in S0}: rows.append((*[3*x for x in y],-6))
for w in np.fromfile(r/'mixed17_seed_flatminus.bin',np.uint16):
 for z in {int(w)^int(h) for h in H0}: rows.append((*[(-4 if z>>i&1 else 4) for i in range(16)],-4))
rows.extend([tuple([0]*16+[12]),tuple([0]*16+[-12])])
fixed=len(rows)
for idx in ids:
 if idx<128:
  w=int(F[idx]); orb={w^int(h) for h in H}; assert len(orb)==16
  for z in orb: rows.append((*[(-4 if z>>i&1 else 4) for i in range(16)],4))
 else:
  orb={tuple(V[idx-128]*s) for s in S}; assert len(orb)==16
  for y in orb: rows.append((*[3*x for x in y],6))
X=np.asarray(rows,np.int32)
assert len(X)==fixed+16*len(ids)
assert len(set(map(tuple,X.tolist())))==len(X), 'duplicate vectors'
assert np.all((X[:,:16]**2).sum(1)+2*X[:,16]**2==288), 'incorrect norm'
maxdot=-10**9
for i in range(0,len(X),128):
 G=X[i:i+128,:16]@X[:,:16].T+2*np.outer(X[i:i+128,16],X[:,16])
 for j in range(len(G)): G[j,i+j]=-10**9
 maxdot=max(maxdot,int(G.max()))
assert maxdot<=144, f'INVALID code: maximum scaled dot product {maxdot} > 144'
print(f'PASS {len(X)} distinct points, norm² 288, maximum off-diagonal dot {maxdot}')
