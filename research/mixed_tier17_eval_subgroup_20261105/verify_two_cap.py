"""Independently expand and check a two-cap search output."""
import numpy as np,sys
r='research/mixed_tier17_eval_subgroup_20261105/';p=r+'cache/';H=np.fromfile(p+'mixed17_eval4_group.bin',np.uint16);S=np.array([[1-2*(int(h)>>i&1) for i in range(16)] for h in H],np.int16);V=np.fromfile(p+'mixed17_eval4_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);F=np.fromfile(p+'mixed17_eval4_flat_reps.bin',np.uint16);H0=np.fromfile(r+'seed/mixed17_rank4_group.bin',np.uint16);S0=np.array([[1-2*(int(h)>>i&1) for i in range(16)] for h in H0],np.int16);E=np.fromfile(r+'seed/mixed17_eq_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);rows=[]
for e in E:
 for y in {tuple(e*s) for s in S0}:rows.append((*[3*x for x in y],0))
assert len(rows)==4320
rows.extend([tuple([0]*16+[12]),tuple([0]*16+[-12])])
ids=np.load(sys.argv[1]);assert len(np.unique(ids))==len(ids)
for i in ids:
 side=1 if i<3968 else -1;idx=int(i)%3968
 if idx<128:
  for z in {int(F[idx])^int(h) for h in H}:rows.append((*[(-4 if z>>k&1 else 4) for k in range(16)],4*side))
 else:
  for y in {tuple(V[idx-128]*s) for s in S}:rows.append((*[3*x for x in y],6*side))
X=np.array(rows,np.int32);assert len(X)==4322+16*len(ids);assert len(set(map(tuple,X)))==len(X);assert np.all((X[:,:16]**2).sum(1)+2*X[:,16]**2==288)
maxdot=0
for i in range(0,len(X),128):
 G=X[i:i+128,:16]@X[:,:16].T+2*np.outer(X[i:i+128,16],X[:,16]);
 for j in range(len(G)):G[j,i+j]=-1
 maxdot=max(maxdot,int(G.max()))
assert maxdot<=144,f'INVALID max dot {maxdot}'
print('PASS',len(X),'points max scaled dot',maxdot)
