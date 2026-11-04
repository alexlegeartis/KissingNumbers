import numpy as np
R='research/mixed_tier17_20261104/';H=np.fromfile(R+'mixed17_rank4_group.bin',np.uint16);S=np.array([[(-1 if int(h)>>i&1 else 1) for i in range(16)] for h in H],np.int16);E=np.fromfile(R+'mixed17_eq_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);T=np.fromfile(R+'mixed17_seed_tier.bin',np.int8).reshape(-1,16).astype(np.int16);V=[]
for e in E:
 for y in {tuple(e*s) for s in S}:V.append((*[3*a for a in y],0))
assert len(V)==4320
for h in (-6,6):
 for t in T:
  for y in {tuple(t*s) for s in S}:V.append((*[3*a for a in y],h))
for side in ('plus','minus'):
 f=np.fromfile(R+f'mixed17_seed_flat{side}.bin',np.uint16);height=4 if side=='plus' else -4
 for w in f:
  for z in {int(w)^int(h) for h in H}:V.append((*[(-4 if z>>i&1 else 4) for i in range(16)],height))
V.extend([tuple([0]*16+[12]),tuple([0]*16+[-12])]);X=np.asarray(V,np.int16);assert X.shape==(5730,17) and len(set(map(tuple,X.tolist())))==5730;assert np.all((X[:,:16].astype(np.int32)**2).sum(1)+2*X[:,16].astype(np.int32)**2==288)
maxdot=-999
for i in range(0,len(X),256):
 G=X[i:i+256,:16]@X[:,:16].T+2*np.outer(X[i:i+256,16],X[:,16]);
 for j in range(len(G)):G[j,i+j]=-32768
 maxdot=max(maxdot,int(G.max()))
assert maxdot<=144;print('PASS 5730 distinct points, exact scaled norm 288, maximum off-diagonal scaled inner product',maxdot)
