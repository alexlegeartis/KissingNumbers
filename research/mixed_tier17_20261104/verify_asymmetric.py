"""Independently expand an asymmetric 17D exchange and check all pairs."""
import sys,numpy as np
R='research/mixed_tier17_20261104/';d=np.load(sys.argv[1]);chosen=set(map(int,d['new_indices']));removed=set(map(int,d['removed_indices']));gain=int(round(float(d['gain'])));H=np.fromfile(R+'mixed17_rank4_group.bin',np.uint16);S=np.array([[(-1 if int(h)>>i&1 else 1) for i in range(16)] for h in H],np.int16);E=np.fromfile(R+'mixed17_eq_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);T=np.fromfile(R+'mixed17_seed_tier.bin',np.int8).reshape(-1,16).astype(np.int16);C=np.fromfile(R+'mixed17_six_all.bin',np.int8).reshape(-1,16).astype(np.int16);F=np.fromfile(R+'mixed17_equator_frontier.bin',np.int32).reshape(-1,6);P=np.fromfile(R+'mixed17_seed_flatplus.bin',np.uint16);M=np.fromfile(R+'mixed17_seed_flatminus.bin',np.uint16);n=len(F);V=[]
for i,e in enumerate(E):
 if i not in removed:
  for y in {tuple(e*s) for s in S}:V.append((*[3*a for a in y],0))
for side,height in enumerate((6,-6)):
 for i,t in enumerate(T):
  if 480+32*side+i not in removed:
   for y in {tuple(t*s) for s in S}:V.append((*[3*a for a in y],height))
for k,flat in enumerate((P,M)):
 for i,w in enumerate(flat):
  if 544+12*k+i not in removed:
   for z in {int(w)^int(h) for h in H}:V.append((*[(-4 if z>>j&1 else 4) for j in range(16)],4 if k==0 else -4))
for i in chosen:
 j=i%n;c=C[int(F[j,0])];height=6 if i<n else -6
 for y in {tuple(c*s) for s in S}:V.append((*[3*a for a in y],height))
V.extend([tuple([0]*16+[12]),tuple([0]*16+[-12])]);X=np.asarray(V,np.int16);assert len(X)==5730+gain and len(set(map(tuple,X.tolist())))==len(X);assert np.all((X[:,:16].astype(np.int32)**2).sum(1)+2*X[:,16].astype(np.int32)**2==288)
mx=-999
for i in range(0,len(X),256):
 G=X[i:i+256,:16]@X[:,:16].T+2*np.outer(X[i:i+256,16],X[:,16]);
 for j in range(len(G)):G[j,i+j]=-32768
 mx=max(mx,int(G.max()))
assert mx<=144;np.savez_compressed(R+'mixed17_verified_asymmetric.npz',coordinates_scaled6=X,metric_last=2,norm_squared=288);print('PASS:',len(X),'points; exact maximum scaled inner product',mx)
