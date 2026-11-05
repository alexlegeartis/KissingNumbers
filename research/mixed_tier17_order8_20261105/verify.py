"""Independently expand and check the full 17D construction."""
import sys,numpy as np
p='research/mixed_tier17_eval_subgroup_20261105/cache/';r='research/mixed_tier17_eval_subgroup_20261105/seed/';own='research/mixed_tier17_order8_20261105/cache/';N=7936
old=np.fromfile(p+'mixed17_eval4_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);H=np.fromfile(p+'mixed17_eval4_group.bin',np.uint16);F=np.fromfile(p+'mixed17_eval4_flat_reps.bin',np.uint16);K=[int(h) for h in H if not(int(h)&2)];h0=next(int(h) for h in H if int(h)&2);assert len(K)==8
H0=np.fromfile(r+'mixed17_rank4_group.bin',np.uint16);S0=np.array([[1-2*(int(h)>>i&1) for i in range(16)] for h in H0],np.int16);E=np.fromfile(r+'mixed17_eq_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);fixed=[]
for e in E:
 for v in {tuple(e*s) for s in S0}:fixed.append((*[3*x for x in v],0))
assert len(fixed)==4320;fixed.extend([tuple([0]*16+[12]),tuple([0]*16+[-12])])
if len(sys.argv)>1:
 z=np.load(sys.argv[1]);ids=np.r_[z['chosen'],z['kept']].astype(int)
else:
 pos=np.fromfile(p+'baseline_seed.bin',np.int32);neg=np.loadtxt(own+'negative_seed.txt',dtype=np.int32);ids=np.array([2*int(v)+j for v in pos for j in (0,1)]+[N+2*int(v)+j for v in neg for j in (0,1)])
assert len(ids)==len(set(map(int,ids))) and np.all((0<=ids)&(ids<2*N));out=[]
for idx in ids:
 sign=1 if idx<N else -1;i=int(idx%N)
 for h in K:
  if i<256:
   word=int(F[i//2])^(h0 if i%2 else 0)^h;v=[-4 if word>>j&1 else 4 for j in range(16)]+[4*sign]
  else:
   q=i-256;row=old[q//2];word=(h0 if q%2 else 0)^h;v=[3*int(row[j])*(-1 if word>>j&1 else 1) for j in range(16)]+[6*sign]
  out.append(v)
Y=np.r_[fixed,np.asarray(out,np.int16)].astype(np.int64);g=np.r_[np.ones(16,dtype=np.int64),2];assert np.all(((Y*Y)*g).sum(1)==288) and len(np.unique(Y,axis=0))==len(Y)
for start in range(0,len(Y),96):
 D=(Y[start:start+96]*g)@Y.T;ix=np.arange(start,min(start+96,len(Y)));D[np.arange(len(ix)),ix]=0;assert D.max()<=144
print('EXACT VALID 17D',len(Y),'points; no novelty claim')
