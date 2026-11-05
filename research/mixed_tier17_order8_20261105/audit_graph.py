"""Compare graph edges with direct 17D integer Gram checks on random pairs."""
import numpy as np
N=7936;R='research/mixed_tier17_eval_subgroup_20261105/cache/';V=np.fromfile(R+'mixed17_eval4_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);H=np.fromfile(R+'mixed17_eval4_group.bin',np.uint16);F=np.fromfile(R+'mixed17_eval4_flat_reps.bin',np.uint16);K=[int(h) for h in H if not int(h)&2];h0=next(int(h) for h in H if int(h)&2)
def expand(i,side):
 out=[]
 for h in K:
  if i<256:
   w=int(F[i//2])^(h0 if i%2 else 0)^h;out.append(([-4 if w>>j&1 else 4 for j in range(16)]+[4*side]))
  else:
   q=i-256;w=(h0 if q%2 else 0)^h;out.append(([3*int(V[q//2,j])*(-1 if w>>j&1 else 1) for j in range(16)]+[6*side]))
 return np.array(out,np.int32)
E=np.fromfile('research/mixed_tier17_order8_20261105/cache/graph_edges.bin',np.int32).reshape(-1,2).astype(np.int64);a=np.minimum(E[:,0],E[:,1]);b=np.maximum(E[:,0],E[:,1]);codes=np.sort(a*N+b);del E,a,b
X=np.fromfile('research/mixed_tier17_order8_20261105/cache/cross_edges.bin',np.int32).reshape(-1,2).astype(np.int64);xcodes=np.sort(X[:,0]*N+X[:,1]);rng=np.random.default_rng(17083);g=np.r_[np.ones(16,dtype=np.int32),2]
for trial in range(2000):
 i=int(rng.integers(N));j=int(rng.integers(N));cross=bool(trial&1)
 if not cross and i==j:continue
 u=expand(i,1);v=expand(j,-1 if cross else 1);bad=bool(np.max((u*g)@v.T)>144)
 if cross:key=i*N+j;at=np.searchsorted(xcodes,key);found=at<len(xcodes) and xcodes[at]==key
 else:key=min(i,j)*N+max(i,j);at=np.searchsorted(codes,key);found=at<len(codes) and codes[at]==key
 assert bad==found,(i,j,cross,bad,found)
print('PASS 2000 random same/opposite cap 8x8 exact Gram comparisons')
