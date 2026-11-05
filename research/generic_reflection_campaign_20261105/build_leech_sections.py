"""Construct 22D and 23D Leech sections using the extended Golay code."""
import itertools,pathlib,numpy as np
root=pathlib.Path('research/generic_reflection_campaign_20261105/cache');root.mkdir(parents=True,exist_ok=True);gen=sum(1<<i for i in [0,1,5,6,7,9,11]);basis=[]
for i in range(12):
 v=gen<<i;v|=(v.bit_count()%2)<<23;basis.append(v)
C=[0]
for v in basis:C += [x^v for x in C]
assert len(C)==4096 and len(set(C))==4096;octads=[x for x in C if x.bit_count()==8];assert len(octads)==759;rows=[]
for i,j in itertools.combinations(range(24),2):
 for a,b in itertools.product((-4,4),repeat=2):
  v=np.zeros(24,np.int8);v[i]=a;v[j]=b;rows.append(v)
for o in octads:
 idx=[i for i in range(24) if o>>i&1]
 for mask in range(256):
  if mask.bit_count()%2:continue
  v=np.zeros(24,np.int8);v[idx]=[2*(1-2*(mask>>j&1)) for j in range(8)];rows.append(v)
for word in C:
 v=np.array([1-2*(word>>j&1) for j in range(24)],np.int8)
 for i in range(24):
  w=v.copy();w[i]*=-3;rows.append(w)
X=np.array(rows,np.int8);assert X.shape==(196560,24) and np.all((X.astype(np.int32)**2).sum(1)==32)
for d in (22,23):
 take=(X[:,0]+X[:,1]==0)
 if d==22:take &= X[:,0]+X[:,2]==0
 V=X[take];Y=np.c_[V[:,0],V[:,3:] if d==22 else V[:,2:]];g=np.r_[3 if d==22 else 2,np.ones(d-1,dtype=np.int64)];assert len(Y)==({22:49896,23:93150}[d]) and np.all(((Y.astype(np.int64)**2)*g).sum(1)==32)
 np.savez_compressed(root/f'seed{d}.npz',X=Y,g=g);print('materialized',d,'D',len(Y),'points; generic search independently verifies')
