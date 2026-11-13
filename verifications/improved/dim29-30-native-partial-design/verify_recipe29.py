import sys
if sys.flags.optimize: raise RuntimeError('Run verification without -O')
from pathlib import Path
from itertools import combinations,product
import json,numpy as np
from parent_recipe import reconstruct
O=Path(__file__).parent/'data'
recipe=json.loads((O/'parent_recipe.json').read_text());parent=reconstruct(recipe)[0]
blocks=[tuple(2*g+s for g,s in zip(groups,signs)) for groups in combinations(range(5),4) for signs in product([0,1],repeat=4)]
un=[]
for b in blocks:
 for signs in product([-2,2],repeat=4):
  w=np.zeros(15,dtype=np.int64);w[list(b)]=signs;z=np.zeros(13,dtype=np.int64)
  for i in range(7):z[i]=w[2*i]+w[2*i+1]
  for i in range(5):z[8+i]=w[2*i]-w[2*i+1]
  un.append(z)
un=np.array(un);old=[i for i,b in enumerate(blocks) if list(b) in recipe['transversal_blocks']]
keys={tuple(r) for r in un.reshape(-1,16,13)[old].reshape(-1,13)}
section=np.delete(parent[np.all(parent[:,13:]==0,axis=1)],[13,14],axis=1);fixed=np.array([r for r in section if tuple(r) not in keys]);assert len(fixed)==568
m=np.load(O/'native29module.npz');selected=m['selected'];assert len(selected)==16 and len(set(selected))==16
triples=[t for j in selected for t in combinations(blocks[j],3)];assert len(set(triples))==len(triples)
a=np.r_[fixed,un.reshape(-1,16,13)[selected].reshape(-1,13)];rot=np.load(O/'orientation30.npz');r=np.eye(13,dtype=np.int64)*4;r[:8,:8]=rot['numerator'][:8,:8];d=4;assert np.array_equal(r@r.T,d*d*np.eye(13,dtype=np.int64));assert np.all(a@r.T%d==0);a=a@r.T//d
assert np.array_equal(a,m['points']);assert len(np.unique(a,axis=0))==824;assert np.all(np.sum(a*a,axis=1)==32)
g=a@a.T;np.fill_diagonal(g,-32);assert g.max()<=16
print(json.dumps({'pass':True,'module_points':len(a),'selected_transversal_blocks':len(selected),'norm':32,'maximum_dot':int(g.max())}))

z=np.load(O/'oriented29_heads.npz');b=a[np.sum(a[:,8:]**2,axis=1)==16];r=z['numerator'];assert np.array_equal(r.T@r,16*np.eye(13,dtype=np.int64));assert np.array_equal(b@r,4*z['points']);print('Exact 29D discovery orientation PASS')
