import sys
if sys.flags.optimize: raise RuntimeError('Run verification without -O')
from pathlib import Path
from itertools import combinations,product
import json,numpy as np
from parent_recipe import reconstruct
O=Path(__file__).parent/'data'
recipe=json.loads((O/'parent_recipe.json').read_text());parent=reconstruct(recipe)[0]
blocks=[tuple(2*g+s for g,s in zip(groups,signs)) for groups in combinations(range(6),4) for signs in product([0,1],repeat=4)]
un=[]
for b in blocks:
 for signs in product([-2,2],repeat=4):
  w=np.zeros(15,dtype=np.int64);w[list(b)]=signs;z=np.zeros(14,dtype=np.int64)
  for i in range(7):z[i]=w[2*i]+w[2*i+1]
  for i in range(6):z[8+i]=w[2*i]-w[2*i+1]
  un.append(z)
un=np.array(un);old=[i for i,b in enumerate(blocks) if list(b) in recipe['transversal_blocks']]
keys={tuple(r) for r in un.reshape(-1,16,14)[old].reshape(-1,14)}
section=np.delete(parent[parent[:,14]==0],14,axis=1);fixed=np.array([r for r in section if tuple(r) not in keys]);assert len(fixed)==876
m=np.load(O/'native30module.npz');selected=m['selected'];assert len(selected)==36 and len(set(selected))==36
triples=[t for j in selected for t in combinations(blocks[j],3)];assert len(set(triples))==len(triples)
a=np.r_[fixed,un.reshape(-1,16,14)[selected].reshape(-1,14)];rot=np.load(O/'orientation30.npz');r=rot['numerator'].astype(np.int64);d=int(rot['denominator']);assert np.array_equal(r@r.T,d*d*np.eye(14,dtype=np.int64));assert np.all(a@r.T%d==0);a=a@r.T//d
assert np.array_equal(a,m['points']);assert len(np.unique(a,axis=0))==1452;assert np.all(np.sum(a*a,axis=1)==32)
g=a@a.T;np.fill_diagonal(g,-32);assert g.max()<=16
print(json.dumps({'pass':True,'module_points':len(a),'selected_transversal_blocks':len(selected),'norm':32,'maximum_dot':int(g.max())}))
