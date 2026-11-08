"""Exact population certificates for two elementary Rankin fiber ceilings."""
from pathlib import Path
from itertools import product,combinations
from collections import Counter
import numpy as np,json,hashlib
ROOT=Path(__file__).resolve().parents[2];out=[]
for d in (27,28):
 p=next((ROOT/('verifications/improved/dim27-31-'+('algebraic-completion' if d==27 else 'joint-refills')+'/data')).glob(f'coordinates{d}_*.npz'))
 X=np.load(p)['points'];norm=np.sum(X[:,24:].astype(np.int64)**2,axis=1);heads=X[norm==64];tail=Counter(map(tuple,heads[:,24:]));alphabet=[]
 for j in range(d-24):
  for s in (-8,8):
   v=[0]*(d-24);v[j]=s;alphabet.append(tuple(v))
 for support in combinations(range(d-24),4):
  for signs in product((-4,4),repeat=4):
   v=[0]*(d-24)
   for j,s in zip(support,signs):v[j]=s
   alphabet.append(tuple(v))
 assert set(tail)==set(alphabet) and set(tail.values())=={48}
 assert np.all(np.sum(heads[:,:24].astype(np.int64)**2,axis=1)==64)
 for t in tail:
  H=heads[np.all(heads[:,24:]==t,axis=1),:24].astype(np.int64);P=H@H.T;np.fill_diagonal(P,0);assert np.max(P)==0
 out.append({'dimension':d,'source':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'baseline':len(X),'tail_count':len(tail),'points_per_tail':48,'head_points':len(heads),'family_upper_bound':len(X)-len(heads)+48*len(tail)})
report={'results':out,'proof':'For a fixed norm64 tail, norm128 points have core squared norm64 and distinct core products at most0. A unit Gram matrix G of m such core points has -1<=Gij<=0 off diagonal, so tr(G^2)<=2m using sum(Gij)>=0. Cauchy-Schwarz gives rank(G)>=m/2, while rank(G)<=24. Hence m<=48. The displayed alphabets are saturated. Deleting any outside points cannot raise this bound.','scope':'Allows arbitrary real core bodies on each of the stated fixed tails and deletion of ANY incumbent points, but no new tail directions or other radii.'}
Path(__file__).with_name('fiber_capacity.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
