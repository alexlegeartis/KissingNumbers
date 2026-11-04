import itertools,time,json
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
for n in [9,11,12,13,14,15,16,17,20,24]:
 for period in [n,n-1,n-2]:
  blocks=list(itertools.combinations(range(n),4));remaining=set(blocks);orbits=[];invoices=[]
  while remaining:
   b=min(remaining);orb=set();cur=b
   while cur not in orb:
    orb.add(cur);cur=tuple(sorted((i+1)%period if i<period else i for i in cur))
   remaining.difference_update(orb)
   ts=[t for s in orb for t in itertools.combinations(s,3)]
   if len(ts)==len(set(ts)):orbits.append(sorted(orb));invoices.append(ts)
  ti={t:i for i,t in enumerate(itertools.combinations(range(n),3))};rr=[];cc=[]
  for j,ts in enumerate(invoices):
   for t in ts:rr.append(ti[t]);cc.append(j)
  M=coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(ti),len(orbits))).tocsc()
  st=time.monotonic();r=milp(-np.array(list(map(len,orbits))),integrality=np.ones(len(orbits)),bounds=Bounds(0,1),constraints=LinearConstraint(M,0,1),options={'time_limit':1.0})
  selected=[] if r.x is None else np.flatnonzero(r.x>.5)
  print(json.dumps({'n':n,'period':period,'supports':sum(len(orbits[i]) for i in selected),'status':int(r.status),'seconds':time.monotonic()-st}),flush=True)
