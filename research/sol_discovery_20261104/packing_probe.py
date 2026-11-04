import itertools, math, time, json
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix
for n in [7,8,9,10,11,12,13,14,15,16]:
 blocks=list(itertools.combinations(range(n),4));triples={a:i for i,a in enumerate(itertools.combinations(range(n),3))}
 rr=[];cc=[]
 for j,b in enumerate(blocks):
  for a in itertools.combinations(b,3):rr.append(triples[a]);cc.append(j)
 A=coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(triples),len(blocks))).tocsc()
 t=time.monotonic()
 r=milp(-np.ones(len(blocks)),integrality=np.ones(len(blocks)),bounds=Bounds(0,1),constraints=LinearConstraint(A,0,1),options={'time_limit':5.0,'mip_rel_gap':0.0})
 sel=np.flatnonzero(r.x>0.5).tolist() if r.x is not None else []
 print(json.dumps({'n':n,'blocks':len(sel),'points':16*len(sel)+2*n,'status':int(r.status),'seconds':time.monotonic()-t,'bound':getattr(r,'mip_dual_bound',None)}),flush=True)
