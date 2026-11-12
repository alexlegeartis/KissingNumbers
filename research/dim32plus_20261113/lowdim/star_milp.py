"""Numerical optimization on an exact fully enumerated coordinate-star pool."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict
import argparse,json,time
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack
ap=argparse.ArgumentParser();ap.add_argument('coordinate',type=int);ap.add_argument('--seconds',type=float,default=45);args=ap.parse_args()
P=Path(__file__).parent;p=args.coordinate
source=[int(x) for x in (P/'seed35.txt').read_text().split()]
old=[x for x in source if x>>p&1];fixed=[x for x in source if not x>>p&1]
V=[int(x) for x in (P/f'star35_{p}.txt').read_text().split()]
by=defaultdict(list)
for i,x in enumerate(V):
 for sub in combinations([j for j in range(35) if x>>j&1],5):by[sub].append(i)
cl=[v for v in by.values() if len(v)>1]
rr=np.repeat(np.arange(len(cl)),[len(x) for x in cl]);cc=np.fromiter((i for row in cl for i in row),np.int64)
A=coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(cl),len(V))).tocsr()
A=vstack([A,coo_matrix(-np.ones((1,len(V))))],format='csr');up=np.r_[np.ones(len(cl)),-float(len(old)+1)]
print('pool',p,len(V),'old',len(old),'rows',len(cl),'nnz',len(rr),flush=True)
t=time.monotonic();res=milp(-np.ones(len(V)),integrality=np.ones(len(V)),bounds=Bounds(0,1),constraints=LinearConstraint(A,-np.inf,up),options={'time_limit':args.seconds,'mip_rel_gap':0.,'threads':1})
chosen=[] if res.x is None else [V[i] for i in np.flatnonzero(res.x>.5)]
report={'coordinate':p,'source_supports':len(source),'removed_star_size':len(old),'complete_pool_size':len(V),'solver_status':int(res.status),'solver_message':res.message,'chosen_supports':len(chosen),'solver_seconds':time.monotonic()-t,'scope':'Solver result on finite coordinate-star neighborhood; only an independently checked feasible witness can imply a kissing lower bound.'}
if chosen:
 assert len(chosen)==len(set(chosen)) and all(x.bit_count()==8 for x in chosen)
 assert all((x&y).bit_count()<=4 for i,x in enumerate(chosen) for y in chosen[:i])
 assert all((x&y).bit_count()<=4 for x in chosen for y in fixed)
 report['verified_support_count']=len(fixed)+len(chosen)
 (P/f'star35_{p}_candidate.txt').write_text('\n'.join(map(str,fixed+chosen))+'\n')
(P/f'star35_{p}_result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
