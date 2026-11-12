"""Bounded numerical MILP discovery on the exact finite latitude conflict graph."""
from pathlib import Path
import json
import warnings
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix

HERE=Path(__file__).resolve().parent
P=HERE/'lifts37'
V=np.fromfile(P/'small.bin',np.int32).reshape(-1,5)
T=np.fromfile(P/'T.bin',np.int8).reshape(253,253,3)
bad=np.zeros((len(V),len(V)),bool)
for t in range(3):
    d=(T[V[:,0,None],V[None,:,0],t]
       +T[V[:,1,None],V[None,:,1],(V[None,:,3]-V[:,3,None]+t)%3]
       +T[V[:,2,None],V[None,:,2],(V[None,:,4]-V[:,4,None]+t)%3])
    bad|=d>2
ii,jj=np.where(np.triu(bad,1));edges=len(ii)
row=np.r_[np.repeat(np.arange(edges),2),np.full(len(V),edges)]
col=np.r_[np.c_[ii,jj].ravel(),np.arange(len(V))]
A=coo_matrix((np.ones(len(row)),(row,col)),shape=(edges+1,len(V))).tocsc()
lower=np.r_[np.full(edges,-np.inf),637]
upper=np.r_[np.ones(edges),np.inf]
with warnings.catch_warnings():
    warnings.simplefilter('ignore',RuntimeWarning)
    result=milp(-np.ones(len(V)),integrality=np.ones(len(V)),bounds=Bounds(0,1),
                constraints=LinearConstraint(A,lower,upper),
                options={'time_limit':60,'threads':1,'mip_rel_gap':0})
report={'vertices':len(V),'edges':edges,'required_orbits':637,
        'status':int(result.status),'message':result.message,
        'scope':'Numerical MILP discovery on the <=4-blocker pool; not a proof of a global construction ceiling.'}
if result.x is not None:
    ids=np.flatnonzero(result.x>.5)
    assert not np.any(np.triu(bad[np.ix_(ids,ids)],1))
    assert len(ids)>=637
    V[ids].tofile(P/'milp_candidate.bin')
    report['candidate_orbits']=len(ids)
(P/'milp_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
