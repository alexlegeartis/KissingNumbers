"""Finite exact-graph feasibility search; MILP statuses are not proof certificates."""
from common import *
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack
import warnings
import time
import argparse

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=300)
    p.add_argument('--target',type=int,default=842);p.add_argument('--threads',type=int,default=8)
    p.add_argument('--pool',type=Path,default=ROOT/'runs'/'trade_pool')
    a=p.parse_args();reserve();folder=a.pool;h=np.fromfile(folder/'header.bin',dtype=np.uint32);n=int(h[0])
    bad=np.fromfile(folder/'weights.bin',dtype=np.uint8).reshape(n,n);i,j=np.nonzero(np.triu(bad,1));m=len(i)
    matrix=coo_matrix((np.ones(2*m),(np.r_[np.arange(m),np.arange(m)],np.r_[i,j])),shape=(m,n)).tocsc()
    matrix=vstack([matrix,coo_matrix(np.ones((1,n)))]).tocsc();lo=np.r_[np.full(m,-np.inf),a.target];hi=np.r_[np.ones(m),np.inf]
    start=time.monotonic()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected')
        result=milp(np.zeros(n),integrality=np.ones(n),bounds=Bounds(0,1),constraints=LinearConstraint(matrix,lo,hi),
                    options=dict(time_limit=a.seconds,threads=a.threads,presolve=True))
    row=dict(target=a.target,pool=str(folder),vertices=n,status=int(result.status),message=str(result.message),
              elapsed_seconds=time.monotonic()-start,found=False)
    if result.x is not None:
        selected=np.flatnonzero(result.x>.5);assert len(selected)>=a.target and not bad[np.ix_(selected,selected)].any()
        row.update(found=True,selected=selected.tolist());np.save(ROOT/'runs'/f'{folder.name}_milp_{a.target}.npy',selected)
    save_json(ROOT/f'{folder.name}_milp_{a.target}.json',row);print(json.dumps(row),flush=True)
