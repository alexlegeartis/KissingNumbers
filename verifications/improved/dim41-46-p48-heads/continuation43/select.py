"""Select a finite packing from the shipped rational proposals.

This optimizer is a discovery/reproduction aid. Verify the resulting explicit
configuration independently; neither an objective value nor solver status is
a kissing-number certificate.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.resolve()==(HERE/'delta.npz').resolve():
        parser.error('write a proposal outside the registered certificate')
    pool = dict(np.load(HERE/'candidate_pool.npz',allow_pickle=False))
    before = dict(np.load(HERE/'previous_points.npz',allow_pickle=False))
    raw = np.load(HERE.parent/'partial43/tiny_bodies.npz')['H'].astype(np.int64)
    bi = pool['body_indices'];old = before['body_indices']
    D = int(pool['head_denominator']);d = int(before['head_denominator'])
    assert D%d==0 and all(v.dtype.kind in 'iu' for v in pool.values())
    Z = pool['head_numerators'].astype(object)
    old_Z = before['head_numerators'].astype(object)*(D//d)
    N = 10*D*D+np.sum(Z*Z,axis=1)
    old_N = 10*D*D+np.sum(old_Z*old_Z,axis=1)
    assert 48*int(np.abs(raw).max())**2<2**62
    products = raw[bi]@raw[old].T
    assert not np.any(products%144)
    beta = products//144
    blockers = []
    for i,z in enumerate(Z):
        values = beta[i].astype(object)*D*D+2*old_Z@z
        blockers.append(np.flatnonzero((values>0)&(values*values>N[i]*old_N)))
    old_ids = sorted({int(j) for ids in blockers for j in ids})
    old_map = {j:i for i,j in enumerate(old_ids)}
    count = len(bi);variables = count+len(old_ids)
    rows=[];cols=[];vals=[];upper=[]
    for i,ids in enumerate(blockers):
        for j in ids:
            r=len(upper);rows.extend([r,r]);cols.extend([i,count+old_map[int(j)]])
            vals.extend([1,-1]);upper.append(0)
    products = raw[bi]@raw[bi].T
    assert not np.any(products%144)
    values = (products//144).astype(object)*D*D+2*Z@Z.T
    conflict = (values>0)&(values*values>N[:,None]*N)
    for i,j in zip(*np.where(np.triu(conflict,1))):
        r=len(upper);rows.extend([r,r]);cols.extend([int(i),int(j)])
        vals.extend([1,1]);upper.append(1)
    A=coo_matrix((vals,(rows,cols)),shape=(len(upper),variables)).tocsc()
    # First maximize the gain; among ties use fewer deletions/insertions.
    objective=np.r_[-np.full(count,variables),np.full(len(old_ids),variables+1)]
    result=milp(objective,integrality=np.ones(variables),bounds=Bounds(0,1),
                constraints=LinearConstraint(A,-np.inf,upper),
                options={'time_limit':120,'mip_rel_gap':0})
    if result.x is None:
        raise RuntimeError(result.message)
    x=np.rint(result.x).astype(int)
    assert np.all((x==0)|(x==1)) and np.all(A@x<=upper)
    added=np.flatnonzero(x[:count]);deleted=np.array(old_ids)[x[count:]>0]
    np.savez_compressed(args.output,deleted_indices=deleted,
        added_body_indices=bi[added],added_head_numerators=Z[added].astype(np.int64),
        head_denominator=np.array(D,dtype=np.int64))
    print(json.dumps(dict(candidate_points=count,inserted_points=len(added),
                         deleted_points=len(deleted),gain=len(added)-len(deleted),
                         status=result.message,output=str(args.output)),indent=2))


if __name__ == '__main__':
    main()
