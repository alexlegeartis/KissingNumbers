"""Exact zero-gain certificate for the 3960 supplied residual-axis candidates.

3928 candidates use twelve axis tails with32 incumbent points each. For any
fixed tail, the core lies inR16, has norm64, and needs nonpositive pair products,
so at most32 core points can remain: added_axis <= deleted_axis.
The remaining32 candidates have an injective set of incompatible incumbent
blockers outside those384 axis-fiber points: added_other <= deleted_other.
The disjoint charges prove totaladded <= totaldeleted for the entire pool.
There is no need for an optimizer or an incomplete numerical solver bound.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

HERE=Path(__file__).resolve().parent
NORMAL=np.array([0,2,4,5,6,10,11,23])
COMPLEMENT=np.array([j for j in range(24) if j not in NORMAL])
SOURCE_SHA256='714c9a0189cd1ff2f9363f6672e84506606ca54cb806a18d155d817481019943'
POOL_SHA256='e6f0e5c46a93f527d1d3b4cb10cfb82a0196ee254b6f8aa3df321fbf6bfcd0b9'


def require(value,label):
    if not value:raise ValueError(label)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('source',type=Path)
    ap.add_argument('--pool',type=Path,default=HERE/'pool.npz')
    ap.add_argument('--report',type=Path,default=HERE/'finite_pool_report.json')
    ap.add_argument('--certificate',type=Path,default=HERE/'finite_pool_certificate.npz')
    ap.add_argument('--verify-only',action='store_true')
    args=ap.parse_args()
    require(hashlib.sha256(args.source.read_bytes()).hexdigest()==SOURCE_SHA256,'Pinned source hash')
    require(hashlib.sha256(args.pool.read_bytes()).hexdigest()==POOL_SHA256,'Pinned candidate pool hash')
    with np.load(args.source,allow_pickle=False) as z:
        X=z['points']
        require(X.dtype.kind in 'iu' and X.shape==(251254,31),'Integer source format')
        require(int(z['squared_norm'])==128,'Source normalization')
        X=X.astype(np.int64)
    with np.load(args.pool,allow_pickle=False) as z:
        C=z['candidates'];ptr=z['blocker_indptr'];blockers=z['blocker_indices']
        require(C.dtype.kind in 'iu' and C.shape==(3960,31),'Integer pool format')
        C=C.astype(np.int64)
    require(np.all(np.sum(X*X,axis=1)==128),'All source norms')
    require(np.all(np.sum(C*C,axis=1)==128),'All candidate norms')
    require(len(np.unique(C,axis=0))==3960,'Candidate distinctness')
    require(not np.any(C[:,NORMAL]),'Candidates lie in complementary core space')
    require(np.all(np.sum(C[:,COMPLEMENT]**2,axis=1)==64),'Core norm64')
    require(np.all(np.sum(C[:,24:]**2,axis=1)==64),'Tail norm64')
    axis=np.count_nonzero(C[:,24:],axis=1)==1
    tails=np.unique(C[axis,24:],axis=0)
    require(axis.sum()==3928 and tails.shape==(12,7),'Axis candidate partition')
    require(np.all(np.max(np.abs(tails),axis=1)==8),'Actual axis tails')
    axis_old=[]
    group_rows=[]
    for t in tails:
        ids=np.flatnonzero(np.all(X[:,24:]==t,axis=1)&~np.any(X[:,NORMAL],axis=1))
        require(len(ids)==32,'Saturated incumbent axis fiber')
        H=X[ids][:,COMPLEMENT]
        require(np.all(np.sum(H*H,axis=1)==64),'Incumbent fiber core norms')
        D=H@H.T;np.fill_diagonal(D,-64)
        require(D.max()<=0 and len(np.unique(H,axis=0))==32,'Incumbent orthoplex geometry')
        axis_old.extend(map(int,ids));group_rows.append(ids)
    axis_old=np.array(sorted(axis_old),np.int64)
    require(len(np.unique(axis_old))==384,'Disjoint axis victim groups')
    exception=np.flatnonzero(~axis)
    require(len(exception)==32,'Remaining candidate count')
    protected=set(map(int,axis_old))
    if args.verify_only:
        with np.load(args.certificate,allow_pickle=False) as z:
            matched=z['matched_old_indices']
            require(np.array_equal(z['exception_candidate_indices'],exception),'Certificate candidate coverage')
            require(np.array_equal(z['axis_old_indices'],axis_old),'Certificate axis groups')
    else:
        p=[0];b=[]
        for i in exception:
            b.extend(int(j) for j in blockers[ptr[i]:ptr[i+1]] if int(j) not in protected)
            p.append(len(b))
        graph=csr_matrix((np.ones(len(b),np.int8),b,p),shape=(32,len(X)))
        matched=maximum_bipartite_matching(graph,perm_type='column')
    require(matched.shape==(32,) and np.all((matched>=0)&(matched<len(X))),'Complete blocker assignment')
    require(len(np.unique(matched))==32,'Injective outside blockers')
    require(not any(int(j) in protected for j in matched),'Separate deletion charges')
    products=np.sum(C[exception]*X[matched],axis=1)
    require(np.all(products>64),'Exact forbidden products for every matched pair')
    report=dict(status='PASS',source_sha256=SOURCE_SHA256,pool_sha256=POOL_SHA256,
                source_count=251254,candidate_count=3960,axis_candidates=3928,
                saturated_axis_fibers=12,old_points_per_axis_fiber=32,
                disjoint_axis_victims=384,exception_candidates=32,
                injectively_matched_outside_victims=32,
                minimum_matched_product=int(products.min()),inner_product_limit=64,
                exact_gain_upper_bound=0,finite_union_population_upper_bound=251254,
                scope='Every feasible subset of this pinned source union this3960-point pool. No exclusion of other candidates, radii, deformations, or kissing configurations.')
    if not args.verify_only:
        np.savez_compressed(args.certificate,exception_candidate_indices=exception,
                            matched_old_indices=matched,axis_old_indices=axis_old,
                            axis_tail_vectors=tails)
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
