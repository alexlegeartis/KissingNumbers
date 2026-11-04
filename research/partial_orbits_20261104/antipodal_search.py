"""Exact coupled selection from a configuration and its full antipodal image."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from time import monotonic

import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.selection import conflict_graph, maximum_independent_bipartite


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('dimension',type=int,choices=(27,29,30,31))
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    start=monotonic();d=args.dimension
    p=next((ROOT/'verifications/improved/dim27-31-algebraic-completion/data').glob(f'coordinates{d}_*.npz'),None)
    if p is None:p=next((ROOT/'verifications/improved/dim27-31-joint-refills/data').glob(f'coordinates{d}_*.npz'))
    X=np.load(p,allow_pickle=False)['points'].astype(np.int64)
    known={tuple(row) for row in X.tolist()}
    Y=np.array([-row for row in X if tuple(-row) not in known],np.int64).reshape(-1,d)
    assert np.all(np.sum(X*X,axis=1)==128) and len(known)==len(X)
    context=RunContext(args.out,{'max_pairs':2000000000,'max_gram_entries':30000000,
                               'max_materialized_points':1000000})
    blockers,backend=conflict_graph(Y.tolist(),X.tolist(),[1]*d,128,context)
    certificate=maximum_independent_bipartite(blockers,len(X))
    keep=certificate['right_indices'];added=certificate['left_indices']
    final=np.r_[X[keep],Y[added]]
    gain=len(final)-len(X)
    # Independent direct-int64 changed pairs; unlike the fast screening path,
    # these products use a different backend and no matching assumptions.
    for at in range(0,len(added),32):
        dots=Y[added[at:at+32]]@final.T
        for j in range(len(dots)):dots[j,len(keep)+at+j]=-128
        assert int(dots.max(initial=-128))<=64
    report={'dimension':d,'baseline':len(X),'result':len(final),'gain':gain,
        'source':str(p.relative_to(ROOT)),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
        'antipodal_candidates':len(Y),'screening_pairs':len(Y)*len(X),
        'conflicts':sum(map(len,blockers)),'backend':backend,
        'matching_size':len(certificate['matching']),
        'selected_new':len(added),'removed_old':len(X)-len(keep),
        'exact_finite_optimality':'Matching and equal-size vertex cover checked against every exact conflict',
        'changed_pairs':'PASS direct int64; common norm/distinctness verified',
        'changed_pairs_count':len(added)*(len(final)-1),
        'inherited_assumptions':'Source validity inherited for search; fresh full proof required for any gain',
        'elapsed_seconds':monotonic()-start,'new_record':False,
        'scope':'Maximum cardinality of this source/antipode union only; no global bound'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (args.out/'matching.json').write_text(json.dumps({k:certificate[k] for k in
        ('matching','cover_left','cover_right','optimal_in_supplied_graph')})+'\n')
    if gain>0:np.savez_compressed(args.out/'candidate.npz',points=final.astype(np.int16),squared_norm=128,
                                 inner_product_limit=64,coordinate_denominator=2)
    print(json.dumps(report),flush=True)

if __name__=='__main__':main()
