"""Search nonintegral rational tail rotations absent from the integer image pool."""
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
from discovery.engines.selection import select_coupled_rows


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('dimension',type=int,choices=(27,29,30,31))
    parser.add_argument('--trials',type=int,default=4)
    parser.add_argument('--seed',type=int,default=20261104)
    parser.add_argument('--mixed',action='store_true',help='Mix one core coordinate and rotate the complete auxiliary subset')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    start=monotonic();d=args.dimension
    p=next((ROOT/'verifications/improved/dim27-31-algebraic-completion/data').glob(f'coordinates{d}_*.npz'),None)
    if p is None:p=next((ROOT/'verifications/improved/dim27-31-joint-refills/data').glob(f'coordinates{d}_*.npz'))
    X=np.load(p,allow_pickle=False)['points'].astype(np.int64)
    assert np.all(np.sum(X*X,axis=1)==128)
    rows=X.tolist();rng=np.random.default_rng(args.seed)
    image_rows=X[np.any(X[:,24:]!=0,axis=1)].tolist() if args.mixed else rows
    grammar=[(2,1,1),(2,2,1,1),(3,1,1,1),(1,1,1,1,1,1)]
    grammar=[shape for shape in grammar if len(shape)<=d-24]
    context=RunContext(args.out,{'max_pairs':150000000000,'max_gram_entries':50000000,
        'max_materialized_points':600000,'max_rank_work':10000000000})
    report={'dimension':d,'baseline':len(X),'mixed_core_tail':args.mixed,'image_source_points':len(image_rows),
        'source':str(p.relative_to(ROOT)),
        'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'trials':[],
        'inherited_assumptions':'Source geometry inherited; every matrix checked as an exact isometry; full output proof mandatory for a gain',
        'scope':'Exact maximum in each finite source/image union, without any blocker cutoff; no global or priority claim'}
    for trial in range(args.trials):
        shape=grammar[trial%len(grammar)]
        normal=[0]*d
        indices=(np.r_[rng.integers(24),rng.choice(np.arange(24,d),len(shape)-1,replace=False)]
                 if args.mixed else rng.choice(np.arange(24,d),len(shape),replace=False))
        for i,x in zip(indices,shape):
            normal[int(i)]=int(x*rng.choice([-1,1]))
        data,receipt=select_coupled_rows(rows,image_rows,[1]*d,128,{'normals':[normal]},context)
        keep=('source_cardinality','image_new_points','overlap_points','cross_pairs_checked',
              'common_left_points','conflict_edges','graph_sha256','integer_backend',
              'matrix_denominator','coordinate_gcd_reduction','optimal_in_supplied_graph')
        record={k:receipt[k] for k in keep}
        record.update(trial=trial,normal=normal,result=len(data['coordinates']),
                      gain=len(data['coordinates'])-len(X),elapsed_seconds=monotonic()-start,
                      matching_size=len(receipt['matching']))
        report['trials'].append(record)
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(record),flush=True)
        if record['gain']>0:
            Y=np.array(data['coordinates'],np.int64)
            known={tuple(row) for row in (X*receipt['matrix_denominator']//max(1,receipt['coordinate_gcd_reduction'])).tolist()}
            # Independent changed-pair products, using direct int64 rather than
            # the bounded binary32 cross-graph arithmetic.
            added=[i for i,row in enumerate(Y.tolist()) if tuple(row) not in known]
            for at in range(0,len(added),32):
                dots=Y[added[at:at+32]]@Y.T
                for j,i in enumerate(added[at:at+32]):dots[j,i]=-data['norm_squared']
                assert int(dots.max(initial=-1))<=data['norm_squared']//2
            np.savez_compressed(args.out/'candidate.npz',points=Y,
                squared_norm=data['norm_squared'],inner_product_limit=data['norm_squared']//2,
                coordinate_denominator=2*receipt['matrix_denominator']//max(1,receipt['coordinate_gcd_reduction']))
            (args.out/'candidate_receipt.json').write_text(json.dumps(receipt)+'\n')
            break

if __name__=='__main__':main()
