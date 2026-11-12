"""Exactly price all whole 128-sign three-core/five-tail fibers."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from source_utils import add_source_arguments,load_integer_source

def main():
    parser=argparse.ArgumentParser(description=__doc__);add_source_arguments(parser);parser.add_argument('--output-dir',type=Path,required=True);args=parser.parse_args()
    filename,record,X,_=load_integer_source(31,args)
    profiles,weights=np.unique(np.abs(X),axis=0,return_counts=True)
    assert np.max(np.count_nonzero(profiles[:,24:],axis=1))<=4
    moving=np.any(profiles[:,24:]!=0,axis=1)
    assert np.all(np.sort(profiles[~moving,:24],axis=1)[:,-3:].sum(1)<=16)
    profiles=profiles[moving];weights=weights[moving]
    core=np.array(list(itertools.combinations(range(24),3)),dtype=np.int16)
    tail=np.array(list(itertools.combinations(range(7),5)),dtype=np.int16)
    C=profiles[:,core].sum(2).T;D=profiles[:,24+tail].sum(2).T
    cost=np.array([(C+d>16)@weights for d in D],dtype=np.int64).T
    assert cost.shape==(2024,21)
    normal={0,2,4,5,6,10,11,23}
    inside=np.array([set(map(int,b))<=normal for b in core]);outside=np.array([not(set(map(int,b))&normal)for b in core])
    summary={}
    for label,mask in [('inside_normal_octad',inside),('outside_normal_octad',outside),('mixed',~(inside|outside))]:
        summary[label]={'fibers':int(mask.sum())*21,'minimum_blockers':int(cost[mask].min()),'median_blockers':float(np.median(cost[mask]))}
    report={'source':filename.name,'source_sha256':record['sha256'],'fibers':int(cost.size),'points_per_fiber':128,'minimum_blockers':int(cost.min()),'maximum_blockers':int(cost.max()),'groups':summary}
    args.output_dir.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output_dir/'full_3plus5_bundle_costs.npz',core_triples=core,tail_fives=tail,blocker_counts=cost)
    (args.output_dir/'full_3plus5_bundle_costs.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
