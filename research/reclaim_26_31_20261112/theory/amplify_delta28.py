"""Regenerate complete transplant pools and exact no-gain matchings.

Writes only to the explicitly requested output directory. The recorded MILP
experiments are superseded by the stronger solver-independent matchings.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from source_utils import add_source_arguments,load_integer_source,radical_positive
from verify_amplification_matching import expected_candidates

def run(dimension,args):
    started=time.time();filename,record,X,source=load_integer_source(dimension,args)
    Y=expected_candidates(dimension,X)
    assert dimension*int(np.abs(X).max())*int(np.abs(Y).max())<2**24
    ptr=[0];ids=[];matrix=X.T.astype(np.float32)
    if dimension==26:
        qa=source['quadratic_A'].astype(np.int64);qb=source['quadratic_B'].astype(np.int64)
        den=source['quadratic_denominators'];rad=source['quadratic_radicands']
    for start in range(0,len(Y),64):
        products=Y[start:start+64].astype(np.float32)@matrix
        for offset,row in enumerate(products):
            blockers=np.flatnonzero(row>64).tolist()
            if dimension==26:
                y=Y[start+offset].astype(np.int64)
                for j,(a,b,d,r)in enumerate(zip(qa,qb,den,rad)):
                    if radical_positive(int(y@a)-32*int(d),int(y@b),int(r)):blockers.append(len(X)+j)
            ids.extend(blockers);ptr.append(len(ids))
    ids=np.array(ids,dtype=np.int32);ptr=np.array(ptr,dtype=np.int64)
    columns,inverse=np.unique(ids,return_inverse=True)
    graph=csr_matrix((np.ones(len(ids),dtype=np.int8),inverse,ptr),shape=(len(Y),len(columns)))
    matching=maximum_bipartite_matching(graph,perm_type='column')
    partner=np.full(len(Y),-1,dtype=np.int64);paired=matching>=0
    partner[paired]=columns[matching[paired]]
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/f'amplify{dimension}_pool.npz',candidates=Y,blocker_ptr=ptr,blocker_ids=ids)
    np.savez_compressed(out/f'amplify{dimension}_matching.npz',incumbent_partner=partner)
    report={'dimension':dimension,'source':filename.name,'source_sha256':record['sha256'],'source_count':record['points'],'candidate_count':len(Y),'matched':int(paired.sum()),'all_candidates_matched':bool(np.all(paired)),'minimum_blockers':int(np.diff(ptr).min()),'maximum_blockers':int(np.diff(ptr).max()),'blocker_incidences':len(ids),'distinct_blockers':len(columns),'elapsed_seconds':time.time()-started}
    print(json.dumps(report),flush=True)
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__);add_source_arguments(parser)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--dimensions',nargs='+',type=int,default=list(range(26,32)))
    args=parser.parse_args();assert all(26<=n<=31 for n in args.dimensions)
    reports=[run(n,args)for n in args.dimensions]
    (args.output_dir/'amplification_exact_matching_report.json').write_text(json.dumps(reports,indent=2)+'\n')

if __name__=='__main__':main()
