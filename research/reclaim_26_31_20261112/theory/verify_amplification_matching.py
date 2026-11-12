"""Read-only, solver-free replay of all six amplification obstructions.

Each candidate is matched injectively to a forbidden incumbent point.
All matched inequalities are recomputed directly with integer arithmetic.
"""
import argparse
import itertools
import json
import numpy as np
from source_utils import ROOT,add_source_arguments,load_integer_source,radical_positive

def expected_candidates(dimension,source):
    delta=np.load(ROOT/'sextet28_delta.npz',allow_pickle=False)
    added=delta['added'].astype(np.int16)
    assert added.shape==(56,28)
    assert np.all(np.sum(added.astype(np.int64)**2,axis=1)==128)
    heads=np.unique(np.r_[added[:,:24],-added[:,:24]],axis=0)
    assert len(heads)==36
    k=dimension-24
    tails=[np.zeros(k,dtype=np.int16)]
    for i,j in itertools.combinations(range(k),2):
        for a,b in itertools.product((-4,4),repeat=2):
            tail=np.zeros(k,dtype=np.int16);tail[i]=a;tail[j]=b;tails.append(tail)
    old=set(row.astype(np.int16).tobytes()for row in source)
    candidates=[]
    for head in heads:
        for tail in tails:
            row=np.r_[head,tail].astype(np.int16)
            if int(row.astype(np.int64)@row.astype(np.int64))!=128:continue
            if row.tobytes()not in old:candidates.append(row)
    return np.unique(candidates,axis=0)

def verify(dimension,args):
    filename,record,X,data=load_integer_source(dimension,args)
    Y=np.load(ROOT/f'amplify{dimension}_pool.npz',allow_pickle=False)['candidates'].astype(np.int64)
    partners=np.load(ROOT/f'amplify{dimension}_matching.npz',allow_pickle=False)['incumbent_partner'].astype(np.int64)
    expected=expected_candidates(dimension,X)
    assert np.array_equal(Y,expected),'Candidate pool does not match the advertised complete tail-image family'
    assert partners.shape==(len(Y),) and np.all(partners>=0)
    assert len(np.unique(partners))==len(Y)
    assert np.all(partners<record['points'])
    assert np.all(np.sum(Y*Y,axis=1)==128)
    if dimension!=26:
        products=np.sum(Y*X[partners],axis=1)
        assert np.all(products>64),'A matched pair is not forbidden'
        quadratic_count=0
    else:
        qa=data['quadratic_A'].astype(np.int64);qb=data['quadratic_B'].astype(np.int64)
        den=data['quadratic_denominators'].astype(np.int64);rad=data['quadratic_radicands'].astype(np.int64)
        assert len(qa)==16 and np.all(den>0) and np.all(rad>0)
        assert np.all(np.sum(qa*qa,1)+rad*np.sum(qb*qb,1)==32*den*den)
        assert np.all(np.sum(qa*qb,1)==0)
        assert np.all(np.any(qb!=0,axis=1))
        quadratic_count=0
        for row,partner in zip(Y,partners):
            if partner<len(X):assert int(row@X[partner])>64
            else:
                j=int(partner-len(X));a=int(row@qa[j])-32*int(den[j]);b=int(row@qb[j])
                assert radical_positive(a,b,int(rad[j])),'Quadratic matched pair is not forbidden'
                quadratic_count+=1
    return {'dimension':dimension,'source_filename':filename.name,'source_sha256':record['sha256'],'incumbent_points':record['points'],'candidate_points':len(Y),'distinct_forbidden_partners':len(partners),'quadratic_matched_pairs':quadratic_count,'complete_candidate_regeneration':True,'maximum_gain':0,'scope':'Every subset of this supplied complete candidate family requires at least as many incumbent deletions as insertions. This is not a bound on other families.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimensions',nargs='+',type=int,default=list(range(26,32)))
    add_source_arguments(parser);args=parser.parse_args()
    assert all(26<=n<=31 for n in args.dimensions)
    print(json.dumps([verify(n,args)for n in args.dimensions],indent=2))

if __name__=='__main__':main()
