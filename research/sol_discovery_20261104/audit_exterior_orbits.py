"""Exact E8 exterior-square controls and the reflection-stabilizer exception.

All seeds, roots and generators are constructed here. No result or source
file from another executor is assumed. Controls are far below 28D221988.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations,product
from pathlib import Path
import json
import argparse
import sys
from time import monotonic

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from discovery.core import RunContext,Artifact,digest
from discovery.certification import verify_integer_spherical
from discovery.exact_gram import rational_json,exact_integer_rank

PAIRS=list(combinations(range(8),2))


def require(value,message):
    if not value:raise ValueError(message)


def roots():
    result=[]
    for i,j in PAIRS:
        for a,b in product((-2,2),repeat=2):
            row=np.zeros(8,np.int64);row[i]=a;row[j]=b;result.append(row)
    for signs in product((-1,1),repeat=8):
        if signs.count(-1)%2==0:result.append(np.array(signs,np.int64))
    require(len(result)==240 and all(int(r@r)==8 for r in result),'E8 root generation')
    return result


def simple_roots():
    result=[np.array([1,-1,-1,-1,-1,-1,-1,1],np.int64),np.array([2,2,0,0,0,0,0,0],np.int64)]
    for j in range(6):
        row=np.zeros(8,np.int64);row[j]=-2;row[j+1]=2;result.append(row)
    return result


def row(matrix):return tuple(int(matrix[i,j]) for i,j in PAIRS)


def exact_orbit(seed,generators,cap=20000):
    matrices=[seed];seen={row(seed)}
    for matrix in matrices:
        for root in generators:
            image=matrix.copy();ar=matrix@root
            delta=np.outer(root,ar)-np.outer(ar,root)
            require(np.all(delta%4==0),'Integer orbit scale')
            image+=delta//4
            key=row(image)
            if key not in seen:
                require(len(matrices)<cap,'Exact orbit cap reached')
                seen.add(key);matrices.append(image)
    return matrices


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    out=parser.parse_args().out;out.mkdir(parents=True,exist_ok=False)
    start=monotonic();R=roots();S=simple_roots()
    # Closure of the finite spanning root set proves finiteness of the actual
    # generated orthogonal matrix group, without enumerating all group elements.
    root_keys=set(map(tuple,R))
    for s in S:
        require(int(s@s)==8,'Simple-root norm')
        for r in R:
            numerator=4*r-int(r@s)*s
            require(np.all(numerator%4==0) and tuple(numerator//4) in root_keys,'Reflection preserves complete E8 root set')
    rank,_=exact_integer_rank([r.tolist() for r in R],RunContext(ROOT))
    require(rank==8,'Root span')
    # Verify the entire quadratic class average directly, with common
    # denominator16 on the exterior-square matrices. No irreducibility premise
    # is needed by this explicit arithmetic replay.
    lines=[r for r in R if next(int(x) for x in r if x)>0]
    require(len(lines)==120,'Root reflection count')
    class_matrices=[]
    for r in lines:
        M=4*np.eye(8,dtype=np.int64)-np.outer(r,r)
        require(np.array_equal(M@M.T,16*np.eye(8,dtype=np.int64)),'Exact reflection matrix')
        T=np.array([[M[i,k]*M[j,l]-M[i,l]*M[j,k] for k,l in PAIRS] for i,j in PAIRS],np.int64)
        require(np.array_equal(T@T.T,256*np.eye(28,dtype=np.int64)),'Exact exterior isometry')
        class_matrices.append(T)
    require(np.array_equal(sum(class_matrices),960*np.eye(28,dtype=np.int64)),'Reflection class average is exactly one half')
    # The outer products of E8 roots span Sym8, establishing the scalar
    # quadratic-form consequence only when every reflection inequality is tight.
    outer_rows=[[int(r[i]*r[j]) for i in range(8) for j in range(i,8)] for r in lines]
    span,_=exact_integer_rank(outer_rows,RunContext(ROOT))
    require(span==36,'Root outer products span Sym8')
    results=[]
    for label,blocks,expected in [('rank_two',[(0,1)],7560),('complex_structure',[(0,1),(2,3),(4,5),(6,7)],15120)]:
        seed=np.zeros((8,8),np.int64)
        for i,j in blocks:seed[i,j]=4;seed[j,i]=-4
        v=np.array(row(seed),np.int64);norm=int(v@v)
        fixed=0;class_dots=Counter()
        for T in class_matrices:
            image=T@v
            if np.array_equal(image,16*v):fixed+=1
            dot=Fraction(int(v@image),16*norm);class_dots[str(dot)]+=1
            require(np.array_equal(image,16*v) or dot<=Fraction(1,2),'Reflection compatibility')
        matrices=exact_orbit(seed,S)
        require(len(matrices)==expected,'Independent exact orbit cardinality')
        coordinates=[list(row(A)) for A in matrices]
        ctx=RunContext(ROOT,{'max_pairs':120000000})
        proof=verify_integer_spherical({'coordinates':coordinates,'norm_squared':norm},ctx)
        require(proof.details['facts']['rank_exact']==28,'Exterior orbit span')
        G=seed.T@seed;scalar=bool(np.array_equal(G,G[0,0]*np.eye(8,dtype=np.int64)))
        np.savez_compressed(out/f'{label}.npz',points=np.array(coordinates,np.int16),norm_squared=norm)
        result={'label':label,'cardinality':len(coordinates),'rank_exact':28,'norm_squared':norm,
                'reflection_fixed_actions':fixed,'reflection_normalized_dot_histogram':dict(class_dots),
                'AtA_scalar':scalar,'AtA_diagonal':np.diag(G).tolist(),
                'complete_exact_coordinate_proof':proof.to_dict(),
                'scope':'Exactly closed whole orbit under the generated E8 reflections; known-scale control, no record or novelty claim'}
        results.append(result)
        print(json.dumps({k:result[k] for k in ['label','cardinality','reflection_fixed_actions','AtA_scalar']}),flush=True)
    report={'status':'PASS','seconds':monotonic()-start,'reflection_class_size':120,
            'exact_quadratic_class_average':[1,2],'root_outer_product_span_rank':36,'controls':results,
            'correction':'A whole-orbit kissing seed need not be a complex structure when reflections stabilize it. Rank-two seed has a valid7560-point whole orbit,30 fixing reflections and non-scalar AtA.',
            'free_reflection_sector':'With no fixing reflections, equality forces every reflection dot to norm/2; root outer-product spanning then forces AtA scalar.',
            'comparison_28D':221988,'new_sota':False}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
