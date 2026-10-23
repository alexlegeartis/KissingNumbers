"""Exact moment identities obstruct every codimension 1--3 linear section.

NumPy integer arithmetic is used for the large moment contractions.
All matrix inversion and final inequalities use exact rational arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
import itertools
import json
import math
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from exact_lattice import enumerate_short


def inverse(A):
    n=len(A);B=[[F(v) for v in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if B[i][j])
        B[j],B[pivot]=B[pivot],B[j]
        value=B[j][j];B[j]=[v/value for v in B[j]]
        for i in range(n):
            if i==j:continue
            value=B[i][j];B[i]=[a-value*b for a,b in zip(B[i],B[j])]
    return [row[n:] for row in B]


def main():
    start=time.monotonic()
    cert=json.loads((HERE/'certificate.json').read_text())
    G=np.array(cert['gram_12'],dtype=np.int64);M=np.array(cert['sublattice_matrix'],dtype=np.int64)
    inv=inverse(G.tolist())
    H=np.array([[int(3*x) for x in row] for row in inv],dtype=np.int64)
    assert np.array_equal(H@G,3*np.eye(12,dtype=np.int64))
    short=list(enumerate_short(G.tolist()))
    X=np.array([v for v,n in short],dtype=np.int64)
    norms=np.array([n for v,n in short],dtype=np.int64)
    first,second=cert['label_first'],cert['label_second'];p=3**np.arange(6)
    labels=((X[:,first]+X[:,second])%3)@p
    weights=np.zeros(729,dtype=np.int64);mult=np.ones(729,dtype=np.int64)
    S=np.zeros((729,12,12),dtype=np.int64)
    for c in range(1,729):
        V=X[labels==c];assert len(V)>0
        assert np.array_equal(V.sum(axis=0),np.zeros(12,dtype=np.int64))
        weights[c]=norms[labels==c][0];mult[c]=len(V);S[c]=V.T@V
    base=np.array(cert['base_outer_cosets'],dtype=np.int64)
    digits=(np.arange(729)[:,None]//p)%3
    negative=((-digits)%3)@p
    assert set(map(tuple,negative[base]))==set(map(tuple,base))
    for j in range(3):
        other=[k for k in range(3) if k!=j]
        marginal=np.zeros(729,dtype=np.int64)
        np.add.at(marginal,base[1:,j],np.prod(mult[base[1:][:,other]],axis=1))
        assert np.all(marginal[weights==4]==81)
        assert np.all(marginal[weights==6]==36)
        assert np.all(marginal[weights==8]==3)
    # Nonzero marginal vectors of the fixed 403704-point part in one block.
    rootL=X[norms==4]@M
    V=np.concatenate([rootL,X]);coef=np.concatenate([np.ones(756,dtype=np.int64),np.choose(norms//2,[0,0,81,36,3])])
    moment2=V.T@(coef[:,None]*V)
    assert np.array_equal(moment2,44856*H)
    checked4=0
    for i,j,k,l in itertools.combinations_with_replacement(range(12),4):
        value=int(coef@(V[:,i]*V[:,j]*V[:,k]*V[:,l]))
        expected=6768*int(H[i,j]*H[k,l]+H[i,k]*H[j,l]+H[i,l]*H[j,k])
        assert value==expected,(i,j,k,l,value,expected)
        checked4+=1
    # The paired moments factor exactly. Index pairs are ordered here;
    # the 144x144 identity verifies every coefficient, including repeats.
    mixed=[]
    for j in range(3):
        k=(j+1)%3
        rows=base[(base[:,j]!=0)&(base[:,k]!=0)]
        A=S[rows[:,j]].reshape(len(rows),144)
        B=S[rows[:,k]].reshape(len(rows),144)
        observed=A.T@B
        expected=3528*np.outer(H.ravel(),H.ravel())
        assert np.array_equal(observed,expected)
        mixed.append({'blocks':[j,k],'cells':len(rows),'integer_matrix_entries_checked':144**2})
    # In orthonormal coordinates, for any normal u=(u1,u2,u3),
    # sum (u.x)^2 = 134568 |u|^2 and
    # sum (u.x)^4 = 95256 |u|^4 + 87480 sum_j |uj|^4.
    # For a rank-r orthogonal projector P, contraction gives
    # sum (x^T P x)^2 <= 60912*r*(r+2), while
    # sum (x^T P x)=134568*r.
    # Cauchy--Schwarz therefore forces at least ceil(ratio) nonzero
    # projections among the fixed base points. Add every other point
    # to the section, which is an upper bound regardless of its geometry.
    table=[]
    benchmarks={33:370240,34:390980,35:415308}
    for r in [1,2,3]:
        ratio=F((134568*r)**2,60912*r*(r+2))
        outside=math.ceil(ratio)
        table.append({'dimension':36-r,'codimension':r,
                      'base_points_forced_outside_at_least':outside,
                      'cauchy_schwarz_ratio':str(ratio),
                      'section_points_at_most':cert['points_36']-outside,
                      'cohn_snapshot_lower_bound':benchmarks[36-r]})
    result={'all_checks_pass':True,'base_points':403704,'full_code_points':cert['points_36'],
            'second_moment_coefficient':134568,
            'fourth_moment_global_coefficient':95256,
            'fourth_moment_block_coefficient':87480,
            'symmetric_fourth_moment_entries_checked':checked4,
            'mixed_moment_checks':mixed,'section_upper_bounds':table,
            'affine_section_not_through_origin_points_at_most':403704//2+(cert['points_36']-403704),
            'scope':'sections of this particular 36D configuration; not universal kissing-number upper bounds. Noncentral affine sections obey the antipodal-base bound.',
            'arithmetic':'exact rational inverse; NumPy int64 moment sums, with magnitudes below 2^63',
            'elapsed_seconds':time.monotonic()-start}
    (HERE/'section_obstruction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
