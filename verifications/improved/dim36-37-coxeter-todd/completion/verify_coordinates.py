"""Independent integer-coordinate replay; no imports from search or verifier.

Exact LDL enumeration reconstructs the quotient. Matrix products of expanded
points check every phase/phase and latitude/latitude pair not settled by
different-label lattice distances, and every exceptional latitude/equator pair.
"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import hashlib
import itertools
import json
import math
import time
import numpy as np

HERE=Path(__file__).resolve().parent


def short_vectors(gram):
    n=len(gram)
    lower=[[Fraction(int(i==j)) for j in range(n)] for i in range(n)]
    diagonal=[]
    for i in range(n):
        d=Fraction(gram[i][i])-sum(lower[i][k]**2*diagonal[k] for k in range(i))
        assert d>0
        diagonal.append(d)
        for j in range(i+1,n):
            lower[j][i]=(gram[j][i]-sum(lower[j][k]*lower[i][k]*diagonal[k]
                                      for k in range(i)))/d
    coefficient=[0]*n
    def enumerate_level(level,spent):
        if level<0:
            if spent:yield tuple(coefficient),int(spent)
            return
        offset=sum((lower[j][level]*coefficient[j] for j in range(level+1,n)),Fraction(0))
        available=(8-spent)/diagonal[level]
        radius=math.isqrt(available.numerator//available.denominator)+2
        center=math.floor(-offset)
        for x in range(center-radius,center+radius+1):
            cost=diagonal[level]*(x+offset)**2
            if spent+cost<=8:
                coefficient[level]=x
                yield from enumerate_level(level-1,spent+cost)
    yield from enumerate_level(n-1,Fraction(0))


def exact_det(matrix):
    a=[[Fraction(x) for x in row] for row in matrix]
    determinant=Fraction(1)
    for k in range(len(a)):
        pivot=next(j for j in range(k,len(a)) if a[j][k])
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];determinant=-determinant
        value=a[k][k];determinant*=value
        for j in range(k+1,len(a)):
            ratio=a[j][k]/value
            for col in range(k+1,len(a)):a[j][col]-=ratio*a[k][col]
    return determinant


def main():
    started=time.monotonic();path=HERE/'certificate.json'
    c=json.loads(path.read_text())
    assert c['format']=='k12_phase_and_lift_orbits_v2' and c['squared_norm']==12
    G=np.array(c['gram_12'],dtype=np.int64);M=np.array(c['sublattice_matrix'],dtype=np.int64)
    assert G.shape==M.shape==(12,12) and np.array_equal(G,G.T)
    assert np.array_equal(M@G@M.T,3*G) and abs(exact_det(M.tolist()))==729
    phi=np.eye(12,dtype=np.int64)-M
    assert np.array_equal(phi@G@phi.T,G) and np.array_equal(phi@phi@phi,np.eye(12,dtype=np.int64))
    first,second=c['label_first'],c['label_second']
    assert len(first)==len(second)==6 and sorted(first+second)==list(range(12))
    powers=3**np.arange(6,dtype=np.int64)
    def labels(a):return ((a[:,first]+a[:,second])%3)@powers
    assert np.all(labels(M)==0)
    short=list(short_vectors(G.tolist()))
    assert Counter(n for _,n in short)=={4:756,6:4032,8:20412}
    vectors=np.array([v for v,_ in short],dtype=np.int64)
    norms=np.einsum('ij,ij->i',vectors@G,vectors)
    assert np.array_equal(norms,np.array([n for _,n in short]))
    lab=labels(vectors);weights=np.full(729,100,dtype=np.int64);mult=np.zeros(729,dtype=np.int64)
    weights[0],mult[0]=0,1
    for v,n,label in zip(vectors,norms,lab):
        if n<weights[label]:weights[label],mult[label]=n,1
        elif n==weights[label]:mult[label]+=1
    assert Counter(zip(weights.tolist(),mult.tolist()))=={(0,1):1,(4,3):252,(6,18):224,(8,81):252}
    digits=np.arange(729)[:,None]//powers%3
    D=weights[((digits[:,None,:]-digits[None,:,:])%3)@powers]
    roots=np.flatnonzero(weights==4)
    zero=np.array(c['root_phase_zero_vectors'],dtype=np.int64)
    assert zero.shape==(252,12) and np.array_equal(labels(zero),roots)
    phases=np.stack((zero,zero@phi,zero@phi@phi),axis=1)
    for label,row in zip(roots,phases):
        assert set(map(tuple,row))==set(map(tuple,vectors[(lab==label)&(norms==4)]))
    phases=np.concatenate((phases,np.zeros((1,3,12),dtype=np.int64)))
    base=np.array(c['base_outer_cosets'],dtype=np.int64)
    assert base.shape==(1429,3) and np.all(base[0]==0)
    assert len(np.unique(base,axis=0))==len(base)
    assert np.all(np.sum(weights[base[1:]],axis=1)==12)
    assert np.all(np.any(base[1:]==0,axis=1))
    base_dist=sum(D[base[:,k,None],base[None,:,k]] for k in range(3))
    np.fill_diagonal(base_dist,12);assert base_dist.min()>=12
    fixed_count=2268+sum(math.prod(int(mult[x]) for x in row) for row in base[1:])
    assert fixed_count==403704

    def expand(rows,lift=False):
        array=np.array(rows,dtype=np.int64)
        assert array.ndim==2 and array.shape[1]==5
        assert np.all(array>=0) and np.all(array[:,3:]<3)
        assert np.all(array[:,:3]<(253 if lift else 252))
        assert len(np.unique(array,axis=0))==len(array)
        if lift:
            assert np.all(np.sum(array[:,:3]==252,axis=1)==1)
            for a,b,d,p,q in array:
                codes=(a,b,d);angles=(0,p,q)
                nonzero=[k for k in range(3) if codes[k]!=252]
                assert angles[nonzero[0]]==0 and angles[codes.index(252)]==0
        points=[];point_labels=[]
        for a,b,d,p,q in array:
            for t in range(3):
                points.append(np.concatenate((phases[a,t],phases[b,(p+t)%3],phases[d,(q+t)%3])))
                point_labels.append((a,b,d))
        points=np.array(points,dtype=np.int64);point_labels=np.array(point_labels,dtype=np.int64)
        weighted=(points.reshape(-1,3,12)@G).reshape(-1,36)
        # A conservative explicit bound rules out overflow in all dot products.
        assert int(np.abs(points).max())*int(np.abs(weighted).max())*36<2**63
        assert np.all(np.sum(points*weighted,axis=1)==(8 if lift else 12))
        assert len(np.unique(points,axis=0))==len(points)
        cells=np.unique(point_labels,axis=0)
        for cell in cells:
            outer=np.array([0 if k==252 else roots[k] for k in cell])
            distances=sum(D[outer[k],base[:,k]] for k in range(3))
            assert distances.min()>=(8 if lift else 12)
        return points,weighted,point_labels

    equator,weighted,eq_labels=expand(c['phase_orbits'])
    products=0
    for coordinate in range(3):
        for label in range(252):
            indices=np.flatnonzero(eq_labels[:,coordinate]==label)
            gram=weighted[indices]@equator[indices].T
            np.fill_diagonal(gram,6)
            assert gram.max(initial=0)<=6
            products+=len(indices)*(len(indices)-1)//2
    # All other equator pairs have three different coset labels: their
    # lattice-distance lower bound is 4+4+4=12.
    lifts,lift_weighted,lift_labels=expand(c['lift_phase_orbits'],True)
    gram=lift_weighted@lifts.T
    np.fill_diagonal(gram,2);assert gram.max()<=2
    cross_products=0
    for index,cell in enumerate(lift_labels):
        nonzero=np.flatnonzero(cell!=252)
        matching=np.flatnonzero(np.all(eq_labels[:,nonzero]==cell[nonzero],axis=1))
        dots=equator[matching]@lift_weighted[index]
        assert dots.max(initial=0)<=6
        cross_products+=len(dots)
    # Other lift/equator pairs have at most one equal nonzero label, giving
    # body dot <=4+2=6. Different heights give total dot <=8-4=4.
    assert c['lift_heads']==[2,-2]
    n36=fixed_count+len(equator);n37=n36+2*len(lifts)
    assert (n36,n37)==(c['points_36'],c['points_37'])
    result=dict(points_36=n36,points_37=n37,base_points=fixed_count,
                expanded_phase_points=len(equator),expanded_lift_bodies=len(lifts),
                phase_coordinate_pair_products=products,
                lift_coordinate_pair_products=len(lifts)*(len(lifts)-1)//2,
                lift_equator_coordinate_products=cross_products,
                exact_shell_counts=dict(Counter(map(int,norms))),
                all_checks_pass=True,certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                arithmetic='Fraction enumeration and int64 coordinate matrix products',
                elapsed_seconds=time.monotonic()-started)
    (HERE/'coordinate_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
