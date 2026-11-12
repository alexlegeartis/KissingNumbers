#!/usr/bin/env python3
"""Exact nested-normal population and geometry checker; standard library only.

Published premises: Gamma72 has minimum squared norm 8, 6,218,175,600
minimal vectors, and its minimal shell is a spherical 11-design. This checks
an actual integral normal embedding, the complete possible profile domain,
finite orientation averaging, rational weak duality, and the geometry.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement, product
from math import gcd, prod
from pathlib import Path
import json
import sys

import exact_orbits as model

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'embedding.json'
SOURCE_SHA256='5f180d5d6c737010083e72210e16b92926ec2e79ed9a0891d0960a018b2d7842'

def transpose(a): return list(map(list,zip(*a)))
def multiply(a,b): return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]

def inverse(a):
    n=len(a)
    m=[[F(x) for x in r]+[F(i==j) for j in range(n)] for i,r in enumerate(a)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if m[i][j])
        m[j],m[pivot]=m[pivot],m[j]
        t=m[j][j];m[j]=[x/t for x in m[j]]
        for i in range(n):
            if i!=j:
                t=m[i][j];m[i]=[x-t*y for x,y in zip(m[i],m[j])]
    return [r[n:] for r in m]

def embedding(k,data_dir=HERE):
    raw=(data_dir/'embedding.json').read_bytes()
    assert sha256(raw).hexdigest()==SOURCE_SHA256
    e=json.loads(raw)
    assert k==7
    witness=json.loads((data_dir/'d7_embedding.json').read_text())
    normals=witness['normals']
    assert len(normals)==k and all(len(r)==72 for r in normals)
    assert all(type(x) is int for r in normals for x in r)
    gram=multiply(multiply(normals,e['catalogue_gram']),transpose(normals))
    basis=[[1,-1,0,0]+[0]*(k-4),[0,1,-1,0]+[0]*(k-4),
           [0,1,1,0]+[0]*(k-4),[-1,0,0,1]+[0]*(k-4)]
    for i in range(4,k):
        row=[0]*k;row[i-1]=-1;row[i]=1;basis.append(row)
    target=[[4*x for x in r] for r in multiply(basis,transpose(basis))]
    assert gram==target==witness['gram']
    # Rows of this integer matrix generate D_k: start from the D4 root
    # basis of determinant 2 and append one coordinate-difference root.
    bi=inverse(basis)
    dual=[[2*x for x in r] for r in bi]
    assert all(x.denominator==1 for r in dual for x in r)
    assert all(len({int(x)%2 for x in c})==1 for c in transpose(dual))
    return normals,gram

def check(path,enumerate_profiles=False,data_dir=HERE):
    if not __debug__: raise RuntimeError('Run without -O; assertions are part of the checker.')
    c=json.loads(Path(path).read_text())
    assert c['schema']=='gamma72-d7-extended-moments-v1'
    k=c['k'];dimension=c['dimension']
    assert k==7 and dimension in (68,69)
    assert F(c['height_scale_squared'])==1
    assert c['population']==model.N
    embedding(k,data_dir)
    profiles=[tuple(h) for h in c['profiles']]
    assert all(len(h)==4 and all(type(x) is int for x in h) for h in profiles)
    assert len(profiles)==len(set(profiles))
    marginal=model.d4_marginal_weights(profiles)
    reps=model.representatives(k)
    sizes=[model.orbit_size(h) for h in reps]
    assert c['representatives']==[list(h) for h in reps]
    assert c['orbit_sizes']==sizes
    roots=list(model.normal_roots(k))
    assert len(roots)==2*k*(k-1)
    # h=4p is integral of common parity by normal duality. Every minimal
    # vector except these roots satisfies |h_i|+|h_j|<=8, and |h|^2<=128.
    assert all(sum(hh*hh for hh in h)<=128 for h in reps)
    obj=[model.objective(h,marginal) for h in reps]
    known=sum((marginal.get(tuple(sorted(map(abs,h[:4]))),F(0)) for h in roots),F(0))
    assert known==F(c['known_objective'])==0
    assert obj==list(map(F,c['objective']))
    rows=[];rhs=[]
    for pp in c['tags']:
        assert len(pp)==k and all(type(x) is int and x>=0 and x%2==0 for x in pp)
        degree=sum(pp)//2;assert 0<=degree<=5
        row=[model.orbit_moment(h,tuple(pp)) for h in reps]
        b=F(model.N*128**degree*prod(prod(range(1,p,2)) for p in pp),
            prod(72+2*j for j in range(degree)))
        b-=sum(prod(x**p for x,p in zip(h,pp)) for h in roots)
        row=[x*b.denominator for x in row]
        bn=b.numerator;g=gcd(*row,bn)
        rows.append([x//g for x in row]);rhs.append(bn//g)
    assert rows==c['rows'] and rhs==c['rhs']
    if enumerate_profiles:
        # Independent direct sums check the compact combinatorial formulas.
        for j,rep in enumerate(reps):
            orbit=list(model.orbit(rep))
            assert len(orbit)==len(set(orbit))==sizes[j]
            assert sum((marginal.get(tuple(sorted(map(abs,h[:4]))),F(0)) for h in orbit),F(0))==obj[j]
            for pp in c['tags']:
                active=[(i,p) for i,p in enumerate(pp) if p]
                direct=sum(prod(h[i]**p for i,p in active) for h in orbit)
                assert direct==model.orbit_moment(rep,tuple(pp))
    y=list(map(F,c['dual_exact']))
    assert len(y)==len(rows)
    slack=[obj[j]-sum(a*r[j] for a,r in zip(y,rows)) for j in range(len(reps))]
    assert all(x>=0 for x in slack)
    assert slack==list(map(F,c['slack_exact']))
    value=sum(a*b for a,b in zip(y,rhs))+known
    assert value==F(c['exact_value'])
    ceiling=-((-value.numerator)//value.denominator)
    assert ceiling==c['exact_bound']
    # First average the original selection over D4 orientations inside
    # the first four normal coordinates, then average D_k orientations.
    # Its exact expected population is bounded above by some actual
    # integer population. Extra k-4 normal coordinates are RETAINED in
    # the 68-dimensional body; no k-dimensional projection is performed.
    heights=list(map(F,c['heights']))
    assert len(heights)==len(profiles)
    inner=lambda a,b:F(sum(x*y for x,y in zip(a,b)),16)
    q=[inner(h,h) for h in profiles]
    assert all(0<=x<8 for x in q)
    alpha2=[(8-t*t)/(8-x) for t,x in zip(heights,q)]
    assert all(x>0 for x in alpha2)
    assert all(a*(8-x)+t*t==8 for a,x,t in zip(alpha2,q,heights))
    if dimension==68: assert all(t==0 for t in heights)
    gaps=[]
    for i,j in combinations_with_replacement(range(len(profiles)),2):
        d=4-inner(profiles[i],profiles[j]);r=4-heights[i]*heights[j]
        assert d>=0 and r>=0
        gap=r*r-alpha2[i]*alpha2[j]*d*d
        assert gap>=0,(i,j,gap)
        gaps.append(gap)
    poles=0
    if dimension==69:
        assert all(t*t<=2 for t in heights)
        # Pure head points (0,+/-sqrt(8)) have dot product <=4 with
        # every image, and each image has a nonzero body component.
        poles=2
        assert c['with_two_poles']==ceiling+poles
    assert c['bound']==ceiling+poles
    assert c['prior_bound']=={68:471884288,69:743522304}[dimension]
    return dict(status='PASS',dimension=dimension,normal_rank=k,
                prior_bound=c['prior_bound'],gain=c['bound']-c['prior_bound'],
                exact_mean_lower_bound=str(value),population_ceiling=ceiling,
                added_poles=poles,bound=ceiling+poles,
                complete_possible_profiles=sum(sizes),known_normal_roots=len(roots),
                averaging_orbits=len(reps),moment_rows=len(rows),
                exact_geometric_pair_types=len(gaps),
                minimum_dual_slack=str(min(slack)),
                direct_enumeration_cross_check=enumerate_profiles,
                published_inputs='Gamma72 integrality, minimum 8, shell population, and spherical 11-design.')

if __name__=='__main__':
    args=[a for a in sys.argv[1:] if a!='--enumerate']
    paths=[Path(a) for a in args] if args else [HERE/'certificate68.json',HERE/'certificate69.json']
    results=[check(path,'--enumerate' in sys.argv) for path in paths]
    print(json.dumps(results,indent=2))
