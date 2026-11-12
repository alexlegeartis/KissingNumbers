#!/usr/bin/env python3
"""Independent exact proof in orthonormal D4 coordinates, standard library only.

The 72-dimensional lattice minimum, shell population and spherical 11-design
property are published premises. This verifies their application, the actual
normal embedding, finite orientation averaging, rational weak duality and
every projected pair type. It does not enumerate the ambient shell.
"""
from fractions import Fraction as F
from itertools import product, permutations, combinations_with_replacement
from pathlib import Path
from math import prod, gcd
import hashlib, json, sys

HERE=Path(__file__).resolve().parent
POPULATION=6218175600
EMBEDDING_SHA256='5f180d5d6c737010083e72210e16b92926ec2e79ed9a0891d0960a018b2d7842'

def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]

def transpose(a): return list(map(list,zip(*a)))

def hadamard(h):
    z=(h[0]+h[1]+h[2]+h[3],h[0]+h[1]-h[2]-h[3],
       h[0]-h[1]+h[2]-h[3],h[0]-h[1]-h[2]+h[3])
    assert all(x%2==0 for x in z)
    return tuple(x//2 for x in z)

def domain_and_groups():
    # Dual of 2D4: h=4p is all even or all odd. For any minimum v
    # outside the normal roots, all root products have absolute value <=4.
    cells=[]
    for parity in [0,1]:
        for h in product(range(-8+parity,9,2),repeat=4):
            if sum(sorted(map(abs,h))[-2:])<=8: cells.append(h)
    cells=set(cells)
    roots=set()
    for i in range(4):
        for j in range(i):
            for a,b in product([-8,8],repeat=2):
                h=[0]*4;h[i]=a;h[j]=b;roots.add(tuple(h))
    assert len(roots)==24 and len(cells)==1681
    assert {hadamard(h) for h in roots}==roots
    assert {hadamard(h) for h in cells}==cells
    # Signed coordinate permutations and this Hadamard isometry preserve
    # the root lattice. Their group is finite: it acts faithfully on the
    # finite spanning set of 24 roots. Thus averaging real normal tuples
    # is legitimate without any ambient lattice symmetry assumption.
    reps=sorted({tuple(sorted(map(abs,h))) for h in cells},key=lambda h:(sum(x*x for x in h),h))
    index={r:i for i,r in enumerate(reps)}
    parent=list(range(len(reps)))
    def root(j):
        while parent[j]!=j:
            parent[j]=parent[parent[j]];j=parent[j]
        return j
    for h in cells:
        a=root(index[tuple(sorted(map(abs,h)))])
        b=root(index[tuple(sorted(map(abs,hadamard(h))))])
        if a!=b:parent[max(a,b)]=min(a,b)
    identifiers=sorted({root(i) for i in range(len(reps))})
    groups=[sorted(h for h in cells if root(index[tuple(sorted(map(abs,h)))])==r) for r in identifiers]
    assert len(groups)==16
    return cells,roots,groups

def check(path):
    if not __debug__: raise RuntimeError('Run without -O: exact assertions are part of verification.')
    c=json.loads(Path(path).read_text())
    assert c['schema']=='gamma72-d4-head-compression-v1'
    assert c['population']==POPULATION and c['dimension']==69
    raw=(HERE/'embedding.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==EMBEDDING_SHA256
    e=json.loads(raw);ambient=e['catalogue_gram'];normals=e['normals'][:4]
    assert len(ambient)==72 and all(len(row)==72 for row in ambient)
    assert ambient==transpose(ambient)
    assert all(isinstance(x,int) for row in normals for x in row)
    gram=matmul(matmul(normals,ambient),transpose(normals))
    basis=[[1,-1,0,0],[0,1,-1,0],[0,1,1,0],[-1,0,0,1]]
    target=[[4*x for x in row] for row in matmul(basis,transpose(basis))]
    assert gram==target==[[8,-4,-4,-4],[-4,8,0,0],[-4,0,8,0],[-4,0,0,8]]
    # det(basis)=2. The normal span is an actual copy of 2D4 in Gamma72.
    cells,roots,groups=domain_and_groups()
    assert [[list(h) for h in gg] for gg in groups]==c['groups']
    profiles=list(map(tuple,c['profiles']))
    assert len(profiles)==10 and len(set(profiles))==len(profiles)
    assert all(h in cells for h in profiles)
    powers=c['moment_powers'];rows=[];rhs=[]
    for pp in powers:
        assert len(pp)==4 and all(isinstance(x,int) and x>=0 and x%2==0 for x in pp)
        degree=sum(pp)//2;assert 0<=degree<=5
        value=lambda h:prod(x**p for x,p in zip(h,pp))
        expectation=F(POPULATION*128**degree*prod(prod(range(1,p,2)) for p in pp),
                      prod(72+2*j for j in range(degree)))
        b=expectation-sum(value(h) for h in roots)
        row=[sum(value(h) for h in gg)*b.denominator for gg in groups]
        bn=b.numerator;g=gcd(*row,bn)
        rows.append([x//g for x in row]);rhs.append(bn//g)
    assert rows==c['rows'] and rhs==c['rhs']
    objective=[sum(h in profiles for h in gg) for gg in groups]
    assert objective==c['objective']
    denominator=c['dual_denominator'];multipliers=c['dual_multipliers']
    assert isinstance(denominator,int) and denominator>0
    assert len(multipliers)==len(rows)==10 and all(isinstance(x,int) for x in multipliers)
    slack=[denominator*objective[j]-sum(a*row[j] for a,row in zip(multipliers,rows)) for j in range(16)]
    assert min(slack)>=0
    numerator=sum(a*b for a,b in zip(multipliers,rhs))
    exact=F(numerator,denominator)
    bound=-((-numerator)//denominator)
    assert bound==c['bound']==743522304
    # The objective is the mean population over the finite orientations.
    # At least one actual orientation has integer population >=ceil(exact).
    heights=list(map(F,c['heights']));assert len(heights)==len(profiles)
    inner=lambda a,b:F(sum(x*y for x,y in zip(a,b)),16)
    q=[inner(h,h) for h in profiles]
    alpha2=[(8-t*t)/(8-r) for t,r in zip(heights,q)]
    assert all(0<=r<8 for r in q) and all(a>0 for a in alpha2)
    assert all(t*t<=r for t,r in zip(heights,q))
    assert all(a*(8-r)+t*t==8 for a,r,t in zip(alpha2,q,heights))
    gaps=[];tight=[]
    for i,j in combinations_with_replacement(range(len(profiles)),2):
        d=4-inner(profiles[i],profiles[j]);r=4-heights[i]*heights[j]
        assert d>=0 and r>=0
        gap=r*r-alpha2[i]*alpha2[j]*d*d
        assert gap>=0,(i,j,gap)
        gaps.append(gap)
        if gap==0:tight.append([i,j])
    assert len(gaps)==55
    return dict(status='PASS',dimension=69,bound=bound,exact_dual_value=str(exact),
                prior_bound=c['prior_bound'],gain=bound-c['prior_bound'],
                complete_profile_domain=len(cells),normal_roots=len(roots),
                averaging_orbits=len(groups),independent_moment_rows=len(rows),
                selected_slices=len(profiles),exact_geometric_pair_types=len(gaps),
                tight_pairs=tight,minimum_positive_squared_margin=str(min(x for x in gaps if x)),
                dual_slack_min=min(slack),actual_embedding_verified=True,
                verification_scope='Exact existence proof using published Gamma72 minimum, shell population and 11-design premises; no full shell enumeration.')

if __name__=='__main__':
    print(json.dumps(check(sys.argv[1] if len(sys.argv)>1 else HERE/'orthonormal_certificate.json'),indent=2))
