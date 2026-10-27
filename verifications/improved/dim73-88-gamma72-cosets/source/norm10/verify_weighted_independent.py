"""Independent reconstruction and integer-threshold audit of the new bounds."""
from pathlib import Path
from math import isqrt, prod
from itertools import permutations, combinations
from fractions import Fraction as F
from collections import Counter
import json
import numpy as np

root=Path(__file__).resolve().parent
if not __debug__:
    raise SystemExit('Run without -O: certificate assertions must be enabled.')
report=json.loads((root/'weighted_norm10_certificate.json').read_text())
old=np.load(root/'norm10_tails.npz', allow_pickle=False)

def sphere(n,q):
    if n==0:
        if q==0: yield ()
        return
    for x in range(-isqrt(q),isqrt(q)+1):
        for rest in sphere(n-1,q-x*x): yield (x,)+rest

def shell(k,q):
    if k==5:
        return np.array([tuple(2*x for x in v)+(0,)*3 for v in sphere(5,q//4)
                         if sum(v)%2==0],dtype=np.int64)
    m=9-k
    rows=[]
    for x in range(-isqrt(q//m),isqrt(q//m)+1):
        for rest in sphere(8-m,q-m*x*x):
            row=(x,)*m+rest
            if all(v%2==x%2 for v in row) and sum(row)%4==0: rows.append(row)
    return np.array(rows,dtype=np.int64)

# The reported distributions must satisfy every design identity. There are
# six remaining unknown populations and six independent moment equations.
small={int(j):v for j,v in report['minimal_shell_absolute_products_against_norm10'].items()}
heads={int(j):v for j,v in report['norm10_absolute_product_distribution'].items()}
for dist,pop,rt in [(small,6218175600,80),(heads,15281788354560,100)]:
    for t in range(6):
        assert sum(v*j**(2*t) for j,v in dist.items())==F(pop*prod(range(1,2*t,2))*rt**t,prod(range(72,72+2*t,2)))
assert heads[6]==small[4] and heads[10]==2

products=0
for k in [5,6,7]:
    row=report['dimensions'][str(72+k)]
    qs=[32,40] if k==5 else [8,40] if k==6 else [32]
    zs=[];rs=[]
    roots={tuple(2*v) for v in old['poles%d'%k]}
    for q in qs:
        z=shell(k,q)
        if k==5:
            rates=[1050 if q==32 and np.count_nonzero(v)==2 else 818 if q==32 else
                   882 if np.max(np.abs(v))==6 else 1000 for v in z]
        elif k==6:rates=[1015 if q==8 else 1000]*len(z)
        else:rates=[1000 if tuple(v) in roots else 311 for v in z]
        if q==(40 if k in [5,6] else 32):
            assert {tuple(v) for v in z}=={tuple(v) for v in old['tails%d'%k]}
        zs.extend(z);rs.extend(rates)
    z=np.array(zs);rates=np.array(rs);norms=np.sum(z*z,axis=1)
    groups=sorted(set(zip(norms.tolist(),rates.tolist())))
    s=F(row['head_scale_squared']);p,d=s.numerator,s.denominator
    profile=Counter()
    for v,rate,q in zip(z,rates,norms):
        dots=z@v
        closed=int(rates.sum())
        for r,weight in groups:
            selected=(norms==r)&(rates==weight)
            for j in [4,5,6]:
                # dot is integral. dot>sqrt(C) is exactly dot>floor(sqrt(C));
                # this tests boundaries independently of the squared-mask implementation.
                cutoff=isqrt(((4*d-j*p)**2*int(q)*r)//(8*d-10*p)**2)
                closed+=(heads[j]//2)*weight*int(np.count_nonzero(selected&(dots>cutoff)))
        profile[(int(rate),closed)]+=1
    assert [[r,c,n] for (r,c),n in sorted(profile.items())]==row['weighted_closed_profiles']
    alpha=(15281788354560//2)*sum((F(r*n,c) for (r,c),n in profile.items()),F(0))
    assert alpha==F(int(row['priority_expectation_numerator']),int(row['priority_expectation_denominator']))
    pairs=-(-alpha.numerator//alpha.denominator)
    assert row['bound']==6218175600+2*pairs+len(old['poles%d'%k])
    products+=len(z)**2
    print('PASS independent tails, strict integer cutoffs, every neighborhood, exact expectation:',72+k,row['bound'],flush=True)
assert products==22848388

# Enumerate every graph on four vertices and every exponential-clock ordering.
# This checks the weighted local-minimum expectation independently by summing
# exact Plackett--Luce order probabilities, and tests the independence property.
edges=list(combinations(range(4),2));weights=[1,2,3,5]
for bits in range(1<<len(edges)):
    neighbors=[{i} for i in range(4)]
    for b,(i,j) in enumerate(edges):
        if bits>>b&1:neighbors[i].add(j);neighbors[j].add(i)
    expected=sum((F(weights[i],sum(weights[j] for j in neighbors[i])) for i in range(4)),F(0))
    enumeration=F(0);probabilities=F(0)
    for order in permutations(range(4)):
        prob=F(1);remaining=sum(weights)
        for i in order:prob*=F(weights[i],remaining);remaining-=weights[i]
        positions={v:i for i,v in enumerate(order)}
        kept={i for i in range(4) if all(positions[i]<=positions[j] for j in neighbors[i])}
        assert all(not(kept&neighbors[i]-{i}) for i in kept)
        enumeration+=len(kept)*prob;probabilities+=prob
    assert probabilities==1 and enumeration==expected
print('PASS 64 graphs x 24 clock orders; 22,848,388 tail products independently checked',flush=True)
