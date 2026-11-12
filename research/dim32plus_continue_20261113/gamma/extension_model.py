"""Exact compact orbit moments for nested D_k -> D4 population averaging.

Model building uses the standard library. The optional search command uses
NumPy, SciPy and SymPy; its results are not certificates until checked exactly.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial, gcd, prod
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
N = 6218175600

def representatives(k):
    out = []
    for parity in (0,1):
        for h in combinations_with_replacement(range(parity,9,2),k):
            if h[-2]+h[-1] <= 8 and sum(x*x for x in h) <= 128:
                out.append(h)
    return sorted(out,key=lambda h:(sum(x*x for x in h),h))

def orbit_size(h):
    return factorial(len(h))*2**sum(x!=0 for x in h)//prod(factorial(t) for t in Counter(h).values())

def multiset_permutations(h):
    counts = Counter(h)
    def rec(prefix):
        if len(prefix)==len(h):
            yield tuple(prefix)
            return
        for x in sorted(counts):
            if counts[x]:
                counts[x]-=1
                yield from rec(prefix+[x])
                counts[x]+=1
    yield from rec([])

def orbit(h):
    for p in multiset_permutations(h):
        nonzero=[i for i,x in enumerate(p) if x]
        for signs in product((-1,1),repeat=len(nonzero)):
            out=list(p)
            for i,s in zip(nonzero,signs): out[i]*=s
            yield tuple(out)

def hadamard(h):
    z=(h[0]+h[1]+h[2]+h[3],h[0]+h[1]-h[2]-h[3],
       h[0]-h[1]+h[2]-h[3],h[0]-h[1]-h[2]+h[3])
    assert all(x%2==0 for x in z)
    return tuple(x//2 for x in z)

def d4_marginal_weights(profiles):
    """Average the indicator of selected profiles under full D4 orientations."""
    reps=representatives(4)
    canonical=lambda h:tuple(sorted(map(abs,h)))
    parent={h:h for h in reps}
    def root(h):
        if parent[h]!=h: parent[h]=root(parent[h])
        return parent[h]
    for r in reps:
        for h in orbit(r):
            a,b=root(r),root(canonical(hadamard(h)))
            if a!=b: parent[max(a,b)]=min(a,b)
    sizes=Counter()
    for r in reps: sizes[root(r)]+=orbit_size(r)
    counts=Counter(root(canonical(h)) for h in profiles)
    return {r:F(counts[root(r)],sizes[root(r)]) for r in reps}

def submultisets(h, size):
    counts=tuple(sorted(Counter(h).items()))
    def rec(i,left,chosen,weight):
        if i==len(counts):
            if left==0: yield tuple(chosen),weight
            return
        value,total=counts[i]
        for take in range(min(left,total)+1):
            yield from rec(i+1,left-take,chosen+[value]*take,weight*comb(total,take))
    yield from rec(0,size,[],1)

def objective(h, marginal):
    size=len(next(iter(marginal)))
    return orbit_size(h)*sum((F(weight,comb(len(h),size))*marginal.get(sub,F(0))
                             for sub,weight in submultisets(h,size)),F(0))

def a2_marginal_weights(offsets):
    selected=set(map(tuple,offsets))
    result={}
    reps=set(representatives(3))|{tuple(sorted(map(abs,h))) for h in normal_roots(3)}
    for rep in reps:
        hits=sum(((h[0]-h[1])//2,(h[0]-h[2])//2) in selected for h in orbit(rep))
        result[rep]=F(hits,orbit_size(rep))
    return result

@lru_cache(None)
def orbit_moment(h,powers):
    """Sum a fixed even monomial over a signed coordinate-permutation orbit."""
    assert all(p%2==0 and p>=0 for p in powers)
    active=tuple(p for p in powers if p)
    counts=tuple(sorted(Counter(h).items()))
    @lru_cache(None)
    def rec(i,remaining):
        if i==len(active): return 1
        out=0
        for j,((value,total),available) in enumerate(zip(counts,remaining)):
            if available:
                rest=list(remaining);rest[j]-=1
                out+=available*value**active[i]*rec(i+1,tuple(rest))
        return out
    numer=orbit_size(h)*rec(0,tuple(t for v,t in counts))
    denom=prod(range(len(h)-len(active)+1,len(h)+1))
    assert numer%denom==0
    return numer//denom

def partitions(m,k,lo=0):
    if k==0:
        if m==0: yield ()
        return
    for x in range(lo,m+1):
        for rest in partitions(m-x,k-1,x): yield (x,)+rest

def normal_roots(k):
    for i,j in combinations(range(k),2):
        for a,b in product((-8,8),repeat=2):
            h=[0]*k;h[i]=a;h[j]=b
            yield tuple(h)

def model(k, profiles, marginal=None):
    reps=representatives(k)
    if marginal is None: marginal=d4_marginal_weights(profiles)
    marginal_size=len(next(iter(marginal)))
    obj=[objective(h,marginal) for h in reps]
    roots=list(normal_roots(k))
    # These roots are handled separately from the unknown profile domain.
    known_objective=sum((marginal.get(tuple(sorted(map(abs,h[:marginal_size]))),F(0)) for h in roots),F(0))
    rows=[];rhs=[];tags=[]
    for degree in range(6):
        for part in partitions(degree,k):
            powers=tuple(2*x for x in part)
            row=[orbit_moment(h,powers) for h in reps]
            b=F(N*128**degree*prod(prod(range(1,p,2)) for p in powers),
                prod(72+2*j for j in range(degree)))
            b-=sum(prod(x**p for x,p in zip(h,powers)) for h in roots)
            row=[x*b.denominator for x in row]
            bn=b.numerator;g=gcd(*row,bn)
            rows.append([x//g for x in row]);rhs.append(bn//g);tags.append(list(powers))
    return dict(k=k,population=N,representatives=[list(h) for h in reps],
                orbit_sizes=[orbit_size(h) for h in reps],rows=rows,rhs=rhs,tags=tags,
                objective=[str(x) for x in obj],known_objective=str(known_objective))

def selected_profiles(dim):
    if dim==69:
        c=json.loads((REPO/'verifications/improved/dim69-gamma72-head-compression/orthonormal_certificate.json').read_text())
        return c['profiles'],c['heights']
    if dim==70:
        c=json.loads((REPO/'verifications/improved/dim70-71-gamma72-cross-sections/affine70/certificate.json').read_text())
        return c['offsets'],['0']*len(c['offsets'])
    assert dim==68
    c=json.loads((REPO/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/certificate4.json').read_text())
    profiles=[]
    for a,b,c,d in c['offsets']:
        profiles.append([2*a+b+c,b+c,c-b,2*d+2*a+b+c])
    return profiles,['0']*len(profiles)

def solve(k,dim,profiles=None,heights=None,height_scale_squared=F(1),label=None):
    import numpy as np
    import sympy as sp
    from scipy.optimize import linprog
    if profiles is None: profiles,heights=selected_profiles(dim)
    assert heights is not None
    M=model(k,profiles,a2_marginal_weights(profiles) if dim>=70 else None)
    A=sp.Matrix(M['rows']);b=sp.Matrix(M['rhs'])
    choose=list(A.T.rref()[1])
    A=A[choose,:];b=b[choose,:]
    M['rows']=[M['rows'][i] for i in choose]
    M['rhs']=[M['rhs'][i] for i in choose]
    M['tags']=[M['tags'][i] for i in choose]
    sizes=np.asarray(M['orbit_sizes'],dtype=float)
    An=np.asarray(A,dtype=float)/sizes[None,:]
    bn=np.asarray(b,dtype=float).ravel()/N
    scales=np.maximum(np.abs(bn),1e-15)
    cn=np.asarray([float(F(x)) for x in M['objective']])/sizes
    res=linprog(cn,A_eq=An/scales[:,None],b_eq=bn/scales,bounds=(0,None),
                method='highs',options=dict(threads=1,primal_feasibility_tolerance=1e-9,
                                           dual_feasibility_tolerance=1e-9,presolve=False))
    assert res.success,res.message
    active=[i for i,x in enumerate(res.x) if x>1e-11]
    c=sp.Matrix([sp.Rational(x) for x in M['objective']])
    sols=list(sp.linsolve((A[:,active].T,c[active,:])))
    exact=None
    if sols:
        y=list(sols[0]);free=set().union(*(x.free_symbols for x in y))
        if free:
            y=[x.subs({f:0 for f in free}) for x in y]
        y=sp.Matrix(y);slack=c-A.T*y
        if all(x>=0 for x in slack):
            exact=(b.T*y)[0]+sp.Rational(M['known_objective'])
            M.update(dual_exact=[str(x) for x in y],slack_exact=[str(x) for x in slack],
                     exact_value=str(exact),exact_bound=int(sp.ceiling(exact)))
    M.update(dimension=dim,profiles=profiles,heights=heights,value=float(res.fun*N+F(M['known_objective'])),
             height_scale_squared=str(height_scale_squared),
             primal=(res.x/sizes).tolist(),active_columns=active,
             scaled_residual=float(max(abs((An@res.x-bn)/scales))))
    if exact is not None and dim in (69,71) and all(F(t)**2*height_scale_squared<=2 for t in heights):
        M['with_two_poles']=M['exact_bound']+2
    out=HERE/(label or f'd{k}_compact{dim}.json')
    out.write_text(json.dumps(M,indent=2)+'\n')
    print(out.name,'groups',len(sizes),'rank',len(choose),'float',M['value'],
          'exact',exact,'ceiling',M.get('exact_bound'),'poles',M.get('with_two_poles'),
          'residual',M['scaled_residual'],flush=True)
    return M

if __name__=='__main__':
    import sys
    for k in map(int,sys.argv[1:] or [5,6,7]):
        for dim in (68,69,70): solve(k,dim)
