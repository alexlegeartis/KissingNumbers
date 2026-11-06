#!/usr/bin/env python3
"""Replay exact 11-design identities for each E8 norm0/2/4/6 fibre.

Published mathematical inputs are the extremal48 minimum6, shell size
52,416,000 and11-design property. No ambient Weyl automorphisms are assumed:
the norm6 polynomial is averaged over an abstract stabilizer of its target.
"""
import json,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from flint import fmpq_mat
from old_e8_first_shells import e8

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'reports'
N=52416000
def poch(x,n):
    y=F(1)
    for i in range(n):y*=x+i
    return y
def expectation(a,b,c,v2,w2):
    m=a+b+c
    return F(72**m)*poch(F(4+b+c),a)*poch(F(1,2),b)*poch(F(1,2),c)*v2**b*w2**c/poch(F(24),m)
def rhs(label,v,w,roots):
    a,b,c=label;v2=int(v@v);w2=int(w@w)
    pure=sum(72**a*int(3*r@v)**(2*b)*int(3*r@w)**(2*c)for r in roots)
    return N*expectation(a,b,c,v2,w2)-pure
def main():
    start=time.monotonic();us=e8();q=np.sum(us*us,axis=1);roots=us[q==8]
    double={tuple(2*r)for r in roots}
    domain=np.array([u for u in us if tuple(u)not in double],dtype=np.int64)
    q=np.sum(domain*domain,axis=1)
    assert len(domain)==26401
    # Exhaust the normal projection bounds for every non-pure ambient minimal
    # vector. Pair roots bound the largest2 magnitudes by6; half roots and
    # squared norm<=72 imply l1<=18. No shell counts enter this enumeration.
    exhaustive=set()
    for vals in (range(-6,7,2),range(-5,6,2)):
        def visit(prefix,square,l1,largest,second):
            if len(prefix)==8:
                if sum(prefix)%4:return
                maximum=l1
                if all(prefix)and sum(x<0 for x in prefix)%2:maximum-=2*min(abs(x)for x in prefix)
                if maximum<=12:exhaustive.add(tuple(prefix))
                return
            for x in vals:
                a=abs(x);top,nxt=(a,largest)if a>=largest else(largest,max(second,a))
                if square+x*x<=72 and l1+a<=18 and top+nxt<=6:
                    visit(prefix+[x],square+x*x,l1+a,top,nxt)
        visit([],0,0,0,0)
    assert exhaustive==set(map(tuple,domain))
    # Exact pair bounds against embedded E8 roots justify this full projection domain.
    assert np.max(np.abs(domain@roots.T))==12
    reports=[]
    for v,expected in [(np.zeros(8,dtype=np.int64),1092000),
                       (np.array([2,2,0,0,0,0,0,0]),116235)]:
        v2=int(v@v);t=(domain@v)**2
        keys=sorted(set(zip(map(int,q),map(int,t))))
        labels=[(a,b,0)for a in range(6)for b in range(6-a)] if v2 else [(a,0,0)for a in range(5)]
        mat=[[r**a*s**b for a,b,_ in labels]+[int(r==v2 and s==v2*v2)]for r,s in keys]
        rr,rank=fmpq_mat(mat).rref();coef=[F(0)]*len(labels)
        for i in range(rank):
            pivot=next(j for j in range(rr.ncols())if rr[i,j])
            assert pivot<len(labels),('No degree10 identity',v.tolist(),len(keys),rank)
            coef[pivot]=F(str(rr[i,rr.ncols()-1]))
        assert all(sum(c*x for c,x in zip(coef,row[:-1]))==row[-1]for row in mat)
        value=sum(c*rhs(l,v,np.zeros(8,dtype=np.int64),roots)for c,l in zip(coef,labels))
        divisor=1 if not v2 else 2
        assert value==divisor*expected,(v.tolist(),value)
        reports.append({'lambda_norm':v2//4,'fibre_count':expected,
            'target':v.tolist(),'certificate':[(list(l),str(c))for l,c in zip(labels,coef)if c],
            'target_pair_moment_value':str(value)})
    # Norm4 target2e1. Its D7 stabilizer permutes the final7 coordinates and
    # flips an even number of their signs. Include its degree8 parity invariant.
    def d7key(u):
        sg=int(np.prod(np.sign(u[1:])))
        rest=tuple(sorted(abs(int(x))for x in u[1:]))
        return min((int(u[0]),rest,sg),(-int(u[0]),rest,-sg))
    groups={}
    for u in domain:groups.setdefault(d7key(u),[]).append(u)
    groups=list(groups.values());v=np.array([4,0,0,0,0,0,0,0]);target=[]
    for group in groups:target.append(sum(np.array_equal(u,v)or np.array_equal(u,-v)for u in group))
    def partitions(n,k=7,largest=5):
        if k==0:
            if n==0:yield ()
            return
        for x in range(min(n,largest),-1,-1):
            for tail in partitions(n-x,k-1,x):yield (x,)+tail
    exponents=[]
    for a0 in range(6):
        for n in range(6-a0):
            for tail in partitions(n):exponents.append((2*a0,)+tuple(2*x for x in tail))
    exponents.extend([(1,)*8,(3,)+(1,)*7,(1,3)+(1,)*6])
    rows=[]
    for alpha in exponents:
        rows.append([sum(int(np.prod(u**np.array(alpha,dtype=np.int64)))for u in group)for group in groups])
    mat=[[rows[j][i]for j in range(len(rows))]+[target[i]]for i in range(len(groups))]
    rr,rank=fmpq_mat(mat).rref();coef=[F(0)]*len(rows)
    for i in range(rank):
        pivot=next(j for j in range(rr.ncols())if rr[i,j]);assert pivot<len(rows)
        coef[pivot]=F(str(rr[i,rr.ncols()-1]))
    assert all(sum(c*x for c,x in zip(coef,row[:-1]))==row[-1]for row in mat)
    def mono_rhs(alpha):
        degree=sum(alpha)
        if any(a%2 for a in alpha):moment=F(0)
        else:
            moment=F(N*72**(degree//2),1)/poch(F(24),degree//2)
            for a in alpha:moment*=poch(F(1,2),a//2)
        pure=sum(int(np.prod((3*r)**np.array(alpha,dtype=np.int64)))for r in roots)
        return moment-pure
    value=sum(c*mono_rhs(a)for c,a in zip(coef,exponents))
    assert value==18360,value
    reports.append({'lambda_norm':4,'fibre_count':9180,'target_pair_moment_value':str(value),
      'stabilizer_orbits':len(groups),'certificate':[(list(a),str(c))for a,c in zip(exponents,coef)if c]})
    cert=json.loads((ROOT/'data/norm6_simple_certificate.json').read_text())
    orb=np.load(ROOT/'data/norm6_orbits.npz');X=orb['points'].astype(np.int64);ends=orb['ends']
    assert set(map(tuple,X))==set(map(tuple,domain)) and len(X)==len(domain)
    labels_by_point=np.repeat(np.arange(len(ends)-1),np.diff(ends))
    keys=np.ascontiguousarray(X).view('V64').ravel();order=np.argsort(keys);sortedkeys=keys[order]
    v=np.array([4,2,2,0,0,0,0,0]);w=np.array([0,2,-2,0,0,0,0,0]);assert v@w==0
    parents=list(range(len(X)))
    def find(i):
        while parents[i]!=i:
            parents[i]=parents[parents[i]];i=parents[i]
        return i
    def unite(i,j):
        i=find(i);j=find(j)
        if i!=j:parents[max(i,j)]=min(i,j)
    for r in roots[roots@v==0]:
        dots=X@r;assert np.all(dots%4==0)
        images=np.ascontiguousarray(X-(dots//4)[:,None]*r).view('V64').ravel()
        loc=np.searchsorted(sortedkeys,images);assert np.all(loc<len(X))
        assert np.all(sortedkeys[loc]==images)
        assert np.array_equal(labels_by_point[order[loc]],labels_by_point)
        for i,j in enumerate(order[loc]):unite(i,int(j))
    images=np.ascontiguousarray(-X).view('V64').ravel();loc=np.searchsorted(sortedkeys,images)
    assert np.all(sortedkeys[loc]==images)
    assert np.array_equal(labels_by_point[order[loc]],labels_by_point)
    for i,j in enumerate(order[loc]):unite(i,int(j))
    assert len({find(i)for i in range(len(X))})==len(ends)-1
    # Verify the abstract Weyl transitivity used to move each target vector.
    simple=np.load(ROOT/'data/prep_e8.npz')['A']
    for norm in (8,16,24):
        family={tuple(u)for u in domain[np.sum(domain*domain,axis=1)==norm]}
        seen={next(iter(family))};todo=list(seen)
        while todo:
            u=todo.pop()
            for root in simple:
                dot=sum(int(x)*int(y)for x,y in zip(u,root));assert dot%4==0
                image=tuple(int(x)-(dot//4)*int(y)for x,y in zip(u,root))
                assert image in family
                if image not in seen:seen.add(image);todo.append(image)
        assert seen==family
    q=np.sum(X*X,axis=1);t=(X@v)**2;s=(X@w)**2
    target=np.add.reduceat(np.all(X==v,axis=1).astype(np.int64)+np.all(X==-v,axis=1).astype(np.int64),ends[:-1])
    assert target.tolist()==cert['target']
    rows=[]
    for a,b,c in cert['labels']:
        rows.append(np.add.reduceat(q**a*t**b*s**c,ends[:-1]).tolist())
    assert rows==cert['rows']
    dual=[F(s)for s in cert['dual']]
    assert all(sum(dual[i]*rows[i][j]for i in range(len(rows)))==target[j]for j in range(len(target)))
    value=sum(c*rhs(l,v,w,roots)for c,l in zip(dual,cert['labels']))
    assert value==990,value
    reports.append({'lambda_norm':6,'fibre_count':495,'target_pair_moment_value':str(value),
         'stabilizer_orbits':len(target),'certificate_source':'record/norm6_simple_certificate.json'})
    result={'status':'PASS','uniform_fibres':reports,
       'projection_domain_count':26401,'pure_E8_roots':240,
       'published_inputs':['ambient minimum squared norm6','52,416,000 minimal vectors','minimal shell is an11design'],
       'all_shell_transitivity_scope':'Abstract W(E8) is transitive on each norm0,2,4,6 projection shell; no ambient automorphism assumption.',
       'total_count_with_explicit_norm8_heads':1092000+240*116235+2160*9180+6720*495+9*30240+240,
       'seconds':time.monotonic()-start}
    assert result['total_count_with_explicit_norm8_heads']==52416000
    (OUT/'e8_uniform_mass_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
